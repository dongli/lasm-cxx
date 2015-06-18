#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "QuadraturePoints.h"
#include "SkeletonPoints.h"
#include "Parcel.h"
#include "Tracers.h"

namespace lasm {

AdvectionManager::
AdvectionManager()
: geomtk::AdvectionManagerInterface<2, Domain, Mesh, Field, VelocityField>() {
    REPORT_ONLINE;
}

AdvectionManager::
~AdvectionManager() {
    REPORT_OFFLINE
}

void AdvectionManager::
init(const ConfigManager &configManager, const Mesh &mesh) {
    // Read in parameters.
    filamentLimit = configManager.getValue("lasm", "filament_limit", 5.0);
    minBiasLimit  = configManager.getValue("lasm", "min_bias_limit", 0.1);
    maxBiasLimit  = configManager.getValue("lasm", "max_bias_limit", 0.1);
    radialMixing  = configManager.getValue("lasm", "radial_mixing",  1.0);
    lateralMixing = configManager.getValue("lasm", "lateral_mixing", 1000.0);
    restoreFactor = configManager.getValue("lasm", "restore_factor", 0.001);
    reshapeFactor = configManager.getValue("lasm", "reshape_factor", 0.5);
    connectScale  = configManager.getValue("lasm", "connect_scale",  1.5);
    // Record domain and mesh.
    domain = &mesh.domain();
    this->mesh = &mesh;
    // Initialize objects.
    ShapeFunction::init(*domain);
    QuadraturePoints::init(*domain);
    SkeletonPoints::init(*domain);
    Parcel::init(mesh);
    Regrid::init(mesh);
    parcelManager.init(mesh);
    meshAdaptor.init(mesh);
    // Initialize tree objects.
    SpaceCoord x(domain->numDim());
#ifdef LASM_USE_RLL_MESH
    gridCoords.reshape(3, mesh.totalNumGridWithUniquePoleGrid(CENTER));
#else
    gridCoords.reshape(3, mesh.totalNumGrid(CENTER));
#endif
    int k = 0;
    for (uword i = 0; i < mesh.totalNumGrid(CENTER); ++i) {
#ifdef LASM_USE_RLL_MESH
        // In regular lat-lon mesh, there are whole zonal grids located on the
        // Poles. This will cause many computational problems, so we just use
        // the first grid to represent others.
        auto spanIdx = mesh.unwrapIndex(CENTER, i);
        if ((spanIdx[1] == mesh.js(FULL) && spanIdx[0] != mesh.is(FULL)) ||
            (spanIdx[1] == mesh.je(FULL) && spanIdx[0] != mesh.is(FULL))) {
            continue;
        }
#endif
        for (uword j = 0; j < 3; ++j) {
            gridCoords(j, k) = mesh.gridCoord(CENTER, i).cartCoord()[j];
        }
        k++;
    }
    MetricType::domain = domain;
    gridTree = new TreeType(gridCoords);
    // Connect parcels and grids Initially.
    TimeLevelIndex<2> timeIdx;
    connectParcelsAndGrids(timeIdx);
    refVolume = 0;
    for (auto parcel : parcelManager.parcels()) {
        refVolume += parcel->volume(timeIdx);
    }
    refVolume /= parcelManager.parcels().size();
} // init

void AdvectionManager::
addTracer(const string &name, const string &unit, const string &comment) {
    Tracers::add(name, comment, unit);
    for (auto parcel : parcelManager.parcels()) {
        parcel->tracers().add();
    }
    meshAdaptor.addTracer(name, comment, unit);
}

void AdvectionManager::
advance(double dt, const TimeLevelIndex<2> &newIdx, const VelocityField &velocityField) {
    integrate(dt, newIdx, velocityField);
    connectParcelsAndGrids(newIdx);
    remapDensityFromParcelsToGrids(newIdx);
    mixParcels(newIdx);
    statistics(newIdx);
} // advance

void AdvectionManager::
input(const TimeLevelIndex<2> &timeIdx, double *q) {
    // Copy input tracer density onto internal mesh grids.
    Tracers::resetTotalMasses();
    int l = 0;
    for (int t = 0; t < Tracers::numTracer(); ++t) {
        for (uword i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
            meshAdaptor.density(timeIdx, t, i) = q[l++];
            Tracers::totalMass(t) += meshAdaptor.mass(timeIdx, t, i);
        }
    }
    remapDensityFromGridsToParcels(timeIdx);
    remapDensityFromParcelsToGrids(timeIdx);
} // input

void AdvectionManager::
output(const TimeLevelIndex<2> &timeIdx, int ncId) const {
    // Append global attributes.
    nc_redef(ncId);
    nc_put_att(ncId, NC_GLOBAL, "filament_limit", NC_DOUBLE, 1, &filamentLimit);
    nc_put_att(ncId, NC_GLOBAL, "min_bias_limit", NC_DOUBLE, 1, &minBiasLimit);
    nc_put_att(ncId, NC_GLOBAL, "max_bias_limit", NC_DOUBLE, 1, &maxBiasLimit);
    nc_put_att(ncId, NC_GLOBAL, "radial_mixing",  NC_DOUBLE, 1, &radialMixing);
    nc_put_att(ncId, NC_GLOBAL, "lateral_mixing", NC_DOUBLE, 1, &lateralMixing);
    nc_put_att(ncId, NC_GLOBAL, "restore_factor", NC_DOUBLE, 1, &restoreFactor);
    nc_put_att(ncId, NC_GLOBAL, "reshape_factor", NC_DOUBLE, 1, &reshapeFactor);
    string str = to_iso_extended_string(boost::gregorian::day_clock::universal_day());
    nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, str.length(), str.c_str());
    nc_enddef(ncId);
    // Append tracer data.
    parcelManager.output(timeIdx, ncId);
} // output

void AdvectionManager::
statistics(const TimeLevelIndex<2> &timeIdx) {
    double minVolume = 1.0e33;
    double maxVolume = -1.0e33;
    double maxVolumeChange = -1.0e33;
    double maxFilament = -1.0e33;
    int maxFilamentParcelId;
    for (auto parcel : parcelManager.parcels()) {
        if (minVolume > parcel->volume(timeIdx)) minVolume = parcel->volume(timeIdx);
        if (maxVolume < parcel->volume(timeIdx)) maxVolume = parcel->volume(timeIdx);
        double dv = parcel->volume(timeIdx)-parcel->volume(timeIdx-1);
        if (maxVolumeChange < dv) maxVolumeChange = dv;
        if (maxFilament < parcel->filament()) {
            maxFilament = parcel->filament();
            maxFilamentParcelId = parcel->id();
        }
    }
    cout << setw(120) << setfill('-') << "-" << setfill(' ') << endl;
    cout << setw(30) << "minVolume";
    cout << setw(30) << "maxVolume";
    cout << setw(30) << "maxVolumeChange";
    cout << setw(30) << "maxFilament" << endl;
    cout << setw(30) << setprecision(10) << minVolume/refVolume;
    cout << setw(30) << setprecision(10) << maxVolume/refVolume;
    cout << setw(30) << setprecision(10) << maxVolumeChange/refVolume;
    cout << setw(30) << setprecision(10) << maxFilament << " (" << maxFilamentParcelId << ")" << endl;
} // statistics

void AdvectionManager::
integrate(double dt, const TimeLevelIndex<2> &newIdx,
          const VelocityField &velocityField) {
    TimeLevelIndex<2> oldIdx = newIdx-1;
    TimeLevelIndex<2> halfIdx = newIdx-0.5;
    double dt05 = 0.5*dt;
    const auto &divergence = velocityField.divergence();
    Tracers::resetTotalMasses();
    for (auto parcel : parcelManager.parcels()) {
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div, k1, k2, k3, k4, volume;
        SpaceCoord &x0 = parcel->x(oldIdx);
        SpaceCoord &x1 = parcel->x(newIdx);
        MeshIndex &I0 = parcel->meshIndex(oldIdx);
        MeshIndex &I1 = parcel->meshIndex(newIdx);
#ifdef LASM_IN_SPHERE
        // TODO: Should we hide the following codes? Because they are
        //       related to sphere domain and RLL mesh.
        if (I0.isOnPole()) {
            I0.setMoveOnPole(true);
            I1.setMoveOnPole(true);
            x0.transformToPS(*domain);
        } else {
            I0.setMoveOnPole(false);
            I1.setMoveOnPole(false);
        }
#endif
        // Stage 1.
        regrid.run(LINEAR, oldIdx, velocityField, x0, v1, &I0);
        regrid.run(LINEAR, oldIdx, divergence, x0, div, &I0);
        mesh->move(x0, dt05, v1, I0, x1);
        if (!domain->isValid(x1)) {
            REPORT_ERROR("Parcel moved out of the domain!");
        }
        I1.locate(*mesh, x1);
        k1 = parcel->volume(oldIdx)*div;
        volume = parcel->volume(oldIdx)+dt05*k1;
        // Stage 2.
        regrid.run(LINEAR, halfIdx, velocityField, x1, v2, &I1);
        regrid.run(LINEAR, halfIdx, divergence, x1, div, &I1);
        mesh->move(x0, dt05, v2, I0, x1);
        if (!domain->isValid(x1)) {
            REPORT_ERROR("Parcel moved out of the domain!");
        }
        I1.locate(*mesh, x1);
        k2 = volume*div;
        volume = parcel->volume(oldIdx)+dt05*k2;
        // Stage 3.
        regrid.run(LINEAR, halfIdx, velocityField, x1, v3, &I1);
        regrid.run(LINEAR, halfIdx, divergence, x1, div, &I1);
        mesh->move(x0, dt, v3, I0, x1);
        if (!domain->isValid(x1)) {
            REPORT_ERROR("Parcel moved out of the domain!");
        }
        I1.locate(*mesh, x1);
        k3 = volume*div;
        volume = parcel->volume(oldIdx)+dt*k3;
        // Stage 4.
        regrid.run(LINEAR, newIdx, velocityField, x1, v4, &I1);
        regrid.run(LINEAR, newIdx, divergence, x1, div, &I1);
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, I0, x1);
        if (!domain->isValid(x1)) {
            REPORT_ERROR("Parcel moved out of the domain!");
        }
        I1.locate(*mesh, x1);
        k4 = volume*div;
        parcel->updateVolume(newIdx, parcel->volume(oldIdx)+dt*(k1+2*k2+2*k3+k4)/6);
        // Update the skeleton points.
        field<SpaceCoord> &xs0 = parcel->skeletonPoints().spaceCoords(oldIdx);
        field<SpaceCoord> &xs1 = parcel->skeletonPoints().spaceCoords(newIdx);
        field<MeshIndex> &Is0 = parcel->skeletonPoints().meshIndexs(oldIdx);
        field<MeshIndex> &Is1 = parcel->skeletonPoints().meshIndexs(newIdx);
        for (uword i = 0; i < xs0.size(); ++i) {
            if (!Is0[i].isValid()) continue;
#ifdef LASM_IN_SPHERE
            if (Is0[i].isOnPole()) {
                Is0[i].setMoveOnPole(true);
                Is1[i].setMoveOnPole(true);
                xs0[i].transformToPS(*domain);
            } else {
                Is0[i].setMoveOnPole(false);
                Is1[i].setMoveOnPole(false);
            }
#endif
            // Stage 1.
            regrid.run(LINEAR, oldIdx, velocityField, xs0[i], v1, &Is0[i]);
            mesh->move(xs0[i], dt05, v1, Is0[i], xs1[i]);
            if (!domain->isValid(xs1[i])) {
                Is1[i].reset();
                continue;
            }
            Is1[i].locate(*mesh, xs1[i]);
            // Stage 2.
            regrid.run(LINEAR, halfIdx, velocityField, xs1[i], v2, &Is1[i]);
            mesh->move(xs0[i], dt05, v2, Is0[i], xs1[i]);
            if (!domain->isValid(xs1[i])) {
                Is1[i].reset();
                continue;
            }
            Is1[i].locate(*mesh, xs1[i]);
            // Stage 3.
            regrid.run(LINEAR, halfIdx, velocityField, xs1[i], v3, &Is1[i]);
            mesh->move(xs0[i], dt, v3, Is0[i], xs1[i]);
            if (!domain->isValid(xs1[i])) {
                Is1[i].reset();
                continue;
            }
            Is1[i].locate(*mesh, xs1[i]);
            // Stage 4.
            regrid.run(LINEAR, newIdx, velocityField, xs1[i], v4, &Is1[i]);
            v = (v1+v2*2.0+v3*2.0+v4)/6.0;
            mesh->move(xs0[i], dt, v, Is0[i], xs1[i]);
            if (!domain->isValid(xs1[i])) {
                Is1[i].reset();
                continue;
            }
            Is1[i].locate(*mesh, xs1[i]);
        }
        parcel->updateDeformMatrix(newIdx);
        for (uword t = 0; t < Tracers::numTracer(); ++t) {
            parcel->tracers().mass(t) += parcel->tracers().tendency(t)*dt;
            parcel->tracers().density(t) = parcel->tracers().mass(t)/parcel->volume(newIdx);
            Tracers::totalMass(t) += parcel->tracers().mass(t);
        }
    }
#ifndef NDEBUG
    if (domain->numDim() == 2) {
        double totalVolume = 0;
        for (auto parcel : parcelManager.parcels()) {
            totalVolume += parcel->volume(newIdx);
        }
#ifdef LASM_IN_SPHERE
        double trueTotalVolume = 4*PI*pow(domain->radius(), 2);
#elif defined LASM_IN_CARTESIAN
        double trueTotalVolume = totalVolume;
#endif
        double error = (totalVolume-trueTotalVolume)/trueTotalVolume;
        assert(fabs(error) < 1.0e-12);
    }
#endif
} // integrate

void AdvectionManager::
connectParcelAndGrids(const TimeLevelIndex<2> &timeIdx, Parcel *parcel) {
    meshAdaptor.containParcel(parcel->meshIndex(timeIdx).cellIndex(*mesh, CENTER), parcel);
    parcel->updateShapeSize(timeIdx);
    double longAxisSize = connectScale*parcel->shapeSize(timeIdx)[0];
    SearchType search(gridTree, NULL, gridCoords, parcel->x(timeIdx).cartCoord(), true);
    mlpack::math::Range r(0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    search.Search(r, neighbors, distances);
    BodyCoord y(domain->numDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int cellIdx = meshAdaptor.cellIndex(neighbors[0][i]);
        const SpaceCoord &x = meshAdaptor.coord(cellIdx);
        parcel->calcBodyCoord(timeIdx, x, y);
        // NOTE: Increase influence radius of parcel to avoid fluctuation.
        y() /= connectScale;
        double w = parcel->shapeFunction(timeIdx, y);
        if (w == 0) continue;
        meshAdaptor.connectParcel(cellIdx, parcel, w);
        parcel->connectCell(cellIdx);
    }
    meshAdaptor.connectParcel(parcel->hostCellIndex(), parcel, 1);
    parcel->connectCell(parcel->hostCellIndex());
} // connectParcelAndGrids

void AdvectionManager::
connectParcelsAndGrids(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetConnectedAndContainedParcels();
#ifdef LASM_USE_DIAG
    Diagnostics::resetMetric<Field<int> >("ncp1");
    Diagnostics::resetMetric<Field<int> >("ncp2");
#endif
    for (auto parcel : parcelManager.parcels()) {
        parcel->resetConnectedCells();
        connectParcelAndGrids(timeIdx, parcel);
    }
} // connectParcelsAndGrids

void AdvectionManager::
mixParcels(const TimeLevelIndex<2> &timeIdx) {
#ifdef LASM_USE_DIAG
    Diagnostics::resetMetric<Field<int> >("nmp");
#endif
    const int maxNumNeighborParcel = 1000;
    int numNeighborParcel;
    Parcel *neighborParcels[maxNumNeighborParcel];
    for (auto parcel : parcelManager.parcels()) {
        // Check the validity of the linear assumption.
        field<SpaceCoord> &x = parcel->skeletonPoints().spaceCoords(timeIdx);
        const field<BodyCoord> &y = SkeletonPoints::bodyCoords();
        SpaceCoord X(domain->numDim());
        double bias = -1.0e+33;
        for (uword i = 0; i < x.size(); ++i) {
            parcel->calcSpaceCoord(timeIdx, y[i], X);
            double bias0 = domain->calcDistance(x[i], X)/parcel->shapeSize(timeIdx).max();
            if (bias0 > bias) bias = bias0;
        }
        // Set the bias limit based on the filament of the parcel and its volume
        // compared with the reference volume.
        double ratio = parcel->volume(timeIdx)/refVolume;
        double biasLimit = transitionFunction(1, maxBiasLimit,
                                              filamentLimit, minBiasLimit,
                                              ratio*parcel->filament());
        // Check parcel bias.
        bool isDegenerated = true;
        // TODO: Redesign the thresholds.
        if (bias < biasLimit && parcel->filament() < filamentLimit &&
            !parcel->meshIndex(timeIdx).atBoundary(*mesh)) {
            isDegenerated = false;
#ifndef LASM_MAX_MIX
            continue;
#endif
        }
        numNeighborParcel = getNeighborParcels(parcel, neighborParcels,
                                               maxNumNeighborParcel);
        if (numNeighborParcel == 0) {
            parcel->dump(timeIdx, meshAdaptor);
            REPORT_ERROR("Parcel " << parcel->id() << " has no neighbor!");
        }
        // Calcuate the mixing weights.
        vec x0(domain->numDim());
#ifdef LASM_IN_SPHERE
        domain->project(geomtk::SphereDomain::STEREOGRAPHIC, parcel->x(timeIdx),
                        parcel->longAxisVertexSpaceCoord(), x0);
#else
        x0 = parcel->longAxisVertexSpaceCoord()()-parcel->x(timeIdx)();
#endif
        vec x1(domain->numDim());
        vec weights(numNeighborParcel, arma::fill::zeros);
        for (uword i = 0; i < numNeighborParcel; ++i) {
#ifdef LASM_IN_SPHERE
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC, parcel->x(timeIdx),
                            neighborParcels[i]->x(timeIdx), x1);
#else
            x1 = domain->diffCoord(neighborParcels[i]->x(timeIdx), parcel->x(timeIdx));
#endif
            double cosTheta = fmin(1, fmax(-1, norm_dot(x0, x1)));
            double sinTheta = sqrt(1-cosTheta*cosTheta);
            double n = norm(x1, 2);
            double d1 = n*cosTheta;
            double d2 = n*sinTheta;
            weights[i] = radialMixing*d1*d1+lateralMixing*d2*d2;
        }
        weights /= max(weights);
        for (uword i = 0; i < numNeighborParcel; ++i) {
            weights[i] = exp(-10*weights[i]);
        }
        double sumWeights = sum(weights);
        if (sumWeights < 1.0e-14) {
            for (uword i = 0; i < numNeighborParcel; ++i) {
                cout << setw(20) << neighborParcels[i]->id() << setw(30) << setprecision(15) << weights[i] << endl;
            }
            parcel->dump(timeIdx, meshAdaptor);
            cout << "Parcel " << parcel->id() << " failed to mix!" << endl;
            exit(-1);
            continue;
        }
#ifdef LASM_USE_DIAG
        Diagnostics::metric<Field<int> >("nmp")(parcel->hostCellIndex())++;
#endif
        weights /= sumWeights;
        // Caclulate total mass.
        double totalMass[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            totalMass[t] = parcel->tracers().mass(t);
            for (uword i = 0; i < numNeighborParcel; ++i) {
                if (weights[i] == 0) continue;
                totalMass[t] += neighborParcels[i]->tracers().mass(t);
            }
        }
        // Caclulate weighted total volume.
        double weightedTotalVolume = parcel->volume(timeIdx);
        for (uword i = 0; i < numNeighborParcel; ++i) {
            if (weights[i] == 0) continue;
            weightedTotalVolume += neighborParcels[i]->volume(timeIdx)*weights[i];
        }
        // Calculate weighted mean tracer density.
        double weightedMeanDensity[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            weightedMeanDensity[t] = parcel->tracers().mass(t);
            for (uword i = 0; i < numNeighborParcel; ++i) {
                if (weights[i] == 0) continue;
                weightedMeanDensity[t] += neighborParcels[i]->tracers().mass(t)*weights[i];
            }
            weightedMeanDensity[t] /= weightedTotalVolume;
        }
        // Restore tracer density to mean density.
        double newTotalMass[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            parcel->tracers().density(t) += restoreFactor*(weightedMeanDensity[t]-parcel->tracers().density(t));
            parcel->tracers().mass(t) = parcel->tracers().density(t)*parcel->volume(timeIdx);
            newTotalMass[t] = parcel->tracers().mass(t);
        }
        for (uword i = 0; i < numNeighborParcel; ++i) {
            if (weights[i] == 0) continue;
            double c = restoreFactor*weights[i];
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                neighborParcels[i]->tracers().density(t) += c*(weightedMeanDensity[t]-neighborParcels[i]->tracers().density(t));
                neighborParcels[i]->tracers().mass(t) = neighborParcels[i]->tracers().density(t)*neighborParcels[i]->volume(timeIdx);
                newTotalMass[t] += neighborParcels[i]->tracers().mass(t);
            }
        }
        // Fix mass inconservation due to floating point inaccuracy.
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            if (totalMass[t] == 0) {
                if (newTotalMass[t] != 0) {
                    REPORT_ERROR("totalMass[" << t << "] is zero, but " <<
                                 "newTotalMass[" << t << "] is not!");
                }
                continue;
            }
            double fixer = totalMass[t]/newTotalMass[t];
            if (fixer == 1) continue;
            if (fixer < 0.99) {
                cout << numNeighborParcel << endl;
                cout << setw(30) << setprecision(15) << totalMass[t] << endl;
                cout << setw(30) << setprecision(15) << newTotalMass[t] << endl;
                cout << setw(30) << setprecision(15) << weightedMeanDensity[t] << endl;
                parcel->dump(timeIdx, meshAdaptor);
                REPORT_ERROR("Mass error (" << fixer <<
                             ") is too large when mixing parcels (" <<
                             parcel->id() << ")!");
            }
            parcel->tracers().mass(t) *= fixer;
            parcel->tracers().density(t) = parcel->tracers().mass(t)/parcel->volume(timeIdx);
            for (uword i = 0; i < numNeighborParcel; ++i) {
                neighborParcels[i]->tracers().mass(t) *= fixer;
                neighborParcels[i]->tracers().density(t) =
                    neighborParcels[i]->tracers().mass(t)/neighborParcels[i]->volume(timeIdx);
            }
        }
        // Change problematic parcel shape (make parcel more uniform).
        if (!isDegenerated) continue;
        auto S = parcel->S();
        if (parcel->meshIndex(timeIdx).atBoundary(*mesh)) {
            S.fill(pow(parcel->detH(timeIdx), 1.0/domain->numDim()));
        } else {
            double filament = parcel->filament();
            double factor = -1;
            while (filament > filamentLimit || factor == -1) {
                factor = transitionFunction(1, 0, filamentLimit, reshapeFactor,
                                            filament);
                double a = sqrt(factor+(1-factor)/filament);
                if (domain->numDim() == 2) {
                    S[0] *= a;
                    S[1] /= a;
                } else if (domain->numDim() == 3) {
                    double b = (S[0]-S[1])/(S[0]-S[2]);
                    double S0 = S[0]*a;
                    double S2 = S[2]/a;
                    double S1 = S0-b*(S0-S2);
                    double s = pow(parcel->volume(timeIdx)/(S0*S1*S2), 1.0/3.0);
                    S[0] = S0*s;
                    S[1] = S1*s;
                    S[2] = S2*s;
                }
                filament = S[0]/S.min();
            }
        }
        parcel->updateDeformMatrix(timeIdx, S);
        parcel->resetSkeletonPoints(timeIdx, *mesh);
    }
} // mixParcels

void AdvectionManager::
remapDensityFromParcelsToGrids(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetTracers(timeIdx);
    vector<int> voidCellIdxs;
    for (uword i = 0; i < gridCoords.n_cols; ++i) {
        int cellIdx = meshAdaptor.cellIndex(i);
        if (meshAdaptor.numConnectedParcel(cellIdx) == 0) {
            voidCellIdxs.push_back(cellIdx);
            continue;
        }
        double totalWeight = 0;
        for (uword j = 0; j < meshAdaptor.numConnectedParcel(cellIdx); ++j) {
            totalWeight += meshAdaptor.remapWeight(cellIdx, j);
        }
        for (uword j = 0; j < meshAdaptor.numConnectedParcel(cellIdx); ++j) {
            Parcel *parcel = meshAdaptor.connectedParcels(cellIdx)[j];
            double weight = meshAdaptor.remapWeight(cellIdx, parcel)/totalWeight;
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                meshAdaptor.density(timeIdx, t, cellIdx) += parcel->tracers().density(t)*weight;
            }
        }
    }
    // Fill void grids.
    for (int i = 0; i < voidCellIdxs.size(); ++i) {
        int cellIdx = voidCellIdxs[i];
        SearchType search(gridTree, NULL, gridCoords, meshAdaptor.coord(cellIdx).cartCoord(), true);
        double searchRadius = mesh->cellSize(CENTER, cellIdx).max();
        vector<vector<size_t> > neighbors;
        vector<vector<double> > distances;
        while (true) {
            mlpack::math::Range r(0, searchRadius);
            search.Search(r, neighbors, distances);
            if (neighbors[0].size() != 0) {
                vector<int> neighborCellIdxs;
                for (int j = 0; j < neighbors[0].size(); ++j) {
                    int neighborCellIdx = meshAdaptor.cellIndex(neighbors[0][j]);
                    // check if the cell is not a void one
                    if (find(voidCellIdxs.begin(), voidCellIdxs.end(),
                             neighborCellIdx) != voidCellIdxs.end() ||
                        meshAdaptor.numConnectedParcel(neighborCellIdx) == 0) continue;
                    neighborCellIdxs.push_back(neighborCellIdx);
                }
                if (neighborCellIdxs.size() == 0) {
                    searchRadius *= 2;
                    continue;
                }
                vec weights(neighborCellIdxs.size());
                for (int i = 0; i < neighborCellIdxs.size(); ++i) {
                    double d = domain->calcDistance(meshAdaptor.coord(cellIdx),
                                                    meshAdaptor.coord(neighborCellIdxs[i]));
                    weights[i] = 1/d;
                }
                weights /= sum(weights);
                for (int i = 0; i < neighborCellIdxs.size(); ++i) {
                    for (int t = 0; t < Tracers::numTracer(); ++t) {
                        meshAdaptor.density(timeIdx, t, cellIdx) +=
                        meshAdaptor.density(timeIdx, t, neighborCellIdxs[i])*weights[i];
                    }
                }
                break;
            } else {
                // Increase the search radius and search again.
                searchRadius *= 2;
            }
        }
    }
#ifdef LASM_USE_RLL_MESH
    for (uword k = mesh->ks(FULL); k <= mesh->ke(FULL); ++k) {
        for (uword i = mesh->is(FULL)+1; i <= mesh->ie(FULL); ++i) {
            uword j1 = mesh->js(FULL);
            uword j2 = mesh->je(FULL);
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                meshAdaptor.density(t)(timeIdx, i, j1, k) =
                meshAdaptor.density(t)(timeIdx, mesh->is(FULL), j1, k);
                meshAdaptor.density(t)(timeIdx, i, j2, k) =
                meshAdaptor.density(t)(timeIdx, mesh->is(FULL), j2, k);
            }
        }
    }
#endif
    // Correct the total mass on the grids.
    vec scale(Tracers::numTracer(), arma::fill::zeros);
    for (uword i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            scale[t] += meshAdaptor.mass(timeIdx, t, i);
        }
    }
    scale = Tracers::totalMasses()/scale;
    for (uword i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            meshAdaptor.density(timeIdx, t, i) *= scale[t];
        }
    }
} // remapDensityFromParcelsToGrids

void AdvectionManager::
remapDensityFromGridsToParcels(const TimeLevelIndex<2> &timeIdx) {
    parcelManager.resetDensities();
    for (auto parcel : parcelManager.parcels()) {
        double totalWeight = 0;
        for (uword i = 0; i < parcel->numConnectedCell(); ++i) {
            int cellIdx = parcel->connectedCellIndexs()[i];
            totalWeight += meshAdaptor.remapWeight(cellIdx, parcel);
        }
        for (uword i = 0; i < parcel->numConnectedCell(); ++i) {
            int cellIdx = parcel->connectedCellIndexs()[i];
            double weight = meshAdaptor.remapWeight(cellIdx, parcel)/totalWeight;
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                parcel->tracers().density(t) += meshAdaptor.density(timeIdx, t, cellIdx)*weight;
            }
        }
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            parcel->tracers().mass(t) = parcel->tracers().density(t)*parcel->volume(timeIdx);
        }
    }
    // Correct the total mass on the parcels.
    vec scale(Tracers::numTracer(), arma::fill::zeros);
    for (auto parcel : parcelManager.parcels()) {
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            scale[t] += parcel->tracers().mass(t);
        }
    }
    scale = Tracers::totalMasses()/scale;
    for (auto parcel : parcelManager.parcels()) {
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            parcel->tracers().mass(t) *= scale[t];
        }
        // Reset the parcel volume, because we would like to respect the density.
        parcel->volume(timeIdx) = parcel->tracers().mass(0)/parcel->tracers().density(0);
        parcel->updateDeformMatrix(timeIdx);
        parcel->resetSkeletonPoints(timeIdx, *mesh);
    }
} // remapDensityFromGridsToParcels

void AdvectionManager::
remapTendencyFromGridsToParcels(const TimeLevelIndex<2> &timeIdx) {
    parcelManager.resetTendencies();
    for (uword i = 0; i < gridCoords.n_cols; ++i) {
        int cellIdx = meshAdaptor.cellIndex(i);
        double totalWeight = 0;
        for (uword j = 0; j < meshAdaptor.numConnectedParcel(cellIdx); ++j) {
            totalWeight += meshAdaptor.remapWeight(cellIdx, j);
        }
        for (uword j = 0; j < meshAdaptor.numConnectedParcel(cellIdx); ++j) {
            Parcel *parcel = meshAdaptor.connectedParcels(cellIdx)[j];
            double weight = meshAdaptor.remapWeight(cellIdx, parcel)/totalWeight;
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                parcel->tracers().tendency(t) += meshAdaptor.tendency(t, cellIdx)*weight;
            }
        }
    }
} // remapTendencyFromGridsToParcels

int AdvectionManager::
getNeighborParcels(Parcel *parcel, Parcel **neighborParcels,
                   int maxNumNeighborParcel) const {
    const vector<int> &connectedCellIndexs = parcel->connectedCellIndexs();
    int numNeighborParcel = 0;
    for (uword i = 0; i < parcel->numConnectedCell(); ++i) {
        for (uword j = 0; j < meshAdaptor.numConnectedParcel(connectedCellIndexs[i]); ++j) {
            Parcel *parcel1 = meshAdaptor.connectedParcels(connectedCellIndexs[i])[j];
            bool alreadyAdded = false;
            if (parcel1 != parcel) {
                for (uword k = 0; k < numNeighborParcel; ++k) {
                    if (parcel1 == neighborParcels[k]) {
                        alreadyAdded = true;
                        break;
                    }
                }
                if (!alreadyAdded) {
                    if (numNeighborParcel == maxNumNeighborParcel) {
                        REPORT_ERROR("Exceed maximum number of neightbor parcels for parcel " << parcel->id() << "!");
                    }
                    neighborParcels[numNeighborParcel] = parcel1;
                    numNeighborParcel++;
                }
            }
        }
    }
    return numNeighborParcel;
} // getNeighborParcels

} // lasm