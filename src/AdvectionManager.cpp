#include "AdvectionManager.h"
#include "ShapeFunction.h"
#include "QuadraturePoints.h"
#include "SkeletonPoints.h"
#include "Parcel.h"
#include "Tracers.h"

namespace lasm {

AdvectionManager::
AdvectionManager() : geomtk::AdvectionManagerInterface<2, Domain, Mesh, VelocityField>() {
    minBiasLimit = 0.1;
    maxBiasLimit = 0.5;
}

AdvectionManager::
~AdvectionManager() {

}

void AdvectionManager::
init(const ConfigManager &configManager, const Mesh &mesh) {
    domain = &mesh.domain();
    this->mesh = &mesh;
    // Initialize objects.
    ShapeFunction::init(*domain);
    QuadraturePoints::init(*domain);
    SkeletonPoints::init(*domain);
    Parcel::init(*domain);
    Regrid::init(mesh);
    parcelManager.init(mesh);
    meshAdaptor.init(mesh);
    // Initialize tree structure of mesh for searching.
#if defined LASM_IN_SPHERE
    cellCoords.reshape(3, mesh.totalNumGrid(CENTER, domain->numDim()));
#elif defined LASM_IN_CARTESIAN
    cellCoords.reshape(domain->numDim(), mesh.totalNumGrid(CENTER, domain->numDim()));
#endif
    for (int i = 0; i < mesh.totalNumGrid(CENTER, domain->numDim()); ++i) {
        cellCoords.col(i) = meshAdaptor.coord(i).cartCoord();
    }
    cellTree = new Tree(cellCoords, cellCoordsMap);
    cellCoords = cellTree->Dataset();
    // Connect parcels and grids Initially.
    TimeLevelIndex<2> timeIdx;
    connectParcelsAndGrids(timeIdx);
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
} // advance

void AdvectionManager::
input(const TimeLevelIndex<2> &timeIdx, double *q) {
    // Copy input tracer density onto internal mesh grids.
    int l = 0;
    for (int t = 0; t < Tracers::numTracer(); ++t) {
        for (int i = 0; i < mesh->totalNumGrid(CENTER, domain->numDim()); ++i) {
            meshAdaptor.mass(timeIdx, t, i) = q[l++]*meshAdaptor.volume(i);
        }
    }
    // Transfer tracer mass from grids to tracers.
    remapFromGridsToParcels(timeIdx);
    // Transfer back to grids (this will be different from input).
    remapFromParcelsToGrids(timeIdx);
} // input

void AdvectionManager::
output(const TimeLevelIndex<2> &timeIdx, int ncId) {
    // Append global attributes.
    nc_redef(ncId);
    nc_put_att(ncId, NC_GLOBAL, "filament_limit", NC_DOUBLE, 1, &filamentLimit);
    nc_put_att(ncId, NC_GLOBAL, "radial_mixing", NC_DOUBLE, 1, &radialMixing);
    nc_put_att(ncId, NC_GLOBAL, "lateral_mixing", NC_DOUBLE, 1, &lateralMixing);
    nc_put_att(ncId, NC_GLOBAL, "restore_factor", NC_DOUBLE, 1, &restoreFactor);

    string str = to_iso_extended_string(boost::gregorian::day_clock::universal_day());
    nc_put_att(ncId, NC_GLOBAL, "create_date", NC_CHAR, str.length(), str.c_str());
    nc_enddef(ncId);
    // Append tracer data.
    parcelManager.output(timeIdx, ncId);
} // output

void AdvectionManager::
integrate(double dt, const TimeLevelIndex<2> &newIdx, const VelocityField &velocityField) {
    TimeLevelIndex<2> oldIdx = newIdx-1;
    TimeLevelIndex<2> halfIdx = newIdx-0.5;
    double dt05 = 0.5*dt;
    const auto &divergence = velocityField.divergence();
    for (auto parcel : parcelManager.parcels()) {
        Velocity v1(domain->numDim());
        Velocity v2(domain->numDim());
        Velocity v3(domain->numDim());
        Velocity v4(domain->numDim());
        Velocity v(domain->numDim());
        double div, k1, k2, k3, k4, volume;
        SpaceCoord &x0 = parcel->x(oldIdx);
        SpaceCoord &x1 = parcel->x(newIdx);
        MeshIndex &idx0 = parcel->meshIdx(oldIdx);
        MeshIndex &idx1 = parcel->meshIdx(newIdx);
        // Stage 1.
        regrid.run(BILINEAR, oldIdx, velocityField, x0, v1, &idx0);
        regrid.run(BILINEAR, oldIdx, divergence, x0, div, &idx0);
        mesh->move(x0, dt05, v1, idx0, x1);
        idx1.locate(*mesh, x1);
        k1 = -parcel->volume(oldIdx)*div;
        volume = parcel->volume(oldIdx)+dt05*k1;
        // Stage 2.
        regrid.run(BILINEAR, halfIdx, velocityField, x1, v2, &idx1);
        regrid.run(BILINEAR, halfIdx, divergence, x1, div, &idx1);
        mesh->move(x0, dt05, v2, idx0, x1);
        idx1.locate(*mesh, x1);
        k2 = -volume*div;
        volume = parcel->volume(oldIdx)+dt05*k2;
        // Stage 3.
        regrid.run(BILINEAR, halfIdx, velocityField, x1, v3, &idx1);
        regrid.run(BILINEAR, halfIdx, divergence, x1, div, &idx1);
        mesh->move(x0, dt, v3, idx0, x1);
        idx1.locate(*mesh, x1);
        k3 = -volume*div;
        volume = parcel->volume(oldIdx)+dt*k3;
        // Stage 4.
        regrid.run(BILINEAR, newIdx, velocityField, x1, v4, &idx1);
        regrid.run(BILINEAR, newIdx, divergence, x1, div, &idx1);
        v = (v1+v2*2.0+v3*2.0+v4)/6.0;
        mesh->move(x0, dt, v, idx0, x1);
        idx1.locate(*mesh, x1);
        k4 = -volume*div;
        parcel->updateVolume(newIdx, parcel->volume(oldIdx)+dt*(k1+2*k2+2*k3+k4)/6);
        // Update the skeleton points.
        field<SpaceCoord> &xs0 = parcel->skeletonPoints().spaceCoords(oldIdx);
        field<SpaceCoord> &xs1 = parcel->skeletonPoints().spaceCoords(newIdx);
        field<MeshIndex> &idxs0 = parcel->skeletonPoints().meshIdxs(oldIdx);
        field<MeshIndex> &idxs1 = parcel->skeletonPoints().meshIdxs(newIdx);
        for (int i = 0; i < xs0.size(); ++i) {
            // Stage 1.
            regrid.run(BILINEAR, oldIdx, velocityField, xs0[i], v1, &idxs0[i]);
            mesh->move(xs0[i], dt05, v1, idxs0[i], xs1[i]);
            idxs1[i].locate(*mesh, xs1[i]);
            // Stage 2.
            regrid.run(BILINEAR, halfIdx, velocityField, xs1[i], v2, &idxs1[i]);
            mesh->move(xs0[i], dt05, v2, idxs0[i], xs1[i]);
            idxs1[i].locate(*mesh, xs1[i]);
            // Stage 3.
            regrid.run(BILINEAR, halfIdx, velocityField, xs1[i], v3, &idxs1[i]);
            mesh->move(xs0[i], dt, v3, idxs0[i], xs1[i]);
            idxs1[i].locate(*mesh, xs1[i]);
            // Stage 4.
            regrid.run(BILINEAR, newIdx, velocityField, xs1[i], v4, &idxs1[i]);
            v = (v1+v2*2.0+v3*2.0+v4)/6.0;
            mesh->move(xs0[i], dt, v, idxs0[i], xs1[i]);
            idxs1[i].locate(*mesh, xs1[i]);
        }
        parcel->updateDeformMatrix(newIdx);
    }
} // integrate

void AdvectionManager::
connectParcelAndGrids(const TimeLevelIndex<2> &timeIdx, Parcel *parcel) {
    parcel->updateShapeSize(timeIdx);
    Searcher a(cellTree, NULL, cellCoords,
               parcel->x(timeIdx).cartCoord(), true);
    double longAxisSize = parcel->shapeSize(timeIdx).max();
    mlpack::math::Range r(0.0, longAxisSize);
    vector<vector<size_t> > neighbors;
    vector<vector<double> > distances;
    a.Search(r, neighbors, distances);
    BodyCoord y(domain->numDim());
    for (int i = 0; i < neighbors[0].size(); ++i) {
        int j = cellCoordsMap[neighbors[0][i]];
        parcel->calcBodyCoord(timeIdx, meshAdaptor.coord(j), y);
        // TODO: Change the weight calculation.
        double w = parcel->shapeFunction(timeIdx, y);
        if (w > 0.0) {
            meshAdaptor.connectParcel(j, parcel, w);
            parcel->connectCell(j);
        }
    }
} // connectParcelAndGrids

void AdvectionManager::
connectParcelsAndGrids(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetConnectedParcels();
    for (auto parcel : parcelManager.parcels()) {
        parcel->resetConnectedCells();
        connectParcelAndGrids(timeIdx, parcel);
    }
} // connectParcelsAndGrids

void AdvectionManager::
mixParcels(const TimeLevelIndex<2> &timeIdx) {
    for (auto parcel : parcelManager.parcels()) {
        // Check the validity of the linear assumption.
        field<SpaceCoord> &x = parcel->skeletonPoints().spaceCoords(timeIdx);
        const field<BodyCoord> &y = SkeletonPoints::bodyCoords();
        SpaceCoord X(2);
        double bias = -1.0e+33;
        for (int i = 0; i < x.size(); ++i) {
            parcel->calcSpaceCoord(timeIdx, y[i], X);
            double bias0 = domain->calcDistance(x[i], X)/parcel->shapeSize(timeIdx).max();
            if (bias0 > bias) bias = bias0;
        }
        // Set the bias limit based on the filament of the parcel and its volume
        // compared with the neighbor tracers.
        vector<Parcel*> neighborParcels = getNeighborParcels(parcel);
        double meanVolume = parcel->volume(timeIdx);
        for (int i = 0; i < neighborParcels.size(); ++i) {
            meanVolume += neighborParcels[i]->volume(timeIdx);
        }
        meanVolume /= neighborParcels.size()+1;
        double ratio = parcel->volume(timeIdx)/meanVolume;
        double biasLimit = transitionFunction(1, maxBiasLimit, 5, minBiasLimit,
                                              ratio*parcel->filamentDegree(timeIdx));
        // Check parcel bias.
        bool isDegenerated = true;
        if (bias < biasLimit && parcel->filamentDegree(timeIdx) < filamentLimit) {
            isDegenerated = false;
#ifndef LASM_TEST_ALL_MIX
            continue;
#endif
        }
        // Calcuate the mixing weights.
        vec x0(2);
#ifdef LASM_SPHERE_DOMAIN
        domain->project(geomtk::SphereDomain::STEREOGRAPHIC, parcel->x(timeIdx),
                        parcel->longAxisVertexSpaceCoord(), x0);
#else
        x0 = parcel->longAxisVertexSpaceCoord()()-parcel->x(timeIdx)();
#endif
        double n0 = norm(x0, 2);
        vec x1(2);
        vec weights(neighborParcels.size(), arma::fill::zeros);
        for (int i = 0; i < neighborParcels.size(); ++i) {
#ifdef LASM_SPHERE_DOMAIN
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC, parcel->x(timeIdx),
                            neighborParcels[i]->x(timeIdx), x1);
#else
            x1 = domain->diffCoord(neighborParcels[i]->x(timeIdx), parcel->x(timeIdx));
#endif
            x1 /= n0;
            double cosTheta = norm_dot(x0, x1);
            // NOTE: Ensure cosTheta is in the range [-1,1]!
            cosTheta = fmin(1, fmax(-1, cosTheta));
            double sinTheta = sqrt(1-cosTheta*cosTheta);
            double n1 = norm(x1, 2);
            double d1 = n1*cosTheta;
            double d2 = n1*sinTheta;
            weights[i] = exp(-(radialMixing*d1*d1+lateralMixing*d2*d2));
            if (weights[i] < 1.0e-6) {
                weights[i] = 0;
            }
        }
        double sumWeights = sum(weights);
        if (sumWeights < 1.0e-14) {
            continue;
        }
        weights /= sumWeights;
        // Caclulate total mass.
        double totalMass[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            totalMass[t] = parcel->tracers().mass(timeIdx, t);
            for (int i = 0; i < neighborParcels.size(); ++i) {
                if (weights[i] == 0) continue;
                totalMass[t] += neighborParcels[i]->tracers().mass(timeIdx, t);
            }
        }
        // Caclulate weighted total volume.
        double weightedTotalVolume = parcel->volume(timeIdx);
        for (int i = 0; i < neighborParcels.size(); ++i) {
            if (weights[i] == 0) continue;
            weightedTotalVolume += neighborParcels[i]->volume(timeIdx)*weights[i];
        }
        // Calculate weighted mean tracer density.
        double weightedMeanDensity[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            weightedMeanDensity[t] = parcel->tracers().mass(timeIdx, t);
            for (int i = 0; i < neighborParcels.size(); ++i) {
                if (weights[i] == 0) continue;
                weightedMeanDensity[t] += neighborParcels[i]->tracers().mass(timeIdx, t)*weights[i];
            }
            weightedMeanDensity[t] /= weightedTotalVolume;
        }
        // Restore tracer density to mean density.
        double newTotalMass[Tracers::numTracer()];
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            double rho = parcel->tracers().density(timeIdx, t);
            rho += restoreFactor*(weightedMeanDensity[t]-rho);
            parcel->tracers().mass(timeIdx, t) = rho*parcel->volume(timeIdx);
            newTotalMass[t] = parcel->tracers().mass(timeIdx, t);
        }
        for (int i = 0; i < neighborParcels.size(); ++i) {
            if (weights[i] == 0) continue;
            double c = restoreFactor*weights[i];
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                double rho = neighborParcels[i]->tracers().density(timeIdx, t);
                rho += c*(weightedMeanDensity[t]-rho);
                neighborParcels[i]->tracers().mass(timeIdx, t) = rho*neighborParcels[i]->volume(timeIdx);
                newTotalMass[t] += neighborParcels[i]->tracers().mass(timeIdx, t);
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
            parcel->tracers().mass(timeIdx, t) *= fixer;
            for (int i = 0; i < neighborParcels.size(); ++i) {
                neighborParcels[i]->tracers().mass(timeIdx, t) *= fixer;
            }
        }
        // Change problematic parcel shape (make parcel more uniform).
        if (!isDegenerated) continue;
        const double maxUniformFactor = 0.5;
        double a, b;
        double uniformFactor = transitionFunction(1, 1, 5, maxUniformFactor,
                                                  parcel->filamentDegree(timeIdx));
        a = pow(pow(uniformFactor, domain->numDim()-1), 1.0/domain->numDim());
        b = 1/a;
        auto S = parcel->S();
        S[0] *= a;
        for (int i = 1; i < domain->numDim(); ++i) {
            S[i] *= b;
        }
        parcel->updateDeformMatrix(timeIdx, S);
        parcel->resetSkeletonPoints(timeIdx, *mesh);
    }
} // mixParcels

void AdvectionManager::
remapFromParcelsToGrids(const TimeLevelIndex<2> &timeIdx) {
    meshAdaptor.resetTracers(timeIdx);
    for (auto parcel : parcelManager.parcels()) {
        double totalWeight = 0;
        for (int j = 0; j < parcel->numConnectedCell(); ++j) {
            int i = parcel->connectedCellIdxs()[j];
            totalWeight += meshAdaptor.remapWeight(i, parcel);
        }
        for (int j = 0; j < parcel->numConnectedCell(); ++j) {
            int i = parcel->connectedCellIdxs()[j];
            double weight = meshAdaptor.remapWeight(i, parcel)/totalWeight;
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                meshAdaptor.mass(timeIdx, t, i) += parcel->tracers().mass(timeIdx, t)*weight;
            }
        }
    }
    for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        if (meshAdaptor.numConnectedParcel(i) == 0) {
            REPORT_ERROR("Encounter void cell! Handle it.");
        }
    }
}

void AdvectionManager::
remapFromGridsToParcels(const TimeLevelIndex<2> &timeIdx) {
    parcelManager.resetTracers(timeIdx);
    for (int i = 0; i < mesh->totalNumGrid(CENTER); ++i) {
        double totalWeight = 0;
        for (int j = 0; j < meshAdaptor.numConnectedParcel(i); ++j) {
            totalWeight += meshAdaptor.remapWeight(i, meshAdaptor.connectedParcels(i)[j]);
        }
        for (int j = 0; j < meshAdaptor.numConnectedParcel(i); ++j) {
            Parcel *parcel = meshAdaptor.connectedParcels(i)[j];
            double weight = meshAdaptor.remapWeight(i, parcel)/totalWeight;
            for (int t = 0; t < Tracers::numTracer(); ++t) {
                parcel->tracers().mass(timeIdx, t) += meshAdaptor.mass(timeIdx, t, i)*weight;
            }
        }
    }
}

vector<Parcel*> AdvectionManager::
getNeighborParcels(Parcel *parcel) const {
    const vector<int> &connectedCellIdxs = parcel->connectedCellIdxs();
    vector<Parcel*> neighborParcels;
    for (int i = 0; i < parcel->numConnectedCell(); ++i) {
        for (int j = 0; j < meshAdaptor.numConnectedParcel(connectedCellIdxs[i]); ++j) {
            Parcel *parcel1 = meshAdaptor.connectedParcels(connectedCellIdxs[i])[j];
            bool alreadyAdded = false;
            if (parcel1 != parcel) {
                for (int k = 0; k < neighborParcels.size(); ++k) {
                    if (parcel1 == neighborParcels[k]) {
                        alreadyAdded = true;
                        break;
                    }
                }
                if (!alreadyAdded) {
                    neighborParcels.push_back(parcel1);
                }
            }
        }
    }
    return neighborParcels;
} // getNeighborParcels

} // lasm