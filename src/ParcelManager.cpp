#include "ParcelManager.h"
#include "Parcel.h"
#include "SkeletonPoints.h"
#include "QuadraturePoints.h"
#include "Tracers.h"

namespace lasm {

ParcelManager::
ParcelManager() {

}

ParcelManager::
~ParcelManager() {
    Parcels::iterator p = _parcels.begin();
    for (; p != _parcels.end(); ++p) {
        delete *p;
    }
}

// TODO: Support to only set parcels within limited region.
void ParcelManager::
init(const Mesh &mesh) {
    this->mesh = &mesh;
    TimeLevelIndex<2> timeIdx;
#ifdef LASM_USE_RLL_MESH
    vec areas(mesh.numGrid(1, FULL));
    vec areaWeights(mesh.numGrid(1, FULL));
    for (uword j = mesh.js(FULL); j <= mesh.je(FULL); ++j) {
        int cellIdx = mesh.wrapIndex(CENTER, mesh.is(FULL), j, mesh.ks(FULL));
        areas[j] = mesh.cellVolume(cellIdx);
        areaWeights[j] = exp(-0.1*fabs(mesh.sinLat(FULL, j)));
    }
    double maxArea = areas.max();
    SpaceCoord x(mesh.domain().numDim());
    int id = 0;
    vec size(mesh.domain().numDim());
    for (uword k = mesh.ks(FULL); k <= mesh.ke(FULL); ++k) {
        if (mesh.domain().numDim() == 3) {
            size[2] = mesh.gridInterval(2, HALF, k);
        }
        for (uword j = mesh.js(FULL)+1; j < mesh.je(FULL); ++j) {
            size[1] = mesh.gridInterval(1, HALF, j)*mesh.domain().radius();
            int numReducedLon = mesh.numGrid(0, FULL)/
                std::max(1, static_cast<int>(maxArea*areaWeights[j]/areas[j]));
            double dlon = PI2/numReducedLon;
            for (uword i = 0; i < numReducedLon; ++i) {
                double lon = i*dlon;
                if (mesh.domain().numDim() == 2) {
                    x.set(lon, mesh.lat(FULL, j));
                } else if (mesh.domain().numDim() == 3) {
                    x.set(lon, mesh.lat(FULL, j), mesh.lev(FULL, k));
                }
                x.transformToCart(mesh.domain());
                size[0] = dlon*mesh.domain().radius()*mesh.cosLat(FULL, j);
                Parcel *parcel = new Parcel;
                parcel->init(id++);
                parcel->x(timeIdx) = x;
                parcel->meshIndex(timeIdx).locate(mesh, x);
                parcel->skeletonPoints().init(mesh, size);
                parcel->updateDeformMatrix(timeIdx);
                parcel->tracers().init();
                _parcels.push_back(parcel);
            }
        }
    }
#else
    for (uword i = 0; i < mesh.totalNumGrid(CENTER, mesh.domain().numDim()); ++i) {
        const SpaceCoord &x = mesh.gridCoord(CENTER, i);
        Parcel *parcel = new Parcel;
        parcel->init(i);
        parcel->x(timeIdx) = x;
        parcel->meshIndex(timeIdx).locate(mesh, x);
        auto cellSize = mesh.cellSize(CENTER, i);
        parcel->skeletonPoints().init(mesh, cellSize);
        parcel->volume(timeIdx) = mesh.cellVolume(i);
        parcel->updateDeformMatrix(timeIdx);
        parcel->resetSkeletonPoints(timeIdx, mesh);
        parcel->tracers().init();
        _parcels.push_back(parcel);
    }
#endif
} // init

void ParcelManager::
output(const TimeLevelIndex<2> &timeIdx, int ncId) const {
    int parcelDimId, skel1DimId, dimDimId, tracerDimId;
    int idVarId;
    int cDimIds[2], cVarId;
    int mDimIds[2], mVarId;
    int sDimIds[3], s1VarId;
#define OUTPUT_TRACER_SHAPE
#ifdef OUTPUT_TRACER_SHAPE
    int skel2DimId, s2VarId, numSkel2 = 40;
#endif
    char str[100];
    int l;
    int *intData;
    double *doubleData;

    nc_redef(ncId);

    nc_def_dim(ncId, "parcel", _parcels.size(), &parcelDimId);
    nc_def_dim(ncId, "dim", mesh->domain().numDim(), &dimDimId);
    nc_def_dim(ncId, "tracer", Tracers::numTracer(), &tracerDimId);
    nc_def_dim(ncId, "skel1", SkeletonPoints::numPoint(), &skel1DimId);

#ifdef OUTPUT_TRACER_SHAPE
    if (mesh->domain().numDim() == 2) {
        // Only output skeleton in 2D domain, since in 3D it could be messy.
        nc_def_dim(ncId, "skel2", numSkel2, &skel2DimId);
    }
#endif
    
    nc_def_var(ncId, "id", NC_INT, 1, &parcelDimId, &idVarId);
    sprintf(str, "parcel identifier");
    nc_put_att(ncId, idVarId, "long_name", NC_CHAR, strlen(str), str);

    cDimIds[0] = parcelDimId; cDimIds[1] = dimDimId;
    nc_def_var(ncId, "c", NC_DOUBLE, 2, cDimIds, &cVarId);
    sprintf(str, "parcel centroid coordinates on %s", mesh->domain().brief().c_str());
    nc_put_att(ncId, cVarId, "long_name", NC_CHAR, strlen(str), str);

    mDimIds[0] = parcelDimId; mDimIds[1] = tracerDimId;
    nc_def_var(ncId, "m", NC_DOUBLE, 2, mDimIds, &mVarId);
    sprintf(str, "tracer mass");
    nc_put_att(ncId, mVarId, "long_name", NC_CHAR, strlen(str), str);

    sDimIds[0] = parcelDimId; sDimIds[1] = skel1DimId; sDimIds[2] = dimDimId;
    nc_def_var(ncId, "s1", NC_DOUBLE, 3, sDimIds, &s1VarId);
    sprintf(str, "parcel actual skeleton");
    nc_put_att(ncId, s1VarId, "long_name", NC_CHAR, strlen(str), str);

#ifdef OUTPUT_TRACER_SHAPE
    if (mesh->domain().numDim() == 2) {
        sDimIds[1] = skel2DimId;
        nc_def_var(ncId, "s2", NC_DOUBLE, 3, sDimIds, &s2VarId);
        sprintf(str, "parcel fitted skeleton");
        nc_put_att(ncId, s2VarId, "long_name", NC_CHAR, strlen(str), str);
    }
#endif
    
    nc_enddef(ncId);

    intData = new int[_parcels.size()];
    l = 0;
    for (auto parcel : _parcels) {
        intData[l++] = parcel->id();
    }
    nc_put_var(ncId, idVarId, intData);
    delete [] intData;
    
    doubleData = new double[_parcels.size()*mesh->domain().numDim()];
    l = 0;
    for (auto parcel : _parcels) {
        for (uword m = 0; m < mesh->domain().numDim(); ++m) {
            doubleData[l++] = parcel->x(timeIdx)(m);
        }
    }
    nc_put_var(ncId, cVarId, doubleData);
    delete [] doubleData;
    
    doubleData = new double[_parcels.size()*Tracers::numTracer()];
    l = 0;
    for (auto parcel : _parcels) {
        for (int t = 0; t < Tracers::numTracer(); ++t) {
            doubleData[l++] = parcel->tracers().mass(t);
        }
    }
    nc_put_var(ncId, mVarId, doubleData);
    delete [] doubleData;

    doubleData = new double[_parcels.size()*SkeletonPoints::numPoint()*mesh->domain().numDim()];
    l = 0;
    for (auto parcel : _parcels) {
        for (auto xs : parcel->skeletonPoints().spaceCoords(timeIdx)) {
            for (uword m = 0; m < mesh->domain().numDim(); ++m) {
                doubleData[l++] = xs(m);
            }
        }
    }
    nc_put_var(ncId, s1VarId, doubleData);
    delete [] doubleData;

#ifdef OUTPUT_TRACER_SHAPE
    if (mesh->domain().numDim() == 2) {
        double dtheta = PI2/numSkel2;
        BodyCoord y(2); SpaceCoord x(2);
        doubleData = new double[_parcels.size()*numSkel2*mesh->domain().numDim()];
        l = 0;
        for (auto parcel : _parcels) {
            for (int i = 0; i < numSkel2; ++i) {
                double theta = i*dtheta;
                y(0) = cos(theta);
                y(1) = sin(theta);
                parcel->calcSpaceCoord(timeIdx, y, x);
                for (uword m = 0; m < mesh->domain().numDim(); ++m) {
                    doubleData[l++] = x(m);
                }
            }
        }
        nc_put_var(ncId, s2VarId, doubleData);
        delete [] doubleData;
    }
#endif
} // output

void ParcelManager::
resetDensities() {
    for (auto parcel : _parcels) {
        parcel->tracers().resetDensities();
    }
} // resetDensities

void ParcelManager::
resetTendencies() {
    for (auto parcel : _parcels) {
        parcel->tracers().resetTendencies();
    }
} // resetTendencies

} // lasm