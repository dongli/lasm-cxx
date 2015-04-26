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
    for (uword i = 0; i < mesh.totalNumGrid(CENTER, mesh.domain().numDim()); ++i) {
        const SpaceCoord &x = mesh.gridCoord(CENTER, i);
        Parcel *parcel = new Parcel;
        parcel->init(i);
        parcel->x(timeIdx) = x;
        parcel->meshIndex(timeIdx).locate(mesh, x);
        vec sizes(mesh.domain().numDim());
#if defined LASM_IN_CARTESIAN
        // TODO: Find a better way to specify grid sizes.
        int I, J, K;
        mesh.unwrapIndex(CENTER, i, I, J, K);
        sizes[0] = mesh.gridInterval(0, HALF, I-1);
        sizes[1] = mesh.gridInterval(1, HALF, J-1);
        sizes[2] = mesh.gridInterval(2, HALF, K);
#endif
        parcel->skeletonPoints().init(mesh, sizes);
        parcel->updateDeformMatrix(timeIdx);
        parcel->quadraturePoints().updateSpaceCoords(timeIdx);
        parcel->tracers().init();
        _parcels.push_back(parcel);
    }
} // init

void ParcelManager::
output(const TimeLevelIndex<2> &timeIdx, int ncId) {
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
        intData[l++] = parcel->ID();
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
            doubleData[l++] = parcel->tracers().mass(timeIdx, t);
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
resetTracers(const TimeLevelIndex<2> &timeIdx) {
    for (auto parcel : _parcels) {
        parcel->tracers().reset(timeIdx);
    }
}

} // lasm