#include "MeshAdaptor.h"
#include "Parcel.h"

namespace lasm {

MeshAdaptor::
MeshAdaptor() {
}

MeshAdaptor::
~MeshAdaptor() {
    for (uword t = 0; t < _densities.size(); ++t) {
        delete _densities[t];
    }
}

void MeshAdaptor::
init(const Mesh &mesh) {
    _mesh = &mesh;
    _numConnectedParcel.resize(mesh.totalNumGrid(CENTER));
    _connectedParcels.resize(mesh.totalNumGrid(CENTER));
    _remapWeights.resize(mesh.totalNumGrid(CENTER));
    _numContainedParcel.resize(mesh.totalNumGrid(CENTER));
    _containedParcels.resize(mesh.totalNumGrid(CENTER));
#ifdef LASM_USE_RLL_MESH
    _cellIndexs.resize(mesh.totalNumGridWithUniquePoleGrid(CENTER));
    uword j = 0;
    for (uword i = 0; i < mesh.totalNumGrid(CENTER); ++i) {
        auto spanIdx = mesh.unwrapIndex(CENTER, i);
        if ((spanIdx[1] == mesh.js(FULL) && spanIdx[0] != mesh.is(FULL)) ||
            (spanIdx[1] == mesh.je(FULL) && spanIdx[0] != mesh.is(FULL))) {
            continue;
        }
        _cellIndexs[j++] = i;
    }
#else
    _cellIndexs.resize(mesh.totalNumGrid(CENTER));
    for (uword i = 0; i < mesh.totalNumGrid(CENTER); ++i) {
        _cellIndexs[i] = i;
    }
#endif
} // init

void MeshAdaptor::
addTracer(const string &name, const string &unit, const string &comment) {
    _densities.push_back(new Field<double, 2>);
    _densities.back()->create(name, unit, comment,
                              *_mesh, CENTER, _mesh->domain().numDim());
    _tendencies.push_back(new Field<double>);
    _tendencies.back()->create("d"+name, unit+" s-1", "mass tendency of "+comment,
                               *_mesh, CENTER, _mesh->domain().numDim());
} // addTracer

void MeshAdaptor::
resetTracers(const TimeLevelIndex<2> &timeIdx) {
    for (uword t = 0; t < _densities.size(); ++t) {
        for (uword i = 0; i < _mesh->totalNumGrid(CENTER, _mesh->domain().numDim()); ++i) {
            (*_densities[t])(timeIdx, i) = 0;
        }
    }
} // resetTracers

void MeshAdaptor::
connectParcel(int cellIdx, Parcel *parcel, double weight) {
    for (uword i = 0; i < _numConnectedParcel[cellIdx]; ++i) {
        if (_connectedParcels[cellIdx][i] == parcel) {
            return;
        }
    }
    if (_numConnectedParcel[cellIdx] == _connectedParcels[cellIdx].size()) {
        _connectedParcels[cellIdx].push_back(parcel);
        _remapWeights[cellIdx].push_back(weight);
    } else {
        _connectedParcels[cellIdx][_numConnectedParcel[cellIdx]] = parcel;
        _remapWeights[cellIdx][_numConnectedParcel[cellIdx]] = weight;
    }
    _numConnectedParcel[cellIdx]++;
#ifdef LASM_USE_DIAG
    Diagnostics::metric<Field<int> >("ncp2")(cellIdx)++;
#endif
} // connectParcel

double MeshAdaptor::
remapWeight(int cellIdx, Parcel *parcel) const {
    for (uword i = 0; i < _numConnectedParcel[cellIdx]; ++i) {
        if (_connectedParcels[cellIdx][i] == parcel) {
            return _remapWeights[cellIdx][i];
        }
    }
    REPORT_ERROR("Parcel is not connected!");
} // remapWeight

double& MeshAdaptor::
remapWeight(int cellIdx, Parcel *parcel) {
    for (uword i = 0; i < _numConnectedParcel[cellIdx]; ++i) {
        if (_connectedParcels[cellIdx][i] == parcel) {
            return _remapWeights[cellIdx][i];
        }
    }
    REPORT_ERROR("Parcel is not connected!");
} // remapWeight

void MeshAdaptor::
containParcel(int cellIdx, Parcel *parcel) {
#ifndef NDEBUG
    for (uword i = 0; i < _numContainedParcel[cellIdx]; ++i) {
        if (_containedParcels[cellIdx][i] == parcel) {
            REPORT_ERROR("Parcel (id = " << parcel->id() <<
                         ") has already been contained!");
        }
    }
#endif
    if (_numContainedParcel[cellIdx] == _containedParcels[cellIdx].size()) {
        _containedParcels[cellIdx].push_back(parcel);
    } else {
        _containedParcels[cellIdx][_numContainedParcel[cellIdx]] = parcel;
    }
    parcel->hostCellIndex() = cellIdx;
    _numContainedParcel[cellIdx]++;
#ifdef LASM_USE_DIAG
    Diagnostics::metric<Field<int> >("ncp1")(cellIdx)++;
#endif
} // containParcel

void MeshAdaptor::
resetConnectedAndContainedParcels() {
    for (uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
        _numConnectedParcel[i] = 0;
        _numContainedParcel[i] = 0;
    }
} // resetConnectedAndContainedParcels

} // lasm