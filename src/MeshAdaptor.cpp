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
    this->mesh = &mesh;
    _numConnectedParcel.resize(mesh.totalNumGrid(CENTER));
    _connectedParcels.resize(mesh.totalNumGrid(CENTER));
    _remapWeights.resize(mesh.totalNumGrid(CENTER));
    _numContainedParcel.resize(mesh.totalNumGrid(CENTER));
    _containedParcels.resize(mesh.totalNumGrid(CENTER));
} // init

void MeshAdaptor::
addTracer(const string &name, const string &unit, const string &comment) {
    _densities.push_back(new Field<double, 2>);
    _densities.back()->create(name, unit, comment, *mesh, CENTER, mesh->domain().numDim());
} // addTracer

void MeshAdaptor::
resetTracers(const TimeLevelIndex<2> &timeIdx) {
    for (uword t = 0; t < _densities.size(); ++t) {
        for (uword i = 0; i < mesh->totalNumGrid(CENTER, mesh->domain().numDim()); ++i) {
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
} // containParcel

void MeshAdaptor::
resetConnectedAndContainedParcels() {
    for (uword i = 0; i < mesh->totalNumGrid(CENTER, mesh->domain().numDim()); ++i) {
        _numConnectedParcel[i] = 0;
        _numContainedParcel[i] = 0;
    }
} // resetConnectedAndContainedParcels

} // lasm