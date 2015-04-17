#include "MeshAdaptor.h"
#include "Parcel.h"

namespace lasm {

MeshAdaptor::
MeshAdaptor() {
}

MeshAdaptor::
~MeshAdaptor() {
    for (int t = 0; t < _masses.size(); ++t) {
        delete _masses[t];
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
}

void MeshAdaptor::
addTracer(const string &name, const string &unit, const string &comment) {
    _masses.push_back(new Field<double, 2>);
    _masses.back()->create(name, unit, comment, *mesh, CENTER, mesh->domain().numDim());
}

void MeshAdaptor::
resetTracers(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < _masses.size(); ++t) {
        for (int i = 0; i < mesh->totalNumGrid(CENTER, mesh->domain().numDim()); ++i) {
            (*_masses[t])(timeIdx, i);
        }
    }
}

void MeshAdaptor::
connectParcel(int cellIdx, Parcel *parcel, double weight) {
#ifndef NDEBUG
    for (int i = 0; i < _numConnectedParcel[cellIdx]; ++i) {
        if (_connectedParcels[cellIdx][i] == parcel) {
            REPORT_ERROR("Parcel (ID = " << parcel->ID() <<
                         ") has already been connected!");
        }
    }
#endif
    if (_numConnectedParcel[cellIdx] == _connectedParcels[cellIdx].size()) {
        _connectedParcels[cellIdx].push_back(parcel);
        _remapWeights[cellIdx].push_back(weight);
    } else {
        _connectedParcels[cellIdx][_numConnectedParcel[cellIdx]] = parcel;
        _remapWeights[cellIdx][_numConnectedParcel[cellIdx]] = weight;
    }
    _numConnectedParcel[cellIdx]++;
}

void MeshAdaptor::
resetConnectedParcels() {
    for (int i = 0; i < mesh->totalNumGrid(CENTER, mesh->domain().numDim()); ++i) {
        _numConnectedParcel[i] = 0;
    }
}

double MeshAdaptor::
remapWeight(int cellIdx, Parcel *parcel) const {
    for (int i = 0; i < _numConnectedParcel[cellIdx]; ++i) {
        if (_connectedParcels[cellIdx][i] == parcel) {
            return _remapWeights[cellIdx][i];
        }
    }
    REPORT_ERROR("Parcel is not connected!");
}

void MeshAdaptor::
containParcel(int cellIdx, Parcel *parcel) {
#ifndef NDEBUG
    for (int i = 0; i < _numContainedParcel[cellIdx]; ++i) {
        if (_containedParcels[cellIdx][i] == parcel) {
            REPORT_ERROR("Parcel (ID = " << parcel->ID() <<
                         ") has already been contained!");
        }
    }
#endif
    if (_numContainedParcel[cellIdx] == _containedParcels[cellIdx].size()) {
        _containedParcels[cellIdx].push_back(parcel);
    } else {
        _containedParcels[cellIdx][_numContainedParcel[cellIdx]] = parcel;
    }
    parcel->hostCellIdx() = cellIdx;
    _numContainedParcel[cellIdx]++;
}
    
void MeshAdaptor::
resetContainedParcels() {
    for (int i = 0; i < mesh->totalNumGrid(CENTER, mesh->domain().numDim()); ++i) {
        _numContainedParcel[i] = 0;
    }
}

}