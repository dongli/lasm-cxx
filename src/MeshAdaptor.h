#ifndef __LASM_MeshAdaptor__
#define __LASM_MeshAdaptor__

#include "lasm_commons.h"

namespace lasm {

class MeshAdaptor {
    const Mesh *mesh;
    vector<Field<double, 2>*> _densities;
    vector<Field<double>*> _tendencies;
    vector<uword> _numConnectedParcel;
    vector<vector<Parcel*> > _connectedParcels;
    vector<vector<double> > _remapWeights;
    vector<uword> _numContainedParcel;
    vector<vector<Parcel*> > _containedParcels;
    vector<uword> _cellIndexs;
public:
    MeshAdaptor();
    virtual ~MeshAdaptor();

    void
    init(const Mesh &mesh);

    uword
    cellIndex(int i) const {
        return _cellIndexs[i];
    }

    const SpaceCoord&
    coord(int cellIdx) const {
        return mesh->gridCoord(CENTER, cellIdx);
    }

    double
    volume(int cellIdx) const {
        return mesh->cellVolume(cellIdx);
    }

    const Field<double, 2>&
    density(int tracerIdx) const {
        return *_densities[tracerIdx];
    }

    Field<double, 2>&
    density(int tracerIdx) {
        return *_densities[tracerIdx];
    }

    double
    density(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) const {
        return (*_densities[tracerIdx])(timeIdx, cellIdx);
    }

    double&
    density(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) {
        return (*_densities[tracerIdx])(timeIdx, cellIdx);
    }

    double
    mass(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) const {
        return (*_densities[tracerIdx])(timeIdx, cellIdx)*volume(cellIdx);
    }

    double&
    tendency(int tracerIdx, int cellIdx) {
        return (*_tendencies[tracerIdx])(cellIdx);
    }

    const Field<double>&
    tendency(int tracerIdx) const {
        return *_tendencies[tracerIdx];
    }

    Field<double>&
    tendency(int tracerIdx) {
        return *_tendencies[tracerIdx];
    }

    void
    addTracer(const string &name, const string &unit, const string &comment);

    void
    resetTracers(const TimeLevelIndex<2> &timeIdx);

    void
    connectParcel(int cellIdx, Parcel *parcel, double weight);

    double
    remapWeight(int cellIdx, Parcel *parcel) const;

    double&
    remapWeight(int cellIdx, Parcel *parcel);

    double
    remapWeight(int cellIdx, int parcelIdx) const {
        return _remapWeights[cellIdx][parcelIdx];
    }

    uword
    numConnectedParcel(int cellIdx) const {
        return _numConnectedParcel[cellIdx];
    }

    const vector<Parcel*>& connectedParcels(int cellIdx) const {
        return _connectedParcels[cellIdx];
    }

    void
    containParcel(int cellIdx, Parcel *parcel);

    uword
    numContainedParcel(int cellIdx) const {
        return _numContainedParcel[cellIdx];
    }

    const vector<Parcel*>&
    containedParcels(int cellIdx) const {
        return _containedParcels[cellIdx];
    }

    void
    resetConnectedAndContainedParcels();
}; // MeshAdaptor

} // lasm

#endif // __LASM_MeshAdaptor__