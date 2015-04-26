#ifndef __LASM_MeshAdaptor__
#define __LASM_MeshAdaptor__

#include "lasm_commons.h"

namespace lasm {

class MeshAdaptor {
    const Mesh *mesh;
    vector<Field<double, 2>*> _masses;
    vector<uword> _numConnectedParcel;
    vector<vector<Parcel*> > _connectedParcels;
    vector<vector<uword> > _numContainedQuadraturePoint;
    vector<vector<double> > _remapWeights;
    vector<uword> _numContainedParcel;
    vector<vector<Parcel*> > _containedParcels;
public:
    MeshAdaptor();
    virtual ~MeshAdaptor();

    void
    init(const Mesh &mesh);

    const SpaceCoord&
    coord(int cellIdx) const {
        return mesh->gridCoord(CENTER, cellIdx);
    }

    double
    volume(int cellIdx) const {
        return mesh->cellVolume(cellIdx);
    }

    const vector<Field<double, 2>*>&
    masses() const {
        return _masses;
    }

    double&
    mass(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) const {
        return (*_masses[tracerIdx])(timeIdx, cellIdx);
    }

    void
    addTracer(const string &name, const string &unit, const string &comment);

    void
    resetTracers(const TimeLevelIndex<2> &timeIdx);

    void
    connectParcel(int cellIdx, Parcel *parcel, double weight);

    void
    resetConnectedParcels();

    uword
    numContainedQuadraturePoint(int cellIdx, Parcel *parcel) const;

    double
    remapWeight(int cellIdx, Parcel *parcel) const;

    double&
    remapWeight(int cellIdx, Parcel *parcel);

    uword
    numConnectedParcel(int cellIdx) const {
        return _numConnectedParcel[cellIdx];
    }

    const vector<Parcel*>& connectedParcels(int cellIdx) const {
        return _connectedParcels[cellIdx];
    }

    void
    containParcel(int cellIdx, Parcel *parcel);

    void
    resetContainedParcels();

    uword
    numContainedParcel(int cellIdx) const {
        return _numContainedParcel[cellIdx];
    }

    const vector<Parcel*>&
    containedParcels(int cellIdx) const {
        return _containedParcels[cellIdx];
    }
}; // MeshAdaptor

} // lasm

#endif // __LASM_MeshAdaptor__