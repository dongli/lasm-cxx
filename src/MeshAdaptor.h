#ifndef __LASM_MeshAdaptor__
#define __LASM_MeshAdaptor__

#include "lasm_commons.h"

namespace lasm {

class MeshAdaptor {
    const Mesh *mesh;
    vector<Field<double, 2>*> _masses;
    vector<int> _numConnectedParcel;
    vector<vector<Parcel*> > _connectedParcels;
    vector<vector<double> > _remapWeights;
    vector<int> _numContainedParcel;
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
    mass(const TimeLevelIndex<2> &timeIdx, int speciesIdx, int cellIdx) {
        return (*_masses[speciesIdx])(timeIdx, cellIdx);
    }

    void
    addTracer(const string &name, const string &unit, const string &comment);

    void
    resetTracers(const TimeLevelIndex<2> &timeIdx);

    void
    connectParcel(int cellIdx, Parcel *parcel, double weight);

    void
    resetConnectedParcels();

    double
    remapWeight(int cellIdx, Parcel *parcel) const;

    int
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

    int
    numContainedParcel(int cellIdx) const {
        return _numContainedParcel[cellIdx];
    }

    const vector<Parcel*>& containedParcels(int cellIdx) const {
        return _containedParcels[cellIdx];
    }
};

}

#endif // __LASM_MeshAdaptor__