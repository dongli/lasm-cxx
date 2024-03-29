#ifndef __LASM_AdvectionManager__
#define __LASM_AdvectionManager__

#include "lasm_commons.h"
#include "ParcelManager.h"
#include "MeshAdaptor.h"

namespace lasm {

class AdvectionManager 
: public geomtk::AdvectionManagerInterface<2, Domain, Mesh, Field, VelocityField> {
    ParcelManager parcelManager;
    MeshAdaptor meshAdaptor;
    Regrid regrid;

    double maxFilament;
    double maxShapeBias;
    double radialMixWeight;
    double lateralMixWeight;
    double mixCoef;
    double reshapeCoef;
    double connectScale;

    double refVolume;       // Reference parcel volume.

    TreeType *gridTree;
    mat gridCoords;
public:
    AdvectionManager();
    virtual ~AdvectionManager();

    virtual void
    init(const Mesh &mesh);

    virtual void
    addTracer(const string &name, const string &unit, const string &comment);

    virtual void
    advance(double dt, const TimeLevelIndex<2> &newIdx,
            const VelocityField &velocityField);

    virtual void
    input(const TimeLevelIndex<2> &timeIdx, double *q);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx, int ncId) const;

    virtual double&
    density(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) {
        return meshAdaptor.density(timeIdx, tracerIdx, cellIdx);
    }

    virtual const Field<double, 2>&
    density(int tracerIdx) const {
        return meshAdaptor.density(tracerIdx);
    }

    virtual Field<double, 2>&
    density(int tracerIdx) {
        return meshAdaptor.density(tracerIdx);
    }

    virtual double&
    tendency(int tracerIdx, int cellIdx) {
        return meshAdaptor.tendency(tracerIdx, cellIdx);
    }

    virtual const Field<double>&
    tendency(int tracerIdx) const {
        return meshAdaptor.tendency(tracerIdx);
    }

    virtual Field<double>&
    tendency(int tracerIdx) {
        return meshAdaptor.tendency(tracerIdx);
    }

    Parcels&
    parcels() {
        return parcelManager.parcels();
    }
protected:
    void
    statistics(const TimeLevelIndex<2> &timeIdx);

    void
    integrate(double dt, const TimeLevelIndex<2> &newIdx,
              const VelocityField &velocityField);
    void
    connectParcelAndGrids(const TimeLevelIndex<2> &timeIdx,
                          Parcel *parcel);

    void
    connectParcelsAndGrids(const TimeLevelIndex<2> &timeIdx);

    void
    mixParcels(const TimeLevelIndex<2> &timeIdx);

    void
    remapDensityFromParcelsToGrids(const TimeLevelIndex<2> &timeIdx);

    void
    remapDensityFromGridsToParcels(const TimeLevelIndex<2> &timeIdx);

    void
    remapTendencyFromGridsToParcels(const TimeLevelIndex<2> &timeIdx);

    int
    getNeighborParcels(Parcel *parcel, Parcel **neighborParcels,
                       int maxNumNeighborParcel) const;
};

} // lasm

#endif // __GEOMTK_AdvectionManagerInterface__