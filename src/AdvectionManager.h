#ifndef __LASM_AdvectionManager__
#define __LASM_AdvectionManager__

#include "lasm_commons.h"
#include "ParcelManager.h"
#include "MeshAdaptor.h"

namespace lasm {

class AdvectionManager 
: public geomtk::AdvectionManagerInterface<2, Domain, Mesh, VelocityField> {
    ParcelManager parcelManager;
    MeshAdaptor meshAdaptor;
    Regrid regrid;

    double filamentLimit;   // control the tracer filament degree
    double minBiasLimit;    //
    double maxBiasLimit;    //
    double radialMixing;    // control the radial mixing degree
    double lateralMixing;   // control the lateral mixing degree
    double restoreFactor;   // control the base density restore degree
public:
    AdvectionManager();
    virtual ~AdvectionManager();

    virtual void
    init(const ConfigManager &configManager, const Mesh &mesh);

    virtual void
    addTracer(const string &name, const string &unit, const string &comment);

    virtual void
    advance(double dt, const TimeLevelIndex<2> &newIdx,
            const VelocityField &velocityField);

    virtual void
    input(const TimeLevelIndex<2> &timeIdx, double *q);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx, int ncId);

    virtual double
    density(const TimeLevelIndex<2> &timeIdx, int tracerIdx, int cellIdx) const {
        return meshAdaptor.mass(timeIdx, tracerIdx, cellIdx)/meshAdaptor.volume(cellIdx);
    }
protected:
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
    remapFromParcelsToGrids(const TimeLevelIndex<2> &timeIdx);

    void
    remapFromGridsToParcels(const TimeLevelIndex<2> &timeIdx);

    vector<Parcel*>
    getNeighborParcels(Parcel *parcel) const;
};

} // lasm

#endif // __GEOMTK_AdvectionManagerInterface__