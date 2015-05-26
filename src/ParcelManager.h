#ifndef __LASM_ParcelManager__
#define __LASM_ParcelManager__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class ParcelManager {
    const Mesh *mesh;
    Parcels _parcels;
public:
    ParcelManager();
    virtual ~ParcelManager();

    Parcels&
    parcels() {
        return _parcels;
    }

    void
    init(const Mesh &mesh);

    void
    output(const TimeLevelIndex<2> &timeIdx, int ncId) const;

    void
    addTracer();

    void
    resetDensities();

    void
    resetTendencies();
}; // ParcelManager

} // lasm

#endif // __LASM_ParcelManager__