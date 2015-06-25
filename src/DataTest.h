#ifndef __LASM_DataTest__
#define __LASM_DataTest__

#include "lasm_commons.h"

// TODO: Move this class to GEOMTK.

class DataTest
: public geomtk::AdvectionTestInterface<2, Domain, Mesh, Field, VelocityField, IOManager> {
    arma::uword dataIdx;
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, Field, VelocityField> AdvectionManager;
    typedef geomtk::AdvectionTestInterface<2, Domain, Mesh, Field, VelocityField, IOManager> Interface;

    DataTest();
    virtual ~DataTest();

    virtual void
    init(AdvectionManager &advectionManager);

    virtual void
    setInitialCondition(AdvectionManager &advectionManager);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx,
           const AdvectionManager &advectionManager);
protected:
    virtual void
    setVelocityField(const TimeLevelIndex<2> &timeIdx);
}; // DataTest

#endif // __LASM_DataTest__