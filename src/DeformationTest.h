#ifndef __LASM_DeformationTest__
#define __LASM_DeformationTest__

#include "lasm_commons.h"

// TODO: Move this class to GEOMTK.

class DeformationTest
: public geomtk::AdvectionTestInterface<2, Domain, Mesh, Field, VelocityField, IOManager> {
protected:
    std::string subcase;
    double period;
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, Field, VelocityField> AdvectionManager;
    typedef geomtk::AdvectionTestInterface<2, Domain, Mesh, Field, VelocityField, IOManager> Interface;

    DeformationTest();
    virtual ~DeformationTest();

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
}; // DeformationTest

#endif // __LASM_DeformationTest__
