#ifndef __LASM_BarotropicTest__
#define __LASM_BarotropicTest__

#include "DataTest.h"

class BarotropicTest : public DataTest {
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, Field, VelocityField> AdvectionManager;
    typedef geomtk::AdvectionTestInterface<2, Domain, Mesh, Field, VelocityField, IOManager> Interface;

    BarotropicTest();
    virtual ~BarotropicTest();

    virtual void
    init(AdvectionManager &advectionManager);

    virtual void
    setInitialCondition(AdvectionManager &advectionManager);
};

#endif // __LASM_BarotropicTest__
