#ifndef __LASM_DeformationTest__
#define __LASM_DeformationTest__

#include "lasm_commons.h"

// TODO: Move this class to GEOMTK.

class DeformationTest
: public geomtk::AdvectionTestInterface<2, Domain, Mesh, Field<double>, VelocityField, IOManager> {
protected:
    std::string subcase;
    double period;
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, VelocityField> AdvectionManager;
    typedef geomtk::AdvectionTestInterface<2, Domain, Mesh, Field<double>, VelocityField, IOManager> Interface;

    DeformationTest();
    virtual ~DeformationTest();

    virtual void
    init(const ConfigManager &configManager, AdvectionManager &advectionManager);

    virtual void
    setInitialCondition(AdvectionManager &advectionManager);

    virtual void
    advanceDynamics(const TimeLevelIndex<2> &timeIdx,
                    AdvectionManager &advectionManager);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx,
           AdvectionManager &advectionManager);
}; // DeformationTest

#endif // __LASM_DeformationTest__
