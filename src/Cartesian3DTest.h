#ifndef __LASM_Cartesian3DTest__
#define __LASM_Cartesian3DTest__

#include "lasm_commons.h"

// TODO: Move this class to GEOMTK.

class Cartesian3DTest
: public geomtk::AdvectionTestInterface<2, Domain, Mesh, Field<double>, VelocityField, IOManager> {
    arma::uword dataIdx;
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, VelocityField> AdvectionManager;
    typedef geomtk::AdvectionTestInterface<2, Domain, Mesh, Field<double>, VelocityField, IOManager> Interface;

    Cartesian3DTest();
    virtual ~Cartesian3DTest();

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
}; // Cartesian3DTest

#endif // __LASM_Cartesian3DTest__