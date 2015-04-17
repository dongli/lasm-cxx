#ifndef __LASM_Cartesian3DTest__
#define __LASM_Cartesian3DTest__

#include "lasm_commons.h"

namespace lasm {

class Cartesian3DTest
: public geomtk::AdvectionTestInterface<2, Domain, Mesh, IOManager, VelocityField> {
public:
    typedef geomtk::AdvectionManagerInterface<2, Domain, Mesh, VelocityField> AdvectionManager;

    Cartesian3DTest();
    virtual ~Cartesian3DTest();

    virtual void
    init(const ConfigManager &configManager, TimeManager &timeManager);

    virtual void
    setInitialCondition(AdvectionManager &advectionManager);

    virtual void
    advanceDynamics(AdvectionManager &advectionManager);

    virtual void
    output(const TimeLevelIndex<2> &timeIdx,
           AdvectionManager &advectionManager);
};

}

#endif // __LASM_Cartesian3DTest__