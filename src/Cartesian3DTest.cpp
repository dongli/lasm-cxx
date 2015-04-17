#include "Cartesian3DTest.h"

namespace lasm {

Cartesian3DTest::
Cartesian3DTest() {

}

Cartesian3DTest::
~Cartesian3DTest() {

}

void Cartesian3DTest::
init(const ConfigManager &configManager, TimeManager &timeManager) {
    _domain = new Domain(3);
    _mesh = new Mesh(*_domain);
    // Read in WRF-LES 3D mesh.

    // Initialize velocity field.
    velocityField.create(*_mesh, true, true);
} // init

void Cartesian3DTest::
setInitialCondition(AdvectionManager &advectionManager) {
    // Read in WRF-LES 3D tracer.
} // setInitialCondition

void Cartesian3DTest::
advanceDynamics(AdvectionManager &advectionManager) {
    // Read in WRF-LES 3D flow.
} // advanceDynamics

void Cartesian3DTest::
output(const TimeLevelIndex<2> &timeIdx, AdvectionManager &advectionManager) {

} // output

}