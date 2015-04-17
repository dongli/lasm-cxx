#ifndef __LASM_AdvectionManager_test__
#define __LASM_AdvectionManager_test__

#include "lasm.h"

using namespace lasm;

class AdvectionManagerTest : public ::testing::Test {
protected:
    Domain *domain;
    Mesh *mesh;
    ConfigManager configManager;
    AdvectionManager advectionManager;
    TimeLevelIndex<2> oldIdx, newIdx;

    virtual void
    SetUp() {
        domain = new Domain(3);
        domain->setAxis(0, "x", "x axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(1, "y", "y axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(2, "z", "z axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        mesh = new Mesh(*domain);
        mesh->init(10, 10, 10);
    }

    virtual void
    TearDown() {
        delete mesh;
        delete domain;
    }
};

TEST_F(AdvectionManagerTest, Init) {
    advectionManager.init(configManager, *mesh);
}

#endif // __LASM_AdvectionManager_test__