#ifndef __LASM_MeshAdaptor_test__
#define __LASM_MeshAdaptor_test__

#include "lasm.h"

using namespace lasm;

class MeshAdaptorTest : public ::testing::Test {
protected:
    Domain *domain;
    Mesh *mesh;
    MeshAdaptor meshAdaptor;
    TimeLevelIndex<2> oldIdx;

    virtual void
    SetUp() {
        domain = new Domain(3);
        domain->setAxis(0, "x", "x axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(1, "y", "y axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(2, "z", "z axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        mesh = new Mesh(*domain);
        mesh->init(10, 10, 10);
        meshAdaptor.init(*mesh);
    }

    virtual void
    TearDown() {

    }
};

TEST_F(MeshAdaptorTest, Basics) {
    meshAdaptor.addTracer("q", "1", "1");
    ASSERT_EQ(1, meshAdaptor.masses().size());
}

#endif // __LASM_MeshAdaptor_test__