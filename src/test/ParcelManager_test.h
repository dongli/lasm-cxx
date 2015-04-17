#ifndef __LASM_ParcelManager_test__
#define __LASM_ParcelManager_test__

#include "lasm.h"

using namespace lasm;

class ParcelManagerTest : public ::testing::Test {
protected:
    Domain *domain;
    Mesh *mesh;
    ParcelManager parcelManager;
    TimeLevelIndex<2> oldIdx, newIdx;

    virtual void
    SetUp() {
        domain = new Domain(3);
        domain->setAxis(0, "x", "x axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(1, "y", "y axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        domain->setAxis(2, "z", "z axis", "m", 0, geomtk::PERIODIC, 1, geomtk::PERIODIC);
        mesh = new Mesh(*domain);
        mesh->init(10, 10, 10);
        ShapeFunction::init(*domain);
        QuadraturePoints::init(*domain);
        SkeletonPoints::init(*domain);
        Parcel::init(*domain);
        Tracers::add("q", "test tracer", "1");
    }

    virtual void
    TearDown() {
        delete mesh;
        delete domain;
    }
};

TEST_F(ParcelManagerTest, Init) {
    parcelManager.init(*mesh);
    // for (auto parcel : parcelManager.parcels()) {
    //     parcel->x(oldIdx).print();
    // }
}

#endif // __LASM_ParcelManager_test__