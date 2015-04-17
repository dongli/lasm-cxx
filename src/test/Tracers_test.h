#ifndef __LASM_TracersTest__
#define __LASM_TracersTest__

#include "lasm.h"

using namespace lasm;

class TracersTest : public ::testing::Test {
protected:
    Domain *domain;
    Mesh *mesh;
    Parcel parcel;
    Tracers *tracers;
    TimeLevelIndex<2> oldIdx;

    virtual void
    SetUp() {
        tracers = new Tracers(NULL);
    }

    virtual void
    TearDown() {
    }
};

TEST_F(TracersTest, Basics) {
    Tracers::add("qm", "moisture", "kg m-3");
    ASSERT_DEATH(Tracers::add("qm", "...", "..."), "has been added");
    tracers->init(); // There is nothing happens in init().
    tracers->add();
    ASSERT_EQ(1, tracers->masses(oldIdx).size());
}

#endif // __LASM_TracersTest__