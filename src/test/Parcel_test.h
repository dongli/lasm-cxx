#ifndef __LASM_Parcel_test__
#define __LASM_Parcel_test__

#include "lasm.h"

using namespace lasm;

class ParcelTest : public ::testing::Test {
protected:
    Domain *domain;
    Mesh *mesh;
    Parcel parcel;
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
        parcel.init(0);
        parcel.x(oldIdx).setCoord(0.5, 0.5, 0.5);
        parcel.meshIdx(oldIdx).locate(*mesh, parcel.x(oldIdx));
    }

    virtual void
    TearDown() {
        delete mesh;
        delete domain;
    }
};

TEST_F(ParcelTest, Basics) {
    parcel.x(oldIdx);
    ASSERT_EQ(parcel.id(), 0);
    ASSERT_EQ(parcel.x(oldIdx)(0), 0.5);
    ASSERT_EQ(parcel.x(oldIdx)(1), 0.5);
    ASSERT_EQ(parcel.x(oldIdx)(2), 0.5);
    ASSERT_EQ(parcel.H(oldIdx).n_rows, 3);
    ASSERT_EQ(parcel.H(oldIdx).n_cols, 3);
    ASSERT_EQ(QuadraturePoints::bodyCoords().n_rows, 5);
    ASSERT_EQ(QuadraturePoints::bodyCoords().n_cols, 5);
    ASSERT_EQ(QuadraturePoints::bodyCoords().n_slices, 5);
    ASSERT_EQ(parcel.quadraturePoints().spaceCoords().n_rows, 5);
    ASSERT_EQ(parcel.quadraturePoints().spaceCoords().n_cols, 5);
    ASSERT_EQ(parcel.quadraturePoints().spaceCoords().n_slices, 5);
}

TEST_F(ParcelTest, SkeletonPoints) {
    vec sizes(3);
    sizes(0) = 0.12;
    sizes(1) = 0.12;
    sizes(2) = 0.12;
    parcel.skeletonPoints().init(*mesh, sizes);
    const field<SpaceCoord> &x = parcel.skeletonPoints().spaceCoords(oldIdx);
    const field<SpaceCoord> &xl = parcel.skeletonPoints().localSpaceCoords(oldIdx);
    ASSERT_EQ(0.38, x(0)(0)); ASSERT_EQ(-0.12, xl(0)(0));
    ASSERT_EQ(0.5,  x(0)(1)); ASSERT_EQ(  0.0, xl(0)(1));
    ASSERT_EQ(0.5,  x(0)(2)); ASSERT_EQ(  0.0, xl(0)(2));
    ASSERT_EQ(0.5,  x(1)(0)); ASSERT_EQ(  0.0, xl(1)(0));
    ASSERT_EQ(0.38, x(1)(1)); ASSERT_EQ(-0.12, xl(1)(1));
    ASSERT_EQ(0.5,  x(1)(2)); ASSERT_EQ(  0.0, xl(1)(2));
    ASSERT_EQ(0.62, x(2)(0)); ASSERT_EQ( 0.12, xl(2)(0));
    ASSERT_EQ(0.5,  x(2)(1)); ASSERT_EQ(  0.0, xl(2)(1));
    ASSERT_EQ(0.5,  x(2)(2)); ASSERT_EQ(  0.0, xl(2)(2));
    ASSERT_EQ(0.5,  x(3)(0)); ASSERT_EQ(  0.0, xl(3)(0));
    ASSERT_EQ(0.62, x(3)(1)); ASSERT_EQ( 0.12, xl(3)(1));
    ASSERT_EQ(0.5,  x(3)(2)); ASSERT_EQ(  0.0, xl(3)(2));
    ASSERT_EQ(0.5,  x(4)(0)); ASSERT_EQ(  0.0, xl(4)(0));
    ASSERT_EQ(0.5,  x(4)(1)); ASSERT_EQ(  0.0, xl(4)(1));
    ASSERT_EQ(0.38, x(4)(2)); ASSERT_EQ(-0.12, xl(4)(2));
    ASSERT_EQ(0.5,  x(5)(0)); ASSERT_EQ(  0.0, xl(5)(0));
    ASSERT_EQ(0.5,  x(5)(1)); ASSERT_EQ(  0.0, xl(5)(1));
    ASSERT_EQ(0.62, x(5)(2)); ASSERT_EQ( 0.12, xl(5)(2));
}

TEST_F(ParcelTest, DeformMatrix) {
    vec sizes(3);
    sizes(0) = 0.12;
    sizes(1) = 0.12;
    sizes(2) = 0.12;
    parcel.skeletonPoints().init(*mesh, sizes);
    parcel.updateDeformMatrix(oldIdx);
    const mat &H = parcel.H(oldIdx);
    ASSERT_EQ(H(0, 0), 0.12); ASSERT_EQ(H(0, 1),  0.0); ASSERT_EQ(H(0, 2),  0.0);
    ASSERT_EQ(H(1, 0),  0.0); ASSERT_EQ(H(1, 1), 0.12); ASSERT_EQ(H(1, 2),  0.0);
    ASSERT_EQ(H(2, 0),  0.0); ASSERT_EQ(H(2, 1),  0.0); ASSERT_EQ(H(2, 2), 0.12);
    const vec &S = parcel.S();
    ASSERT_EQ(S(0), 0.12); ASSERT_EQ(S(1), 0.12); ASSERT_EQ(S(2), 0.12);
    ASSERT_NEAR(pow(0.12, 3), parcel.volume(oldIdx), 1.0e-12);
    vec v = domain->diffCoord(parcel.x(oldIdx), parcel.longAxisVertexSpaceCoord());
    ASSERT_NEAR(0.12, norm(v), 1.0e-12);
    ASSERT_NEAR(1, parcel.filamentDegree(oldIdx), 1.0e-12);
}

TEST_F(ParcelTest, QuadraturePoints) {
    vec sizes(3);
    sizes(0) = 0.12;
    sizes(1) = 0.12;
    sizes(2) = 0.12;
    parcel.skeletonPoints().init(*mesh, sizes);
    parcel.updateDeformMatrix(oldIdx);
    parcel.quadraturePoints().updateSpaceCoords(oldIdx);
    const field<double> &weights = parcel.quadraturePoints().weights();
    double sum = 0;
    for (int k = 0; k < weights.n_slices; ++k) {
        for (int j = 0; j < weights.n_cols; ++j) {
            for (int i = 0; i < weights.n_rows; ++i) {
                sum += weights(i, j, k);
            }
        }
    }
    ASSERT_NEAR(1, sum, 1.0e-12);
}

#endif // __LASM_Parcel_test__