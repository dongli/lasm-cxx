#include "SkeletonPoints.h"
#include "Parcel.h"

namespace lasm {

const Domain *SkeletonPoints::domain;
field<BodyCoord> SkeletonPoints::y;

SkeletonPoints::
SkeletonPoints(const Parcel *hostParcel) {
    this->hostParcel = hostParcel;
} // SkeletonPoints

SkeletonPoints::
~SkeletonPoints() {
} // ~SkeletonPoints

void SkeletonPoints::
init(const Domain &_domain) {
    domain = &_domain;
    // Set the fixed body coordinates of the skeleton points.
    y.set_size(domain->numDim()*2);
    for (uword i = 0; i < y.size(); ++i) {
        y[i].init(domain->numDim());
    }
    double d = 1;
    if (domain->numDim() == 2) {
        y[0]() <<  -d << 0.0 << arma::endr;
        y[1]() << 0.0 <<  -d << arma::endr;
        y[2]() <<   d << 0.0 << arma::endr;
        y[3]() << 0.0 <<   d << arma::endr;
    } else if (domain->numDim() == 3) {
        y[0]() <<  -d << 0.0 << 0.0 << arma::endr;
        y[1]() << 0.0 <<  -d << 0.0 << arma::endr;
        y[2]() <<   d << 0.0 << 0.0 << arma::endr;
        y[3]() << 0.0 <<   d << 0.0 << arma::endr;
        y[4]() << 0.0 << 0.0 <<  -d << arma::endr;
        y[5]() << 0.0 << 0.0 <<   d << arma::endr;
    }
} // init

void SkeletonPoints::
init(const Mesh &mesh, const vec &sizes) {
    TimeLevelIndex<2> timeIdx;
#ifndef NDEBUG
    assert(sizes.size() == domain->numDim());
#endif
    // Set the space coordinates of the skeleton points.
    for (int l = 0; l < x.numLevel(); ++l) {
        x.level(l).set_size(y.size());
        xl.level(l).set_size(y.size());
        I.level(l).set_size(y.size());
        for (uword i = 0; i < x.level(l).size(); ++i) {
            x.level(l)[i].init(domain->numDim());
            xl.level(l)[i].init(domain->numDim());
            I.level(l)[i].init(domain->numDim());
        }
    }
#if defined LASM_IN_CARTESIAN
    if (domain->numDim() == 2) {
        x.level(timeIdx)[0].set(x0(0)-sizes(0), x0(1));
        x.level(timeIdx)[1].set(x0(0), x0(1)-sizes(1));
        x.level(timeIdx)[2].set(x0(0)+sizes(0), x0(1));
        x.level(timeIdx)[3].set(x0(0), x0(1)+sizes(1));
    } else if (domain->numDim() == 3) {
        x.level(timeIdx)[0].set(x0(0)-sizes(0), x0(1), x0(2));
        x.level(timeIdx)[1].set(x0(0), x0(1)-sizes(1), x0(2));
        x.level(timeIdx)[2].set(x0(0)+sizes(0), x0(1), x0(2));
        x.level(timeIdx)[3].set(x0(0), x0(1)+sizes(1), x0(2));
        x.level(timeIdx)[4].set(x0(0), x0(1), x0(2)-sizes(2));
        x.level(timeIdx)[5].set(x0(0), x0(1), x0(2)+sizes(2));
    }
#elif defined LASM_IN_SPHERE
    const SpaceCoord &x0 = hostParcel->x(timeIdx);
    if (domain->numDim() == 2) {
        domain->rotateBack(x0, x.level(timeIdx)[0], 0, PI*0.5-sizes(1)/domain->radius());
        domain->rotateBack(x0, x.level(timeIdx)[1], PI*0.5, PI*0.5-sizes(0)/domain->radius());
        domain->rotateBack(x0, x.level(timeIdx)[2], PI, PI*0.5-sizes(1)/domain->radius());
        domain->rotateBack(x0, x.level(timeIdx)[3], PI*1.5, PI*0.5-sizes(0)/domain->radius());
    } else if (domain->numDim() == 3) {
        x.level(timeIdx)[0](2) = x0(2);
        domain->rotateBack(x0, x.level(timeIdx)[0], 0, PI*0.5-sizes(1)/domain->radius());
        x.level(timeIdx)[1](2) = x0(2);
        domain->rotateBack(x0, x.level(timeIdx)[1], PI*0.5, PI*0.5-sizes(0)/domain->radius());
        x.level(timeIdx)[2](2) = x0(2);
        domain->rotateBack(x0, x.level(timeIdx)[2], PI, PI*0.5-sizes(1)/domain->radius());
        x.level(timeIdx)[3](2) = x0(2);
        domain->rotateBack(x0, x.level(timeIdx)[3], PI*1.5, PI*0.5-sizes(0)/domain->radius());
        x.level(timeIdx)[4](0) = x0(0);
        x.level(timeIdx)[4](1) = x0(1);
        x.level(timeIdx)[4](2) = x0(2)-sizes(2);
        x.level(timeIdx)[5](0) = x0(0);
        x.level(timeIdx)[5](1) = x0(1);
        x.level(timeIdx)[5](2) = x0(2)+sizes(2);
    }
#endif
    for (uword i = 0; i < I.level(timeIdx).size(); ++i) {
        if (domain->isValid(x.level(timeIdx)[i])) {
            I.level(timeIdx)[i].locate(mesh, x.level(timeIdx)[i]);
        } else {
            I.level(timeIdx)[i].reset();
        }
    }
    updateLocalSpaceCoords(timeIdx);
} // init

void SkeletonPoints::
updateLocalSpaceCoords(const TimeLevelIndex<2> &timeIdx) {
    for (uword i = 0; i < I.level(timeIdx).size(); ++i) {
        if (domain->isValid(x.level(timeIdx)[i])) {
#if defined LASM_IN_CARTESIAN
            xl.level(timeIdx)[i]() = domain->diffCoord(x.level(timeIdx)[i],
                                                       hostParcel->x(timeIdx));
#elif defined LASM_IN_SPHERE
            domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                           hostParcel->x(timeIdx), x.level(timeIdx)[i],
                           xl.level(timeIdx)[i]());
#endif
        }
    }
}

} // lasm