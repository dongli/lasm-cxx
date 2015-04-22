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
    for (int i = 0; i < y.size(); ++i) {
        y[i].setNumDim(domain->numDim());
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
    const SpaceCoord &x0 = hostParcel->x(timeIdx);
    assert(sizes.size() == domain->numDim());
    // Set the space coordinates of the skeleton points.
    for (int l = 0; l < x.numLevel(); ++l) {
        x.level(l).set_size(y.size());
        xl.level(l).set_size(y.size());
        _meshIdx.level(l).set_size(y.size());
        for (int i = 0; i < x.level(l).size(); ++i) {
            x.level(l)[i].setNumDim(domain->numDim());
            xl.level(l)[i].setNumDim(domain->numDim());
            _meshIdx.level(l)[i].setNumDim(domain->numDim());
        }
    }
#if defined LASM_IN_CARTESIAN
    if (domain->numDim() == 2) {
        x.level(timeIdx)[0].setCoord(x0(0)-sizes(0), x0(1));
        x.level(timeIdx)[1].setCoord(x0(0), x0(1)-sizes(1));
        x.level(timeIdx)[2].setCoord(x0(0)+sizes(0), x0(1));
        x.level(timeIdx)[3].setCoord(x0(0), x0(1)+sizes(1));
    } else if (domain->numDim() == 3) {
        x.level(timeIdx)[0].setCoord(x0(0)-sizes(0), x0(1), x0(2));
        x.level(timeIdx)[1].setCoord(x0(0), x0(1)-sizes(1), x0(2));
        x.level(timeIdx)[2].setCoord(x0(0)+sizes(0), x0(1), x0(2));
        x.level(timeIdx)[3].setCoord(x0(0), x0(1)+sizes(1), x0(2));
        x.level(timeIdx)[4].setCoord(x0(0), x0(1), x0(2)-sizes(2));
        x.level(timeIdx)[5].setCoord(x0(0), x0(1), x0(2)+sizes(2));
    }
#endif
    for (int i = 0; i < _meshIdx.level(timeIdx).size(); ++i) {
        if (domain->isValid(x.level(timeIdx)[i])) {
#if defined LASM_IN_CARTESIAN
            xl.level(timeIdx)[i]() = domain->diffCoord(x.level(timeIdx)[i], hostParcel->x(timeIdx));
#endif
            _meshIdx.level(timeIdx)[i].locate(mesh, x.level(timeIdx)[i]);
        } else {
            _meshIdx.level(timeIdx)[i].reset();
        }
    }
}

} // lasm