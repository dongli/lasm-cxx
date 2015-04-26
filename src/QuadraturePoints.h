#ifndef __LASM_QuadraturePoints__
#define __LASM_QuadraturePoints__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class QuadraturePoints {
private:
    static const Domain *domain;
    static field<double> w;
    static field<BodyCoord> y;
    field<SpaceCoord> x;
    const Parcel *hostParcel;
public:
    QuadraturePoints(const Parcel *hostParcel);
    virtual ~QuadraturePoints();

    static void
    init(const Domain &domain);

    static int
    numPoint() {
        return y.size();
    }

    static const field<BodyCoord>&
    bodyCoords() {
        return y;
    }

    const field<double>&
    weights() const {
        return w;
    }

    double
    weight(int i, int j, int k = 0) const {
        return w(i, j, k);
    }

    double
    weight(int I) const {
        return w(I);
    }

    void
    updateSpaceCoords(const TimeLevelIndex<2> &timeIdx);

    const field<SpaceCoord>&
    spaceCoords() const {
        return x;
    }
}; // QuadraturePoints

} // lasm

#endif // __LASM_QuadraturePoints__