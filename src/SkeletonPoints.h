#ifndef __LASM_SkeletonPoints__
#define __LASM_SkeletonPoints__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class SkeletonPoints {
private:
    static const Domain *domain;
    static field<BodyCoord> y;
    const Parcel *hostParcel;
    TimeLevels<field<SpaceCoord>, 2> x, xl;
    TimeLevels<field<MeshIndex>, 2> I;
public:
    SkeletonPoints(const Parcel *hostParcel);
    virtual ~SkeletonPoints();

    static void
    init(const Domain &domain);

    void
    init();

    void
    init(const Mesh &mesh, const vec &sizes);

    static uword
    numPoint() {
        return y.size();
    }

    static const field<BodyCoord>&
    bodyCoords() {
        return y;
    }

    field<SpaceCoord>&
    spaceCoords(const TimeLevelIndex<2> &timeIdx) {
        return x.level(timeIdx);
    }

    field<MeshIndex>&
    meshIndexs(const TimeLevelIndex<2> &timeIdx) {
        return I.level(timeIdx);
    }

    const SpaceCoord&
    spaceCoord(const TimeLevelIndex<2> &timeIdx, int i) const {
        return x.level(timeIdx)(i);
    }

    SpaceCoord&
    spaceCoord(const TimeLevelIndex<2> &timeIdx, int i) {
        return x.level(timeIdx)(i);
    }

    const field<SpaceCoord>&
    localSpaceCoords(const TimeLevelIndex<2> &timeIdx) const {
        return xl.level(timeIdx);
    }

    const SpaceCoord&
    localSpaceCoord(const TimeLevelIndex<2> &timeIdx, int i) const {
        return xl.level(timeIdx)(i);
    }

    void
    updateLocalSpaceCoords(const TimeLevelIndex<2> &timeIdx);
}; // SkeletonPoints

} // lasm

#endif // __LASM_SkeletonPoints__