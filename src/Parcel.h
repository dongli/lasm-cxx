#ifndef __LASM_Parcel__
#define __LASM_Parcel__

#include "lasm_commons.h"

namespace lasm {

class QuadraturePoints;
class SkeletonPoints;
class Tracers;

class Parcel {
protected:
    static const Domain *domain;
    int _ID;
    TimeLevels<SpaceCoord, 2> _x;
    TimeLevels<mat, 2> _H;
    TimeLevels<double, 2> _detH;
    TimeLevels<mat, 2> _invH;
    TimeLevels<vec, 2> _shapeSize;
    TimeLevels<double, 2> _filamentDegree;
    mat _U, _V; vec _S;
    BodyCoord longAxisVertexY;
    SpaceCoord longAxisVertexX;
    QuadraturePoints *_quadraturePoints;
    SkeletonPoints *_skeletonPoints;
    Tracers *_tracers;

    int _hostCellIdx;
    TimeLevels<MeshIndex, 2> _meshIdx;
    int _numConnectedCell;
    vector<int> _connectedCellIdxs;
public:
    Parcel();
    virtual ~Parcel();

    int
    ID() const {
        return _ID;
    }

    const SpaceCoord&
    x(const TimeLevelIndex<2> &timeIdx) const {
        return _x.level(timeIdx);
    }

    SpaceCoord&
    x(const TimeLevelIndex<2> &timeIdx) {
        return _x.level(timeIdx);
    }

    const mat&
    H(const TimeLevelIndex<2> &timeIdx) const {
        return _H.level(timeIdx);
    }

    mat&
    H(const TimeLevelIndex<2> &timeIdx) {
        return _H.level(timeIdx);
    }

    double
    detH(const TimeLevelIndex<2> &timeIdx) const {
        return _detH.level(timeIdx);
    }

    double
    volume(const TimeLevelIndex<2> &timeIdx) const {
        return detH(timeIdx);
    }

    const mat&
    invH(const TimeLevelIndex<2> &timeIdx) const {
        return _invH.level(timeIdx);
    }

    const vec&
    shapeSize(const TimeLevelIndex<2> &timeIdx) const {
        return _shapeSize.level(timeIdx);
    }

    double
    filamentDegree(const TimeLevelIndex<2> &timeIdx) const {
        return _filamentDegree.level(timeIdx);
    }

    const mat&
    U() const {
        return _U;
    }

    const mat&
    V() const {
        return _V;
    }

    const vec&
    S() const {
        return _S;
    }

    const SpaceCoord&
    longAxisVertexSpaceCoord() const {
        return longAxisVertexX;
    }

    QuadraturePoints&
    quadraturePoints() {
        return *_quadraturePoints;
    }

    SkeletonPoints&
    skeletonPoints() {
        return *_skeletonPoints;
    }

    Tracers&
    tracers() {
        return *_tracers;
    }

    int&
    hostCellIdx() {
        return _hostCellIdx;
    }

    MeshIndex&
    meshIdx(const TimeLevelIndex<2> &timeIdx) {
        return _meshIdx.level(timeIdx);
    }

    static void
    init(const Domain &domain);

    virtual void
    init(int ID);

    void
    updateDeformMatrix(const TimeLevelIndex<2> &timeIdx);

    void
    updateDeformMatrix(const TimeLevelIndex<2> &timeIdx, const vec &S);

    void
    resetSkeletonPoints(const TimeLevelIndex<2> &timeIdx,
                        const Mesh &mesh);

    void
    updateVolume(const TimeLevelIndex<2> &timeIdx, double volume);

    void
    updateShapeSize(const TimeLevelIndex<2> &timeIdx);

    void
    calcSpaceCoord(const TimeLevelIndex<2> &timeIdx,
                   const BodyCoord &y, SpaceCoord &x) const;

    void
    calcBodyCoord(const TimeLevelIndex<2> &timeIdx,
                  const SpaceCoord &x, BodyCoord &y) const;

    double
    shapeFunction(const TimeLevelIndex<2> &timeIdx,
                  const BodyCoord &y) const;

    void
    connectCell(int cellIdx);

    void
    resetConnectedCells();

    const vector<int>&
    connectedCellIdxs() const {
        return _connectedCellIdxs;
    }

    int
    numConnectedCell() const {
        return _numConnectedCell;
    }
}; // Parcel

} // lasm

#endif // __LASM_Parcel__
