#ifndef __LASM_Parcel__
#define __LASM_Parcel__

#include "lasm_commons.h"

#include "MeshAdaptor.h"

namespace lasm {

class QuadraturePoints;
class SkeletonPoints;
class Tracers;

class Parcel {
protected:
    static const Mesh *mesh;
    int _id;
    int _rank;
    TimeLevels<SpaceCoord, 2> _x;
    TimeLevels<mat, 2> _H;
    TimeLevels<double, 2> _detH;
    TimeLevels<mat, 2> _invH;
    TimeLevels<vec, 2> _shapeSize;
    double _filament;
    mat _U, _V; vec _S;
    BodyCoord longAxisVertexY;
    SpaceCoord longAxisVertexX;
    QuadraturePoints *_quadraturePoints;
    SkeletonPoints *_skeletonPoints;
    Tracers *_tracers;

    int _hostCellIndex;
    TimeLevels<MeshIndex, 2> _meshIndex;
    uword _numConnectedCell;
    vector<int> _connectedCellIndexs;
public:
    Parcel();
    virtual ~Parcel();

    int
    id() const {
        return _id;
    }

    int
    rank() const {
        return _rank;
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

    double&
    volume(const TimeLevelIndex<2> &timeIdx) {
        return _detH.level(timeIdx);
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
    filament() const {
        return _filament;
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
    hostCellIndex() {
        return _hostCellIndex;
    }

    MeshIndex&
    meshIndex(const TimeLevelIndex<2> &timeIdx) {
        return _meshIndex.level(timeIdx);
    }

    static void
    init(const Mesh &mesh);

    virtual void
    init(int id, int rank);

    void
    updateDeformMatrix(const TimeLevelIndex<2> &timeIdx);

    void
    updateDeformMatrix(const TimeLevelIndex<2> &timeIdx, const vec &S);

    void
    updateDeformMatrix(const TimeLevelIndex<2> &timeIdx,
                       const mat &U, const vec &S, const mat &V);

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
    connectedCellIndexs() const {
        return _connectedCellIndexs;
    }

    uword
    numConnectedCell() const {
        return _numConnectedCell;
    }

    void
    dump(const TimeLevelIndex<2> &timeIdx, const MeshAdaptor &meshAdaptor,
         const char *fileName = NULL) const;
}; // Parcel

} // lasm

#endif // __LASM_Parcel__
