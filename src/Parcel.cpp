#include "Parcel.h"
#include "QuadraturePoints.h"
#include "SkeletonPoints.h"
#include "Tracers.h"
#include "ShapeFunction.h"

namespace lasm {

const Domain *Parcel::domain = NULL;

Parcel::
Parcel() {
    _quadraturePoints = NULL;
    _skeletonPoints = NULL;
    _tracers = NULL;
}

Parcel::
~Parcel() {
    if (_quadraturePoints != NULL) delete _quadraturePoints;
    if (_skeletonPoints != NULL) delete _skeletonPoints;
    if (_tracers != NULL) delete _tracers;
}

void Parcel::
init(const Domain &_domain) {
    domain = &_domain;
} // init

void Parcel::
init(int id) {
    _id = id;
    for (int l = 0; l < 2; ++l) {
        _x.level(l).init(domain->numDim());
        _H.level(l).set_size(domain->numDim(), domain->numDim());
        _detH.level(l) = -999.0;
        _invH.level(l).set_size(domain->numDim(), domain->numDim());
        _shapeSize.level(l).set_size(domain->numDim());
        _meshIndex.level(l).init(domain->numDim());
    }
    _U.set_size(domain->numDim(), domain->numDim());
    _V.set_size(domain->numDim(), domain->numDim());
    _S.set_size(domain->numDim());
    longAxisVertexY.init(domain->numDim());
    longAxisVertexX.init(domain->numDim());
    _quadraturePoints = new QuadraturePoints(this);
    _skeletonPoints = new SkeletonPoints(this);
    _tracers = new Tracers(this);
} // init

void Parcel::
updateDeformMatrix(const TimeLevelIndex<2> &timeIdx) {
    mat &H = _H.level(timeIdx);
    _skeletonPoints->updateLocalSpaceCoords(timeIdx);
    const field<BodyCoord> &y = SkeletonPoints::bodyCoords();
    const field<SpaceCoord> &x = _skeletonPoints->localSpaceCoords(timeIdx);
    const field<MeshIndex> &meshIndexs = _skeletonPoints->meshIndexs(timeIdx);
    if (domain->numDim() == 2) {
        // Calculate the elements of four matrices.
        double h11_1 = x[0](0)/y[0](0); double h21_1 = x[0](1)/y[0](0);
        double h12_2 = x[1](0)/y[1](1); double h22_2 = x[1](1)/y[1](1);
        double h11_3 = x[2](0)/y[2](0); double h21_3 = x[2](1)/y[2](0);
        double h12_4 = x[3](0)/y[3](1); double h22_4 = x[3](1)/y[3](1);
        // The final matrix is the combination of the four.
        H(0, 0) = (h11_1+h11_3)*0.5; H(0, 1) = (h12_2+h12_4)*0.5;
        H(1, 0) = (h21_1+h21_3)*0.5; H(1, 1) = (h22_2+h22_4)*0.5;
    } else if (domain->numDim() == 3) {
        // Calculate the elements of six matrices.
        double h11_1 = 0, h11_3 = 0, h12_2 = 0, h12_4 = 0, h13_5 = 0, h13_6 = 0;
        double h21_1 = 0, h21_3 = 0, h22_2 = 0, h22_4 = 0, h23_5 = 0, h23_6 = 0;
        double h31_1 = 0, h31_3 = 0, h32_2 = 0, h32_4 = 0, h33_5 = 0, h33_6 = 0;
        double w1 = 0, w2 = 0, w3 = 0, w4 = 0, w5 = 0, w6 = 0;
        if (meshIndexs[0].isValid()) {
            w1 = 1; h11_1 = x[0](0)/y[0](0); h21_1 = x[0](1)/y[0](0); h31_1 = x[0](2)/y[0](0);
        }
        if (meshIndexs[1].isValid()) {
            w2 = 1; h12_2 = x[1](0)/y[1](1); h22_2 = x[1](1)/y[1](1); h32_2 = x[1](2)/y[1](1);
        }
        if (meshIndexs[2].isValid()) {
            w3 = 1; h11_3 = x[2](0)/y[2](0); h21_3 = x[2](1)/y[2](0); h31_3 = x[2](2)/y[2](0);
        }
        if (meshIndexs[3].isValid()) {
            w4 = 1; h12_4 = x[3](0)/y[3](1); h22_4 = x[3](1)/y[3](1); h32_4 = x[3](2)/y[3](1);
        }
        if (meshIndexs[4].isValid()) {
            w5 = 1; h13_5 = x[4](0)/y[4](2); h23_5 = x[4](1)/y[4](2); h33_5 = x[4](2)/y[4](2);
        }
        if (meshIndexs[5].isValid()) {
            w6 = 1; h13_6 = x[5](0)/y[5](2); h23_6 = x[5](1)/y[5](2); h33_6 = x[5](2)/y[5](2);
        }
        // The final matrix is the combination of the six.
        if (w1 == 0 && w3 == 0) {
            REPORT_ERROR("CHECK!");
        } else {
            H(0, 0) = (w1*h11_1+w3*h11_3)/(w1+w3);
            H(1, 0) = (w1*h21_1+w3*h21_3)/(w1+w3);
            H(2, 0) = (w1*h31_1+w3*h31_3)/(w1+w3);
        }
        if (w2 == 0 && w4 == 0) {
            REPORT_ERROR("CHECK!");
        } else {
            H(0, 1) = (w2*h12_2+w4*h12_4)/(w2+w4);
            H(1, 1) = (w2*h22_2+w4*h22_4)/(w2+w4);
            H(2, 1) = (w2*h32_2+w4*h32_4)/(w2+w4);
        }
        if (w5 == 0 && w6 == 0) {
            REPORT_ERROR("CHECK!");
        } else {
            H(0, 2) = (w5*h13_5+w6*h13_6)/(w5+w6);
            H(1, 2) = (w5*h23_5+w6*h23_6)/(w5+w6);
            H(2, 2) = (w5*h33_5+w6*h33_6)/(w5+w6);
        }
    }
    // Set the spatial coordinates of the skeleton points that are out of range.
    for (uword i = 0; i < SkeletonPoints::numPoint(); ++i) {
        if (!meshIndexs[i].isValid()) {
            calcSpaceCoord(timeIdx, y[i], _skeletonPoints->spaceCoord(timeIdx, i));
        }
    }
    if (!svd(_U, _S, _V, H)) {
        REPORT_ERROR("Failed to do SVD on a matrix!");
    }
    // Eliminate the discrepency between detH and volume.
    double detH = arma::prod(_S);
    if (_detH.level(timeIdx) != -999.0) {
        _S *= pow(_detH.level(timeIdx)/detH, 1.0/domain->numDim());
    } else {
       _detH.level(timeIdx) = detH;
    }
    H = _U*diagmat(_S)*_V;
    _invH.level(timeIdx) = inv(H);
    longAxisVertexY() = _invH.level(timeIdx)*H*_V.col(0);
    calcSpaceCoord(timeIdx, longAxisVertexY, longAxisVertexX);
    _filament = _S[0]/_S[1];
} // updateDeformMatrix

void Parcel::
updateDeformMatrix(const TimeLevelIndex<2> &timeIdx, const vec &S) {
    _S = S;
    _H.level(timeIdx) = _U*diagmat(_S)*_V;
    _detH.level(timeIdx) = arma::prod(_S);
    _invH.level(timeIdx) = inv(_H.level(timeIdx));
    longAxisVertexY() = _invH.level(timeIdx)*_H.level(timeIdx)*_V.col(0);
    calcSpaceCoord(timeIdx, longAxisVertexY, longAxisVertexX);
    _filament = _S[0]/_S[1];
} // updateDeformMatrix

void Parcel::
resetSkeletonPoints(const TimeLevelIndex<2> &timeIdx, const Mesh &mesh) {
    const field<BodyCoord> &y = SkeletonPoints::bodyCoords();
    field<SpaceCoord> &x = _skeletonPoints->spaceCoords(timeIdx);
    field<MeshIndex> &meshIndexs = _skeletonPoints->meshIndexs(timeIdx);
    for (uword i = 0; i < x.size(); ++i) {
        calcSpaceCoord(timeIdx, y[i], x[i]);
        meshIndexs[i].reset();
        meshIndexs[i].locate(mesh, x[i]);
    }
} // resetSkeletonPoints

void Parcel::
updateVolume(const TimeLevelIndex<2> &timeIdx, double volume) {
    _detH.level(timeIdx) = volume;
} // updateVolume

void Parcel::
updateShapeSize(const TimeLevelIndex<2> &timeIdx) {
    BodyCoord y(domain->numDim());
    SpaceCoord x(domain->numDim());
    for (uword m = 0; m < domain->numDim(); ++m) {
        y() = _invH.level(timeIdx)*_H.level(timeIdx)*_V.col(m);
        calcSpaceCoord(timeIdx, y, x);
        _shapeSize.level(timeIdx)[m] = domain->calcDistance(x, _x.level(timeIdx));
    }
} // updateShapeSize

void Parcel::
calcSpaceCoord(const TimeLevelIndex<2> &timeIdx,
               const BodyCoord &y, SpaceCoord &x) const {
#ifdef LASM_IN_SPHERE
    x() = H(timeIdx)*y();
    domain->projectBack(geomtk::SphereDomain::STEREOGRAPHIC,
                        this->x(timeIdx), x, x());
#elif defined LASM_IN_CARTESIAN
    x() = this->x(timeIdx)()+H(timeIdx)*y();
    domain->isValid(x);
#endif
} // calcSpaceCoord

void Parcel::
calcBodyCoord(const TimeLevelIndex<2> &timeIdx,
              const SpaceCoord &x, BodyCoord &y) const {
#ifdef LASM_IN_SPHERE
    
    domain->project(geomtk::SphereDomain::STEREOGRAPHIC,
                    this->x(timeIdx), x, y());
    y() = invH(timeIdx)*y();
#elif defined LASM_IN_CARTESIAN
    y() = invH(timeIdx)*domain->diffCoord(x, this->x(timeIdx));
#endif
} // calcBodyCoord

double Parcel::
shapeFunction(const TimeLevelIndex<2> &timeIdx,
              const BodyCoord &y) const {
    double f;
    ShapeFunction::evalFunc(y, f);
    return f;
} // shapeFunction

void Parcel::
connectCell(int cellIdx) {
    for (uword i = 0; i < _numConnectedCell; ++i) {
        if (_connectedCellIndexs[i] == cellIdx) {
            return;
        }
    }
    if (_numConnectedCell == _connectedCellIndexs.size()) {
        _connectedCellIndexs.push_back(cellIdx);
    } else {
        _connectedCellIndexs[_numConnectedCell] = cellIdx;
    }
    _numConnectedCell++;
} // connectCell

void Parcel::
resetConnectedCells() {
    _numConnectedCell = 0;
} // resetConnectedCells

void Parcel::
dump(const TimeLevelIndex<2> &timeIdx,
     const MeshAdaptor &meshAdaptor) const {
    std::ofstream file; file.open("parcel_dump.txt");
    // Parcel centroid.
    file << "centroid = (/" << x(timeIdx)(0)/RAD << ",";
    file << x(timeIdx)(1)/RAD << "/)" << endl;
    // Parcel skeleton points.
    file << "skel_points = new((/4,2/), double)" << endl;
    for (int m = 0; m < 2; ++m) {
        file << "skel_points(:," << m << ") = (/";
        for (int i = 0; i < _skeletonPoints->spaceCoords(timeIdx).size(); ++i) {
            file << _skeletonPoints->spaceCoords(timeIdx)[i](m)/RAD;
            if (i != 3) {
                file << ",";
            } else {
                file << "/)" << endl;
            }
        }
    }
    // Parcel shape.
    int n = 100;
    file << "shape = new((/" << n << ",2/), double)" << endl;
    double dtheta = PI2/(n-1);
    SpaceCoord x(2); BodyCoord y(2);
    vector<vec::fixed<2> > shape(n);
    for (int i = 0; i < shape.size(); ++i) {
        double theta = i*dtheta;
        y(0) = cos(theta);
        y(1) = sin(theta);
        calcSpaceCoord(timeIdx, y, x);
        shape[i] = x();
    }
    for (int m = 0; m < 2; ++m) {
        file << "shape(:," << m << ") = (/";
        for (int i = 0; i < shape.size(); ++i) {
            if (i != shape.size()-1) {
                file << shape[i][m]/RAD << ",";
            } else {
                file << shape[i][m]/RAD << "/)" << endl;
            }
        }
    }
    // Neighbor cells.
    if (_numConnectedCell != 0) {
        file << "ngb_cells = new((/" << _numConnectedCell << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_cells(:," << m << ") = (/";
            for (int i = 0; i < _numConnectedCell; ++i) {
                int j = _connectedCellIndexs[i];
                file << meshAdaptor.coord(j)(m)/RAD;
                if (i != _numConnectedCell-1) {
                    file << ",";
                } else {
                    file << "/)" << endl;
                }
            }
        }
    }
    // Neighbor parcels.
    int numNeighborParcel = 0;
    for (int i = 0; i < _numConnectedCell; ++i) {
        int j = _connectedCellIndexs[i];
        numNeighborParcel += meshAdaptor.numContainedParcel(j);
    }
    if (numNeighborParcel != 0) {
        file << "ngb_tracers = new((/" << numNeighborParcel << ",2/), double)" << endl;
        for (int m = 0; m < 2; ++m) {
            file << "ngb_tracers(:," << m << ") = (/";
            int k = 0;
            for (int i = 0; i < _numConnectedCell; ++i) {
                const vector<Parcel*> &parcels = meshAdaptor.containedParcels(_connectedCellIndexs[i]);
                for (int j = 0; j < meshAdaptor.numContainedParcel(_connectedCellIndexs[i]); ++j) {
                    file << parcels[j]->x(timeIdx)(m)/RAD;
                    if (k != numNeighborParcel-1) {
                        file << ",";
                    } else {
                        file << "/)" << endl;
                    }
                    k++;
                }
            }
        }
    }
    file.close();
}

} // lasm
