#include "QuadraturePoints.h"
#include "Parcel.h"
#include "ShapeFunction.h"

namespace lasm {

const Domain* QuadraturePoints::domain;
field<double> QuadraturePoints::w;
field<BodyCoord> QuadraturePoints::y;

QuadraturePoints::
QuadraturePoints(const Parcel *hostParcel) {
    this->hostParcel = hostParcel;
    if (domain->numDim() == 2) {
        x.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        for (uword j = 0; j < ShapeFunction::nodes().size(); ++j) {
            for (uword i = 0; i < ShapeFunction::nodes().size(); ++i) {
                x(i, j).setNumDim(domain->numDim());
            }
        }
    } else if (domain->numDim() == 3) {
        x.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        for (uword k = 0; k < ShapeFunction::nodes().size(); ++k) {
            for (uword j = 0; j < ShapeFunction::nodes().size(); ++j) {
                for (uword i = 0; i < ShapeFunction::nodes().size(); ++i) {
                    x(i, j, k).setNumDim(domain->numDim());
                }
            }
        }
    }
}

QuadraturePoints::
~QuadraturePoints() {

}

void QuadraturePoints::
init(const Domain &_domain) {
    domain = &_domain;
    // Set the body coordinates of quadrature points.
    if (domain->numDim() == 2) {
        w.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        y.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        for (uword j = 0; j < ShapeFunction::nodes().size(); ++j) {
            for (uword i = 0; i < ShapeFunction::nodes().size(); ++i) {
                y(i, j).setNumDim(domain->numDim());
                w(i, j) = ShapeFunction::weights()(i)*
                          ShapeFunction::weights()(j);
                y(i, j)(0) = ShapeFunction::nodes()(i);
                y(i, j)(1) = ShapeFunction::nodes()(j);
            }
        }
    } else if (domain->numDim() == 3) {
        w.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        y.set_size(ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size(),
                   ShapeFunction::nodes().size());
        for (uword k = 0; k < ShapeFunction::nodes().size(); ++k) {
            for (uword j = 0; j < ShapeFunction::nodes().size(); ++j) {
                for (uword i = 0; i < ShapeFunction::nodes().size(); ++i) {
                    y(i, j, k).setNumDim(domain->numDim());
                    w(i, j, k) = ShapeFunction::weights()(i)*
                                 ShapeFunction::weights()(j)*
                                 ShapeFunction::weights()(k);
                    y(i, j, k)(0) = ShapeFunction::nodes()(i);
                    y(i, j, k)(1) = ShapeFunction::nodes()(j);
                    y(i, j, k)(2) = ShapeFunction::nodes()(k);
                }
            }
        }
    }
} // init

void QuadraturePoints::
updateSpaceCoords(const TimeLevelIndex<2> &timeIdx) {
    for (uword i = 0; i < y.size(); ++i) {
        hostParcel->calcSpaceCoord(timeIdx, y(i), x(i));
    }
} // updateSpaceCoords

}