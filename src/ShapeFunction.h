#ifndef __LASM_ShapeFunction__
#define __LASM_ShapeFunction__

#include "lasm_commons.h"

namespace lasm {

class ShapeFunction {
private:
    static double J;
    // one-dimensional quadrature nodes and weights
    static vec _nodes;
    static vec _weights;
    static const Domain *domain;
    static double _maxValue;
public:
    static void
    init(const Domain &domain);

    static vec&
    nodes() {
        return _nodes;
    }

    static vec&
    weights() {
        return _weights;
    }

    static double
    maxValue() {
        return _maxValue;
    }

    static void
    evalFunc(const BodyCoord& y, double &f);

    static void
    evalDerv(const BodyCoord& y, vec &d);
}; // ShapeFunction

} // lasm

#endif // __LASM_ShapeFunction__