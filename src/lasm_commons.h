#ifndef __LASM_commons__
#define __LASM_commons__

#define LASM_IN_SPHERE
#define LASM_USE_RLL_MESH

#if defined LASM_IN_CARTESIAN
#include "geomtk/Cartesian.h"
#elif defined LASM_IN_SPHERE && defined LASM_USE_RLL_MESH
#include "geomtk/RLLSphere.h"
#endif

namespace lasm {

using arma::vec;
using arma::mat;
using arma::uword;
using arma::field;
using std::list;
using std::vector;
using std::string;
using std::cout;
using std::endl;
using std::setw;

using geomtk::LINEAR;

class Parcel;
typedef list<Parcel*> Parcels;

// Transition function.
static double
transitionFunction(double x0, double y0,
                   double x1, double y1,
                   double x) {
    if (x < x0) {
        return y0;
    } else if (x > x1) {
        return y1;
    } else {
        double dx = x1-x0;
        double dy = y1-y0;
        double t = (x-x0)/dx;
        return dy*(4-3*t)*pow(t, 3)+y0;
    }
} // transitionFunction

} // lasm

#endif // __LASM_commons__
