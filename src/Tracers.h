#ifndef __LASM_Tracers__
#define __LASM_Tracers__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class Tracers {
    static vector<string> _names;
    static vector<string> _units;
    static vector<string> _comments;
    static vec _totalMasses;

    vec _masses;
    vec _densities;
    vec _tendencies;
    const Parcel *hostParcel;
public:
    Tracers(const Parcel *hostParcel);
    virtual ~Tracers();

    const vec&
    masses() const {
        return _masses;
    }

    const vec&
    densities() const {
        return _densities;
    }

    static void
    add(const string &name, const string &unit, const string &comment);

    static const vector<string>&
    names() {
        return _names;
    }

    static const vector<string>&
    comments() {
        return _comments;
    }

    static const vector<string>&
    units() {
        return _units;
    }

    static int
    numTracer() {
        return _names.size();
    }

    static const vec&
    totalMasses() {
        return _totalMasses;
    }

    static double&
    totalMass(int tracerIdx) {
        return _totalMasses[tracerIdx];
    }

    static void
    resetTotalMasses() {
        _totalMasses.fill(0);
    }

    double&
    mass(int tracerIdx) {
        return _masses[tracerIdx];
    }

    double&
    density(int tracerIdx) {
        return _densities[tracerIdx];
    }

    double&
    tendency(int tracerIdx) {
        return _tendencies[tracerIdx];
    }

    void
    init();

    void
    add();

    void
    reset();
}; // Tracers

} // lasm

#endif // __LASM_Tracers__