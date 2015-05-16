#ifndef __LASM_Tracers__
#define __LASM_Tracers__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class Tracers {
    static vector<string> _names;
    static vector<string> _units;
    static vector<string> _comments;

    vec _masses;
    vec _densities;
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

    double&
    mass(int speciesIdx) {
        return _masses(speciesIdx);
    }

    double&
    density(int speciesIdx) {
        return _densities(speciesIdx);
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