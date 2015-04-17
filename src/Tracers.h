#ifndef __LASM_Tracers__
#define __LASM_Tracers__

#include "lasm_commons.h"

namespace lasm {

class Parcel;

class Tracers {
    static vector<string> _names;
    static vector<string> _units;
    static vector<string> _comments;

    TimeLevels<vec, 2> _masses;
    const Parcel *hostParcel;
public:
    Tracers(const Parcel *hostParcel);
    virtual ~Tracers();

    const vec&
    masses(const TimeLevelIndex<2> &timeIdx) const {
        return _masses.level(timeIdx);
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
    mass(const TimeLevelIndex<2> &timeIdx, int speciesIdx) {
        return _masses.level(timeIdx)(speciesIdx);
    }

    double
    density(const TimeLevelIndex<2> &timeIdx, int speciesIdx) const;

    void
    init();

    void
    add();

    void
    reset(const TimeLevelIndex<2> &timeIdx);
}; // Tracers

} // lasm

#endif // __LASM_Tracers__