#include "Tracers.h"
#include "Parcel.h"

namespace lasm {

vector<string> Tracers::_names;
vector<string> Tracers::_comments;
vector<string> Tracers::_units;

Tracers::
Tracers(const Parcel *hostParcel) {
    this->hostParcel = hostParcel;
}

Tracers::
~Tracers() {
}

void Tracers::
add(const string &name, const string &unit, const string &comment) {
    if (find(_names.begin(), _names.end(), name) != _names.end()) {
        REPORT_ERROR("Tracer \"" << name << "\" has been added!");
    }
    _names.push_back(name);
    _units.push_back(unit);
    _comments.push_back(comment);
} // add

double Tracers::
density(const TimeLevelIndex<2> &timeIdx, int speciesIdx) const  {
    return _masses.level(timeIdx)(speciesIdx)/hostParcel->volume(timeIdx);
} // density

void Tracers::
init() {
} // init

void Tracers::
add() {
    for (int l = 0; l < _masses.numLevel(); ++l) {
        _masses.level(l).resize(_names.size());
    }
}

void Tracers::
reset(const TimeLevelIndex<2> &timeIdx) {
    for (int t = 0; t < numTracer(); ++t) {
        _masses.level(timeIdx)[t] = 0;
    }
} // reset

} // lasm