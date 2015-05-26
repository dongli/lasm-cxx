#include "Tracers.h"
#include "Parcel.h"

namespace lasm {

vector<string> Tracers::_names;
vector<string> Tracers::_comments;
vector<string> Tracers::_units;
vec Tracers::_totalMasses;

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
    _totalMasses.resize(_names.size());
} // add

void Tracers::
init() {
} // init

void Tracers::
add() {
    _masses.resize(_names.size());
    _densities.resize(_names.size());
    _tendencies.resize(_names.size());
} // add

void Tracers::
resetDensities() {
    _densities.fill(0);
} // resetDensities

void Tracers::
resetTendencies() {
    _tendencies.fill(0);
} // resetTendencies

} // lasm