#include "BarotropicTest.h"

BarotropicTest::
BarotropicTest() {
    REPORT_ONLINE;
}

BarotropicTest::
~BarotropicTest() {
    REPORT_OFFLINE;
}

void BarotropicTest::
init(AdvectionManager &advectionManager) {
    DataTest::init(advectionManager);
    numTracer = 2;
    advectionManager.addTracer("q0", "m-2", "background tracer");
    advectionManager.addTracer("q1", "m-2", "step test tracer");
    io.file(outputIdx).addField("double", FULL_DIMENSION,
            {&advectionManager.density(0), &advectionManager.density(1)});
} // init

void BarotropicTest::
setInitialCondition(AdvectionManager &advectionManager) {
    double *q = new double[numTracer*_mesh->totalNumGrid(CENTER)];
    int j = 0;
    for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
        q[j++] = 1;
    }
    for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
        const SpaceCoord &x = _mesh->gridCoord(CENTER, i);
        if (x(0) > 160*RAD && x(0) < 200*RAD &&
            x(1) > 10*RAD  && x(1) < 40*RAD) {
            q[j++] = 1;
        } else {
            q[j++] = 0.1;
        }
    }
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
    DataTest::setVelocityField(timeIdx);
} // setInitialCondition
