#include "DeformationTest.h"

#ifdef LASM_IN_SPHERE

DeformationTest::
DeformationTest() {
    period = 5;
    _stepSize = period/600;
    _endTime = _startTime+period;
    REPORT_ONLINE;
}

DeformationTest::
~DeformationTest() {
    REPORT_OFFLINE;
}

void DeformationTest::
init(const ConfigManager &configManager, AdvectionManager &advectionManager) {
    subcase = configManager.getValue("test_case", "subcase",
                                     std::string("case4"));
    // Create domain and mesh objects.
    _domain = new Domain(2);
    _domain->radius() = 1;
    _mesh = new Mesh(*_domain);
    // Initialize mesh.
    int numLon = configManager.getValue("test_case", "num_lon", 240);
    int numLat = configManager.getValue("test_case", "num_lat", 121);
    _mesh->init(numLon, numLat);
    // Initialize time manager.
    _stepSize = configManager.getValue("test_case", "time_step_size_in_seconds", _stepSize);
    _timeManager.init(_startTime, _endTime, _stepSize);
    // Initialize IO manager.
    std::string outputPattern = "lasm.deform."+subcase+"."+
        boost::lexical_cast<std::string>(numLon)+"x"+
        boost::lexical_cast<std::string>(numLat)+".%5s.nc";
    outputPattern = configManager.getValue("test_case", "output_pattern", outputPattern);
    TimeStepUnit freqUnit = geomtk::timeStepUnitFromString(configManager.getValue<std::string>("test_case", "output_frequency_unit"));
    int freq = configManager.getValue<int>("test_case", "output_frequency");
    io.init(_timeManager);
    outputIdx = io.addOutputFile(*_mesh, outputPattern, freqUnit, freq);
    // Initialize velocity field.
    velocityField.create(*_mesh, true, UPDATE_HALF_LEVEL);
    io.file(outputIdx).addField("double", FULL_DIMENSION,
                                {&velocityField(0),
                                 &velocityField(1),
                                 &velocityField.divergence()});
    // Initialize advection manager.
    advectionManager.init(configManager, *_mesh);
    // Initialize density fields for output.
    numTracer = 8;
    advectionManager.addTracer("q0", "1", "backgroud tracer");
    advectionManager.addTracer("q1", "1", "cosine hills tracer");
    advectionManager.addTracer("q2", "1", "q1 correlated tracer");
    advectionManager.addTracer("q3", "1", "slotted cylinders tracer");
    advectionManager.addTracer("q4", "1", "Gaussian hills tracer");
    advectionManager.addTracer("q5", "1", "slotted cylinders tracer");
    advectionManager.addTracer("q6", "1", "displaced slotted cylinders tracer");
    advectionManager.addTracer("q7", "1", "slotted cylinders residual tracer");
    for (int t = 0; t < numTracer; ++t) {
        io.file(outputIdx).addField("double", FULL_DIMENSION,
                                    {&advectionManager.density(t)});
    }
} // init

void DeformationTest::
setInitialCondition(AdvectionManager &advectionManager) {
    // Set initial conditions for each tracer.
    SpaceCoord c0(2), c1(2);
    c0.set(M_PI*5.0/6.0, 0.0); c0.transformToCart(domain());
    c1.set(M_PI*7.0/6.0, 0.0); c1.transformToCart(domain());
    double hmax, r, g, a, b, c;
    double *q = new double[numTracer*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    // - background tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // - Cosine hills tracer
    hmax = 1, r = domain().radius()*0.5, g = 0.1, c = 0.9;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
        if (r0 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r0/r));
        } else if (r1 < r) {
            q[l++] = g+c*hmax*0.5*(1+cos(M_PI*r1/r));
        } else {
            q[l++] = g;
        }
    }
    // - Tracer correlated to cosine hills tracer
    a = -0.8, b = 0.9;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l] = a*pow(q[l-mesh().totalNumGrid(CENTER, 2)], 2)+b; l++;
    }
    // - Slotted cylinders tracer
    b = 0.1, c = 1.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
        if ((r0 <= r && fabs(x(0)-c0(0)) >= r/6.0) ||
            (r1 <= r && fabs(x(0)-c1(0)) >= r/6.0))
            q[l++] = c;
        else if (r0 <= r && fabs(x(0)-c0(0)) < r/6.0 &&
                 x(1)-c0(1) < -5.0/12.0*r)
            q[l++] = c;
        else if (r1 <= r && fabs(x(0)-c1(0)) < r/6.0 &&
                 x(1)-c1(1) > 5.0/12.0*r)
            q[l++] = c;
        else
            q[l++] = b;
    }
    // - Gaussian hills tracer
    hmax = 0.95, b = 5.0;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        arma::vec d0 = x.cartCoord()-c0.cartCoord();
        arma::vec d1 = x.cartCoord()-c1.cartCoord();
        q[l++] = hmax*(exp(-b*dot(d0, d0))+exp(-b*dot(d1, d1)));
    }
    // - Another slotted cylinders tracer
    c0.set(M_PI*3.0/4.0, 0.0); c0.transformToCart(domain());
    c1.set(M_PI*5.0/4.0, 0.0); c1.transformToCart(domain());
    b = 0.1, c = 1.0/3.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
        if ((r0 <= r && fabs(x(0)-c0(0)) >= r/6.0) ||
            (r1 <= r && fabs(x(0)-c1(0)) >= r/6.0))
            q[l++] = c;
        else if (r0 <= r && fabs(x(0)-c0(0)) < r/6.0 &&
                 x(1)-c0(1) < -5.0/12.0*r)
            q[l++] = c;
        else if (r1 <= r && fabs(x(0)-c1(0)) < r/6.0 &&
                 x(1)-c1(1) > 5.0/12.0*r)
            q[l++] = c;
        else
            q[l++] = b;
    }
    // - Displaced slotted cylinders tracer
    c0.set(M_PI*3.0/4.0, M_PI/18.0);
    c1.set(M_PI*5.0/4.0, -M_PI/18.0);
    b = 0.1, c = 2.0/3.0, r = 0.5;
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double r0 = domain().calcDistance(x, c0);
        double r1 = domain().calcDistance(x, c1);
        if ((r0 <= r && fabs(x(0)-c0(0)) >= r/6.0) ||
            (r1 <= r && fabs(x(0)-c1(0)) >= r/6.0))
            q[l++] = c;
        else if (r0 <= r && fabs(x(0)-c0(0)) < r/6.0 &&
                 x(1)-c0(1) < -5.0/12.0*r)
            q[l++] = c;
        else if (r1 <= r && fabs(x(0)-c1(0)) < r/6.0 &&
                 x(1)-c1(1) > 5.0/12.0*r)
            q[l++] = c;
        else
            q[l++] = b;
    }
    // - Slotted cylinders residual tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l] = 1.0-q[l-2*mesh().totalNumGrid(CENTER, 2)]-q[l-mesh().totalNumGrid(CENTER, 2)];
        l++;
    }
    // Propagate initial conditions to advection manager.
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
    delete [] q;
    setVelocityField(timeIdx);
} // setInitialCondition

void DeformationTest::
output(const TimeLevelIndex<2> &timeIdx,
       const AdvectionManager &advectionManager) {
    io.create(outputIdx);
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.output<double, 2>(outputIdx, timeIdx,
                             {&advectionManager.density(t)});
    }
    io.output<double, 2>(outputIdx, timeIdx, {&velocityField(0),
                                              &velocityField(1),
                                              &velocityField.divergence()});
    if (io.isFileActive(outputIdx)) {
        advectionManager.output(timeIdx, io.file(outputIdx).fileId);
    }
    io.close(outputIdx);
} // output

void DeformationTest::
setVelocityField(const TimeLevelIndex<2> &timeIdx) {
    double cosT = cos(PI*_timeManager.seconds()/period);
    double k, R = _domain->radius();
    // advance velocity
    if (subcase == "case1") {
        k = 2.4;
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(0)(timeIdx, i) = k*pow(sin(lon*0.5), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(1)(timeIdx, i) = k*0.5*sin(lon)*cos(lat)*cosT;
        }
    } else if (subcase == "case2") {
        k = 2.0;
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT;
        }
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(1)(timeIdx, i) = k*sin(lon*2.0)*cos(lat)*cosT;
        }
    } else if (subcase == "case3") {
        k = 5.0*R/period;
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(0).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(0)(timeIdx, i) = -k*pow(sin(lon), 2.0)*sin(lat*2.0)*pow(cos(lat), 2.0)*cosT;
        }
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(1).staggerLocation(), i);
            double lon = x(0), lat = x(1);
            velocityField(1)(timeIdx, i) = k*0.5*sin(lon)*pow(cos(lat), 3.0)*cosT;
        }
    } else if (subcase == "case4") {
        k = 10.0*R/period;
        double c1 = PI2*_timeManager.seconds()/period;
        double c2 = PI2*R/period;
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(0).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(0).staggerLocation(), i);
            double lon = x(0)-c1, lat = x(1);
            velocityField(0)(timeIdx, i) = k*pow(sin(lon), 2.0)*sin(lat*2.0)*cosT+c2*cos(lat);
        }
        for (int i = 0; i < _mesh->totalNumGrid(velocityField(1).staggerLocation(), 2); ++i) {
            const SpaceCoord &x = _mesh->gridCoord(velocityField(1).staggerLocation(), i);
            double lon = x(0)-c1, lat = x(1);
            velocityField(1)(timeIdx, i) = k*sin(lon*2.0)*cos(lat)*cosT;
        }
    }
    if (timeIdx.isCurrentIndex()) {
        velocityField.applyBndCond(timeIdx);
    } else {
        velocityField.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
} // setVelocityField

#endif // LASM_IN_SPHERE
