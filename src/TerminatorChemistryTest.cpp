#ifdef LASM_IN_SPHERE

#include "TerminatorChemistryTest.h"

#include "AdvectionManager.h"
#include "Parcel.h"
#include "Tracers.h"

TerminatorChemistryTest::
TerminatorChemistryTest() {
    subcase = "case4";
    period = 12*TimeUnit::DAYS;
    _stepSize = 1800*TimeUnit::SECONDS;
    _endTime = _startTime+period;
    k1Center.init(2);
    k1Center.set(300*RAD, 20*RAD);
    REPORT_ONLINE;
}

TerminatorChemistryTest::
~TerminatorChemistryTest() {
    REPORT_OFFLINE;
}

void TerminatorChemistryTest::
init(const ConfigManager &configManager, AdvectionManager &advectionManager) {
    // Create domain and mesh objects.
    _domain = new Domain(2);
    _domain->radius() = 6.3172e6;;
    _mesh = new Mesh(*_domain);
    // Initialize mesh.
    int numLon = configManager.getValue("test_case", "num_lon", 360);
    int numLat = configManager.getValue("test_case", "num_lat", 181);
    _mesh->init(numLon, numLat);
    // Initialize time manager.
    _stepSize = configManager.getValue("test_case", "time_step_size_in_seconds", _stepSize);
    _timeManager.init(_startTime, _endTime, _stepSize);
    // Initialize IO manager.
    std::string outputPattern = "lasm.terminator_chemistry."+
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
    // Initialize two reactive tracers.
    numTracer = 3;
    advectionManager.addTracer("q0", "1", "background tracer");
    advectionManager.addTracer("q1", "1", "X");
    advectionManager.addTracer("q2", "1", "X2");
    for (int t = 0; t < 3; ++t) {
        io.file(outputIdx).addField("double", FULL_DIMENSION,
                                    {&advectionManager.density(t)});
#ifdef LASM_TENDENCY_ON_MESH
        io.file(outputIdx).addField("double", FULL_DIMENSION,
                                    {&advectionManager.tendency(t)});
#endif
    }
} // init

void TerminatorChemistryTest::
setInitialCondition(AdvectionManager &advectionManager) {
    // Set initial conditions for each tracer.
    double *q = new double[3*mesh().totalNumGrid(CENTER, 2)];
    int l = 0;
    // - background tracer
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        q[l++] = 1.0;
    }
    // - X
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double k1, k2;
        calcChemicalReactionRate(x, k1, k2);
        double r = k1/(4*k2);
        double det = sqrt(r*r+2*r*Xt0);
        q[l++] = det-r;
    }
    // - X2
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double k1, k2;
        calcChemicalReactionRate(x, k1, k2);
        double r = k1/(4*k2);
        double det = sqrt(r*r+2*r*Xt0);
        q[l++] = Xt0*0.5-(det-r)*0.5;
    }
    // Propagate initial conditions to advection manager.
    TimeLevelIndex<2> timeIdx;
    advectionManager.input(timeIdx, q);
    delete [] q;
    setVelocityField(timeIdx);
} // setInitialCondition

void TerminatorChemistryTest::
advancePhysics(const TimeLevelIndex<2> &timeIdx,
               AdvectionManager &advectionManager) {
    // Calculate tracer density tendencies.
#if defined LASM_TENDENCY_ON_MESH
    for (int i = 0; i < mesh().totalNumGrid(CENTER, 2); ++i) {
        const SpaceCoord &x = mesh().gridCoord(CENTER, i);
        double X = advectionManager.density(timeIdx, 1, i)/advectionManager.density(timeIdx, 0, i);
        double X2 = advectionManager.density(timeIdx, 2, i)/advectionManager.density(timeIdx, 0, i);
        calcTendency(x, X, X2, advectionManager.tendency(1, i),
                     advectionManager.tendency(2, i));
        advectionManager.tendency(1, i) *= advectionManager.density(timeIdx, 0, i);
        advectionManager.tendency(2, i) *= advectionManager.density(timeIdx, 0, i);
    }
#elif defined LASM_TENDENCY_ON_PARCEL
    auto _am = static_cast<lasm::AdvectionManager*>(&advectionManager);
    for (auto parcel : _am->parcels()) {
        const SpaceCoord &x = parcel->x(timeIdx);
        double X = parcel->tracers().density(1)/parcel->tracers().density(0);
        double X2 = parcel->tracers().density(2)/parcel->tracers().density(0);
        calcTendency(x, X, X2, parcel->tracers().tendency(1),
                     parcel->tracers().tendency(2));
        parcel->tracers().tendency(1) *= parcel->tracers().density(0)*parcel->volume(timeIdx);
        parcel->tracers().tendency(2) *= parcel->tracers().density(0)*parcel->volume(timeIdx);
    }
#endif
} // advancePhysics

void TerminatorChemistryTest::
output(const TimeLevelIndex<2> &timeIdx,
       const AdvectionManager &advectionManager) {
    io.create(outputIdx);
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.output<double, 2>(outputIdx, timeIdx, {&advectionManager.density(t)});
#if defined LASM_TENDENCY_ON_MESH
        io.output<double>(outputIdx, {&advectionManager.tendency(t)});
#endif
    }
    io.output<double, 2>(outputIdx, timeIdx, {&velocityField(0),
                                              &velocityField(1),
                                              &velocityField.divergence()});
    if (io.isFileActive(outputIdx)) {
        advectionManager.output(timeIdx, io.file(outputIdx).fileId);
    }
    io.close(outputIdx);
} // output

void TerminatorChemistryTest::
calcChemicalReactionRate(const SpaceCoord &x, double &k1, double &k2) {
    k1 = fmax(0, x.sinLat()*k1Center.sinLat()+x.cosLat()*k1Center.cosLat()*
                 cos(x(0)-k1Center(0)));
    k2 = 1;
} // calcChemicalReactionRate

void TerminatorChemistryTest::
calcTendency(const SpaceCoord &x, double X, double X2, double &dX, double &dX2) {
    double k1, k2;

    calcChemicalReactionRate(x, k1, k2);

    double r = k1/(4*k2);
    double Xt = X+2*X2;

    double det = sqrt(r*r+2*r*Xt);
    double expdt = exp(-4*k2*det*stepSize());

    double el;
    if (fabs(det*k2*stepSize()) > 1e-16) {
        el = (1-expdt)/det/stepSize();
    } else {
        el = 4*k2;
    }

    dX = -el*(X-det+r)*(X+det+r)/(1+expdt+stepSize()*el*(X+r));
    dX2 = -dX*0.5;
} // calcTendency

#endif // LASM_IN_SPHERE
