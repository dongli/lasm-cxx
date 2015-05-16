#include "Cartesian3DTest.h"

Cartesian3DTest::
Cartesian3DTest() {
    REPORT_ONLINE;
}

Cartesian3DTest::
~Cartesian3DTest() {
    REPORT_OFFLINE;
}

void Cartesian3DTest::
init(const ConfigManager &configManager, AdvectionManager &advectionManager) {
    // Create domain and mesh objects.
    _domain = new Domain(3);
    _mesh = new Mesh(*_domain);
    // Read in configuration.
    std::string dataRoot = configManager.getValue<std::string>("test_case", "data_root");
    std::string dataPattern = configManager.getValue<std::string>("test_case", "data_pattern");
    std::string outputPattern = configManager.getValue<std::string>("test_case", "output_pattern");
    TimeStepUnit freqUnit = geomtk::timeStepUnitFromString(configManager.getValue<std::string>("test_case", "output_frequency_unit"));
    int freq = configManager.getValue<int>("test_case", "output_frequency");
    // Initialize time manager.
    _stepSize = io.file(dataIdx).getAttribute<float>("time_step_size_in_seconds");
    // Initialize IO manager.
    io.init(_timeManager);
    dataIdx = io.addInputFile(*_mesh, dataRoot+"/"+dataPattern);
    outputIdx = io.addOutputFile(*_mesh, outputPattern, freqUnit, freq);
    // Read in domain from the first data.
    _domain->init(io.file(dataIdx).currentFilePath());
    // Read in mesh from the first data.
    _mesh->init(io.file(dataIdx).currentFilePath());
    // Initialize velocity field.
    velocityField.create(*_mesh, true, true);
    io.file(dataIdx).addField("double", FULL_DIMENSION,
                              {&velocityField(0),
                               &velocityField(1),
                               &velocityField(2)});
    io.file(outputIdx).addField("double", FULL_DIMENSION,
                                {&velocityField(0),
                                 &velocityField(1),
                                 &velocityField(2),
                                 &velocityField.divergence()});
    // Initialize advection manager.
    advectionManager.init(configManager, *_mesh);
    // Initialize density fields for input and output.
    int numTracer = io.file(dataIdx).getAttribute<int>("num_tracer");
    densities.resize(numTracer);
    for (int t = 0; t < numTracer; ++t) {
        std::string name = "q"+std::to_string(t);
        std::string longName = io.file(dataIdx).getAttribute<std::string>(name, "long_name");
        std::string units = io.file(dataIdx).getAttribute<std::string>(name, "units");
        densities[t].create(name, units, longName, *_mesh, CENTER, _domain->numDim());
        io.file(dataIdx).addField("double", FULL_DIMENSION, {&densities[t]});
        io.file(outputIdx).addField("double", FULL_DIMENSION, {&densities[t]});
        advectionManager.addTracer(name, units, longName);
    }
} // init

void Cartesian3DTest::
setInitialCondition(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    io.open(dataIdx);
    io.input<double>(dataIdx, timeIdx, {&velocityField(0),
                                        &velocityField(1),
                                        &velocityField(2)});
    double *q = new double[_mesh->totalNumGrid(CENTER)*densities.size()];
    int j = 0;
    for (arma::uword t = 0; t < densities.size(); ++t) {
        io.input<double>(dataIdx, {&densities[t]});
        for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
            q[j++] = densities[t](i);
        }
    }
    advectionManager.input(timeIdx, q);
    io.close(dataIdx);
} // setInitialCondition

void Cartesian3DTest::
advanceDynamics(const TimeLevelIndex<2> &timeIdx,
                AdvectionManager &advectionManager) {
    io.open(dataIdx);
    io.input<double>(dataIdx, timeIdx, {&velocityField(0),
                                        &velocityField(1),
                                        &velocityField(2)});
    io.close(dataIdx);
    Interface::advanceDynamics(timeIdx, advectionManager);
} // advanceDynamics

void Cartesian3DTest::
output(const TimeLevelIndex<2> &timeIdx, AdvectionManager &advectionManager) {
    io.create(outputIdx);
    for (arma::uword t = 0; t < densities.size(); ++t) {
        for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
            densities[t](i) = advectionManager.density(timeIdx, t, i);
        }
        io.output<double>(outputIdx, {&densities[t]});
    }
    io.output<double, 2>(outputIdx, timeIdx, {&velocityField(0),
        &velocityField(1),  &velocityField(2), &velocityField.divergence()});
    io.close(outputIdx);
} // output