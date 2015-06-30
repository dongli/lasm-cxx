#include "DataTest.h"
#include <stdlib.h>

using namespace boost::xpressive;
using namespace std;

DataTest::
DataTest() {
    REPORT_ONLINE;
}

DataTest::
~DataTest() {
    REPORT_OFFLINE;
}

void DataTest::
init(AdvectionManager &advectionManager) {
    // Read in configuration.
    auto caseName = ConfigManager::getValue<string>("common", "case_name");
    auto dataRoot = ConfigManager::getValue<string>(caseName, "data_root");
    auto dataPattern = ConfigManager::getValue<string>(caseName, "data_pattern");
    auto outputPattern = ConfigManager::getValue<string>(caseName, "output_pattern");
    auto freq = ConfigManager::getValue<string>(caseName, "output_frequency");
    auto domainType = ConfigManager::getValue<string>(caseName, "domain_type");
    // Create domain and mesh objects.
    mark_tag tagDomainType(1), tagNumDomain(2);
    sregex reDomainType = (tagDomainType= +_w) >> ' ' >> (tagNumDomain= +_d) >> "d";
    smatch what;
    if (!regex_match(domainType, what, reDomainType)) {
        REPORT_ERROR("Invalid domain type attribute \"" << domainType << "\"!");
    }
    // Check if the domain type is matched.
    if ((what[tagDomainType] == "Cartesian" && !is_same<Domain, geomtk::CartesianDomain>::value) ||
        (what[tagDomainType] == "Sphere" && !is_same<Domain, geomtk::SphereDomain>::value)) {
        REPORT_ERROR("Domain type is not matched!");
    }
    _domain = new Domain(atoi(what[tagNumDomain].str().c_str()));
    _mesh = new Mesh(*_domain);
    // Initialize IO manager.
    IOManager::init(_timeManager);
    // Initialize time manager.
    vector<string> filePaths = geomtk::SystemTools::getFilePaths(dataRoot, geomtk::StampString::wildcard(dataPattern));
    _startTime = io.getTime(filePaths.front());
    _endTime = io.getTime(filePaths.back());
    _stepSize = ConfigManager::getValue<double>(caseName, "time_step_size_in_seconds");
    _timeManager.init(_startTime, _endTime, _stepSize);
    // Add input and output files.
    dataIdx = io.addInputFile(*_mesh, dataRoot+"/"+dataPattern);
    outputIdx = io.addOutputFile(*_mesh, outputPattern, geomtk::durationFromString(freq));
    // Read in domain from the first data.
    _domain->init(io.file(dataIdx).currentFilePath());
    // Read in mesh from the first data.
    _mesh->init(io.file(dataIdx).currentFilePath());
    // Initialize velocity field.
    auto useStagger = ConfigManager::getValue<bool>(caseName, "velocity_stagger");
    velocityField.create(*_mesh, useStagger, true);
    for (int m = 0; m < _domain->numDim(); ++m) {
        io.file(dataIdx).addField("double", FULL_DIMENSION, {&velocityField(m)});
        io.file(outputIdx).addField("double", FULL_DIMENSION, {&velocityField(m)});
    }
    io.file(outputIdx).addField("double", FULL_DIMENSION, {&velocityField.divergence()});
#ifdef LASM_USE_DIAG
    // Initialize diagnostics.
    Diagnostics::init(*_mesh, io);
    Diagnostics::addMetric<Field<int> >("nmp", "1", "mixed parcel number");
    Diagnostics::addMetric<Field<int> >("ncp1", "1", "contained parcel number");
    Diagnostics::addMetric<Field<int> >("ncp2", "1", "connected parcel number");
#endif
    // Initialize advection manager.
    advectionManager.init(*_mesh);
    // Initialize density fields for input and output.
    if (io.file(dataIdx).hasAttribute("num_tracer")) {
        numTracer = io.file(dataIdx).getAttribute<int>("num_tracer");
        for (int t = 0; t < numTracer; ++t) {
            string name = "q"+to_string(t);
            string longName = io.file(dataIdx).getAttribute<string>(name, "long_name");
            string units = io.file(dataIdx).getAttribute<string>(name, "units");
            advectionManager.addTracer(name, units, longName);
            io.file(dataIdx).addField("double", FULL_DIMENSION,
                                      {&advectionManager.density(t)});
            io.file(outputIdx).addField("double", FULL_DIMENSION,
                                        {&advectionManager.density(t)});
        }
    }
} // init

void DataTest::
setInitialCondition(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    io.open(dataIdx);
    double *q = new double[numTracer*_mesh->totalNumGrid(CENTER)];
    int j = 0;
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.input<double>(dataIdx, timeIdx, {&advectionManager.density(t)});
        for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
            q[j++] = advectionManager.density(t)(timeIdx, i);
        }
    }
    io.close(dataIdx);
    advectionManager.input(timeIdx, q);
    setVelocityField(timeIdx);
} // setInitialCondition

void DataTest::
setVelocityField(const TimeLevelIndex<2> &timeIdx) {
    io.open(dataIdx);
    for (arma::uword m = 0; m < _domain->numDim(); ++m) {
        io.input<double>(dataIdx, timeIdx, {&velocityField(m)});
    }
    if (timeIdx.isCurrentIndex()) {
        velocityField.applyBndCond(timeIdx);
    } else {
        velocityField.applyBndCond(timeIdx, UPDATE_HALF_LEVEL);
    }
    io.close(dataIdx);
} // setVelocityField

void DataTest::
output(const TimeLevelIndex<2> &timeIdx,
       const AdvectionManager &advectionManager) {
    io.create(outputIdx);
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.output<double>(outputIdx, timeIdx,
                             {&advectionManager.density(t)});
    }
    for (arma::uword m = 0; m < _domain->numDim(); ++m) {
        io.output<double>(outputIdx, timeIdx, {&velocityField(m)});
    }
    io.output<double>(outputIdx, timeIdx, {&velocityField.divergence()});
    if (io.isFileActive(outputIdx)) {
        advectionManager.output(timeIdx, io.file(outputIdx).fileId);
    }
    io.close(outputIdx);
#ifdef LASM_USE_DIAG
    Diagnostics::output();
#endif
} // output