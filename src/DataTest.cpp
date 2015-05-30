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
init(const ConfigManager &configManager, AdvectionManager &advectionManager) {
    // Read in configuration.
    auto caseName = configManager.getValue<string>("test_case", "case_name");
    auto dataRoot = configManager.getValue<string>("test_case", "data_root");
    auto dataPattern = configManager.getValue<string>("test_case", "data_pattern");
    auto outputPattern = configManager.getValue<string>("test_case", "output_pattern");
    TimeStepUnit freqUnit = geomtk::timeStepUnitFromString(configManager.getValue<string>("test_case", "output_frequency_unit"));
    int freq = configManager.getValue<int>("test_case", "output_frequency");
    auto domainType = configManager.getValue<string>("test_case", "domain_type");
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
    io.init(_timeManager);
    dataIdx = io.addInputFile(*_mesh, dataRoot+"/"+dataPattern);
    outputIdx = io.addOutputFile(*_mesh, outputPattern, freqUnit, freq);
    // Initialize time manager.
    vector<string> filePaths = geomtk::SystemTools::getFilePaths(dataRoot, geomtk::StampString::wildcard(dataPattern));
    _stepSize = io.file(dataIdx).getAttribute<float>("time_step_size_in_seconds");
    _startTime = io.getTime(filePaths.front());
    _endTime = io.getTime(filePaths.back());
    _timeManager.init(_startTime, _endTime, _stepSize);
    // Read in domain from the first data.
    _domain->init(io.file(dataIdx).currentFilePath());
    // Read in mesh from the first data.
    _mesh->init(io.file(dataIdx).currentFilePath());
    // Initialize velocity field.
    velocityField.create(*_mesh, true, true);
    for (int m = 0; m < _domain->numDim(); ++m) {
        io.file(dataIdx).addField("double", FULL_DIMENSION, {&velocityField(m)});
        io.file(outputIdx).addField("double", FULL_DIMENSION, {&velocityField(m)});
    }
    io.file(outputIdx).addField("double", FULL_DIMENSION, {&velocityField.divergence()});
    // Initialize advection manager.
    advectionManager.init(configManager, *_mesh);
    // Initialize density fields for input and output.
    numTracer = io.file(dataIdx).getAttribute<int>("num_tracer");
    for (int t = 0; t < numTracer; ++t) {
        string name = "q"+to_string(t);
        string longName = io.file(dataIdx).getAttribute<string>(name, "long_name");
        string units = io.file(dataIdx).getAttribute<string>(name, "units");
        advectionManager.addTracer(name, units, longName);
        io.file(dataIdx).addField("double", FULL_DIMENSION,
                                  {&advectionManager.density(t)});
        io.file(outputIdx).addField("double", FULL_DIMENSION,
                                    {&advectionManager.density(t),
                                     &advectionManager.tendency(t)});
    }
} // init

void DataTest::
setInitialCondition(AdvectionManager &advectionManager) {
    TimeLevelIndex<2> timeIdx;
    io.open(dataIdx);
    io.input<double>(dataIdx, timeIdx, {&velocityField(0),
                                        &velocityField(1),
                                        &velocityField(2)});
    double *q = new double[numTracer*_mesh->totalNumGrid(CENTER)];
    int j = 0;
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.input<double>(dataIdx, timeIdx, {&advectionManager.density(t)});
        for (arma::uword i = 0; i < _mesh->totalNumGrid(CENTER); ++i) {
            q[j++] = advectionManager.density(t)(timeIdx, i);
        }
    }
    advectionManager.input(timeIdx, q);
    io.close(dataIdx);
} // setInitialCondition

void DataTest::
setVelocityField(const TimeLevelIndex<2> &timeIdx) {
    io.open(dataIdx);
    io.input<double>(dataIdx, timeIdx, {&velocityField(0),
                                        &velocityField(1),
                                        &velocityField(2)});
    io.close(dataIdx);
} // setVelocityField

void DataTest::
output(const TimeLevelIndex<2> &timeIdx,
       const AdvectionManager &advectionManager) {
    io.create(outputIdx);
    for (arma::uword t = 0; t < numTracer; ++t) {
        io.output<double, 2>(outputIdx, timeIdx,
                             {&advectionManager.density(t)});
    }
    io.output<double, 2>(outputIdx, timeIdx, {&velocityField(0),
        &velocityField(1),  &velocityField(2), &velocityField.divergence()});
    io.close(outputIdx);
} // output