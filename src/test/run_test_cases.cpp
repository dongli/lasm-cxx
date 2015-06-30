#include "lasm.h"

using namespace lasm;

template <class TestType>
void run() {
    TestType test;
    AdvectionManager advectionManager;
    TimeLevelIndex<2> oldIdx, newIdx;

    test.init(advectionManager);
    test.setInitialCondition(advectionManager);
    test.output(oldIdx, advectionManager);

    while (!test.timeManager().isFinished()) {
        newIdx = oldIdx+1;
        test.advance(newIdx, advectionManager);
        test.output(newIdx, advectionManager);
        oldIdx.shift();
    }
}

int main(int argc, char const *argv[])
{
    if (argc != 2) {
        REPORT_ERROR("Configure file path is needed!");
    }

    ConfigManager::parse(argv[1]);

    auto caseName = ConfigManager::getValue<std::string>("common", "case_name");
    if (caseName == "deform") {
#ifdef LASM_IN_SPHERE
        run<DeformationTest>();
#endif
    } else if (caseName == "terminator_chemistry") {
#ifdef LASM_IN_SPHERE
        run<TerminatorChemistryTest>();
#endif
    } else if (caseName == "barotropic") {
#ifdef LASM_IN_SPHERE
        run<BarotropicTest>();
#endif
    } else if (caseName == "data") {
        run<DataTest>();
    }

    return 0;
}