#include "lasm.h"

using namespace lasm;

template <class TestType>
void run(const ConfigManager &configManager) {
    TestType test;
    AdvectionManager advectionManager;
    TimeLevelIndex<2> oldIdx, newIdx;

    test.init(configManager, advectionManager);
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
    ConfigManager configManager;

    configManager.parse(argv[1]);

    std::string caseName = configManager.getValue<std::string>("test_case", "case_name");
    if (caseName == "deform") {
        run<DeformationTest>(configManager);
    } else if (caseName == "terminator_chemistry") {
        run<TerminatorChemistryTest>(configManager);
    }

    return 0;
}