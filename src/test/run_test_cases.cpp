#include "lasm.h"

using namespace lasm;

int main(int argc, char const *argv[])
{
    if (argc != 2) {
        REPORT_ERROR("Configure file path is needed!");
    }
    ConfigManager configManager;
    Cartesian3DTest test;
    AdvectionManager advectionManager;
    TimeLevelIndex<2> oldIdx;

    configManager.parse(argv[1]);

    test.init(configManager, advectionManager);
    test.setInitialCondition(advectionManager);
    test.output(oldIdx, advectionManager);
    return 0;
}