#include "gtest/gtest.h"
#include "Parcel_test.h"
#include "Tracers_test.h"
#include "ParcelManager_test.h"
#include "AdvectionManager_test.h"
#include "MeshAdaptor_test.h"

int main(int argc, char *argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}