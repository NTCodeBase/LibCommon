//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//        __  __        _        __  ___ ____   __  ___
//       / / / /____ _ (_)_____ /  |/  // __ \ /  |/  /
//      / /_/ // __ `// // ___// /|_/ // /_/ // /|_/ /
//     / __  // /_/ // // /   / /  / // ____// /  / /
//    /_/ /_/ \__,_//_//_/   /_/  /_//_/    /_/  /_/
//
//    This file is part of HairMPM - Material Point Method for Hair Simulation.
//    Created: 2019. All rights reserved.
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Array/Array.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/Timer/Timer.h>
#include <LibCommon/ParallelHelpers/AtomicOperations.h>
#include <LibCommon/ParallelHelpers/Scheduler.h>

namespace TestArray {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
using namespace NTCodeBase;
using Real_t = float;

//#define FAST_TEST

#ifdef FAST_TEST
#define PERFORMANCE_TEST_NUM 1
#else
#define PERFORMANCE_TEST_NUM 10
#endif

#define NUM_ELEMENTS         1000000u
#define ARRAY_SIZE           Vec3ui(100, 800, 300)

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
auto generateData = [] {
                        StdVT<Vec3ui> result;
                        result.reserve(NUM_ELEMENTS);
                        for(UInt i = 0; i < NUM_ELEMENTS; ++i) {
                            result.push_back(Vec3ui(rand() % ARRAY_SIZE[0],
                                                    rand() % ARRAY_SIZE[1],
                                                    rand() % ARRAY_SIZE[2]));
                        }
                        return result;
                    };

////////////////////////////////////////////////////////////////////////////////
template<bool USE_Z_ORDER>
double test_function(const StdVT<Vec3ui>& indices) {
    Timer                         timer;
    double                        totalTime = 0;
    Real_t                        data_array0;
    Array<3, Real_t, USE_Z_ORDER> data_array;
    data_array.resize(ARRAY_SIZE);

    for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
        timer.tick();
        data_array.assign(0);
        Scheduler::parallel_for(NUM_ELEMENTS,
                                [&](UInt i) {
                                    Vec3ui idx = indices[i];
                                    AtomicOps::add(data_array(idx[0], idx[1], idx[2]), Real_t(1));
                                });
        data_array0 = data_array(0, 0, 0);
        totalTime  += timer.tock();
    }
    printf("Use z_order: %s, array size = %u-%u-%u, array(0, 0, 0) = %f, time = %f ms (%f s)\n", (USE_Z_ORDER ? "True" : "False"),
           data_array.vsize()[0], data_array.vsize()[1], data_array.vsize()[2],
           data_array0, float(totalTime), float(totalTime / 1000.0));
    return totalTime;
}

double test_function_change_size(const StdVT<Vec3ui>& indices) {
    Timer                   timer;
    double                  totalTime = 0;
    Real_t                  data_array0;
    Array<3, Real_t, false> data_array;
    auto                    maxSize = glm::compMax(ARRAY_SIZE);
    data_array.resize(Vec3ui(NumberHelpers::nextPowerOfTwo(maxSize)));

    for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
        timer.tick();
        data_array.assign(0);
        Scheduler::parallel_for(NUM_ELEMENTS,
                                [&](UInt i) {
                                    Vec3ui idx = indices[i];
                                    AtomicOps::add(data_array(idx[0], idx[1], idx[2]), Real_t(1));
                                });
        data_array0 = data_array(0, 0, 0);
        totalTime  += timer.tock();
    }
    printf("Use z_order: %s, array size = %u-%u-%u, array(0, 0, 0) = %f, time = %f ms (%f s)\n", "False",
           data_array.vsize()[0], data_array.vsize()[1], data_array.vsize()[2],
           data_array0, float(totalTime), float(totalTime / 1000.0));
    return totalTime;
}

////////////////////////////////////////////////////////////////////////////////
void test_array_access() {
    const auto indices              = generateData();
    auto       normalTime           = test_function<false>(indices);
    auto       normalTimeChangeSize = test_function_change_size(indices);
    auto       zTime = test_function<true>(indices);
    printf("Normal time: %f\n",             float(normalTime ));
    printf("Normal time change size: %f\n", float(normalTimeChangeSize ));
    printf("z time: %f\n",                  float(zTime ));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}
