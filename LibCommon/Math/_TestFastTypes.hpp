//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//        __  __        _        __  ___ ____   __  ___
//       / / / /____ _ (_)_____ /  |/  // __ \ /  |/  /
//      / /_/ // __ `// // ___// /|_/ // /_/ // /|_/ /
//     / __  // /_/ // // /   / /  / // ____// /  / /
//    /_/ /_/ \__,_//_//_/   /_/  /_//_/    /_/  /_/
//
//    This file is part of HairMPM - Material Point Method for Hair Simulation.
//    Created: 2018. All rights reserved.
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#include <catch2/catch.hpp>
#include <LibCommon/CommonSetup.h>
#include <LibCommon/Utils/NumberHelpers.h>
#include <LibCommon/Timer/ScopeTimer.h>

#include <LibCommon/Array/Array.h>
#include <LibCommon/Grid/Grid.h>
#include <LibCommon/Grid/FastGrid.h>
#include <LibCommon/ParallelHelpers/AtomicOperations.h>
#include <LibCommon/ParallelHelpers/Scheduler.h>

#include <LibCommon/Math/FastVec3.h>
#include <LibCommon/Math/FastMat3.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
using Real = float;

//#define FAST_TEST
#define PARALLEL_TEST

#define PRECISION            1e-7
#  define DATA_SIZE          100000

#ifdef FAST_TEST
#define PERFORMANCE_TEST_NUM 1
#else
#define PERFORMANCE_TEST_NUM 100
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//#define TEST_CORRECTNESS_INITIALIZATION
//#define TEST_PERFORMANCE_INITIALIZATION
//#define TEST_FAST_VEC3_OPS
//#define TEST_FAST_MAT3_OPS
#define TEST_PERFORMANCE_FAST_VEC3_FAST_MAT3

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
auto init_glmVec3 = [] (auto& v3_data) {
                        v3_data.resize(0);
                        v3_data.reserve(DATA_SIZE);
                        for(int i = 0; i < DATA_SIZE; ++i) {
                            v3_data.push_back(NumberHelpers::fRand<Real>::vrnd<Vec3<Real>>());
                        }
                    };

auto init_FastVec3 = [](auto& fv3_data) {
                         fv3_data.resize(0);
                         fv3_data.reserve(DATA_SIZE);
                         for(int i = 0; i < DATA_SIZE; ++i) {
                             fv3_data.push_back(NumberHelpers::fRand<Real>::vrnd<Vec3<Real>>());
                         }
                     };

auto init_glmMat3 = [] (auto& m3_data) {
                        m3_data.resize(0);
                        m3_data.reserve(DATA_SIZE);
                        for(int i = 0; i < DATA_SIZE; ++i) {
                            m3_data.push_back(NumberHelpers::fRand<Real>::vrnd<Mat3x3<Real>>());
                        }
                    };

auto init_FastMat3 = [](auto& fm3_data) {
                         fm3_data.resize(0);
                         fm3_data.reserve(DATA_SIZE);
                         for(int i = 0; i < DATA_SIZE; ++i) {
                             fm3_data.push_back(NumberHelpers::fRand<Real>::vrnd<Mat3x3<Real>>());
                         }
                     };

auto overwrite_glmMat3 = [] (auto& m3, auto& m3_data) {
                             for(int i = 0; i < DATA_SIZE; ++i) {
                                 m3_data[i] = m3;
                             }
                         };

auto overwrite_FastMat3ToMat3 = [](auto& fm3, auto& m3_data) {
                                    for(int i = 0; i < DATA_SIZE; ++i) {
                                        fm3.copyToMat3x3r(m3_data[i]);
                                    }
                                };

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
auto check_vec3 = [](const auto& v3, const auto& fv3) {
                      bool bResult =
                          ((std::abs(v3[0] - fv3[0]) < PRECISION) &&
                           (std::abs(v3[1] - fv3[1]) < PRECISION) &&
                           (std::abs(v3[2] - fv3[2]) < PRECISION) &&
                           //                           (std::abs(fv3[3]) < PRECISION) &&
                           //                           (std::abs(fv3.dummy) < PRECISION) &&
                           //
                           (std::abs(v3.x - fv3.x) < PRECISION) &&
                           (std::abs(v3.y - fv3.y) < PRECISION) &&
                           (std::abs(v3.z - fv3.z) < PRECISION)
                           //                           (std::abs(fv3.w) < PRECISION)
                          );

                      if(!bResult) {
                          auto err0 = std::abs(v3[0] - fv3[0]);
                          auto err1 = std::abs(v3[1] - fv3[1]);
                          auto err2 = std::abs(v3[2] - fv3[2]);

                          printf("Not correct: glm_v3 = %f, %f, %f  | FastVec3: %f, %f, %f, %f - %f, %f, %f, %f | err = %10.8e, %10.8e, %10.8e\n",
                                 v3.x, v3.y, v3.z,
                                 fv3.x, fv3.y, fv3.z, fv3.w,
                                 fv3[0], fv3[1], fv3[2], fv3[3],
                                 err0, err1, err2
                                 );
                      }

                      return bResult;
                  };
auto check_mat3 = [](const auto& m3, const auto& fm3) {
                      bool bResult_c0 = check_vec3(m3[0], fm3[0]);
                      bool bResult_c1 = check_vec3(m3[1], fm3[1]);
                      bool bResult_c2 = check_vec3(m3[2], fm3[2]);
                      bool bResult    = bResult_c0 && bResult_c1 && bResult_c2;

                      if(!bResult) {
                          printf("Not correct: glm_m3 = \n  %f, %f, %f\n  %f, %f, %f\n  %f, %f, %f  \n"
                                 "FastMat3: \n  %f, %f, %f, %f\n  %f, %f, %f, %f\n  %f, %f, %f, %f\n",
                                 m3[0].x, m3[0].y, m3[0].z,
                                 m3[1].x, m3[1].y, m3[1].z,
                                 m3[2].x, m3[2].y, m3[2].z,
                                 fm3[0].x, fm3[0].y, fm3[0].z, fm3[0].w,
                                 fm3[1].x, fm3[1].y, fm3[1].z, fm3[1].w,
                                 fm3[2].x, fm3[2].y, fm3[2].z, fm3[2].w
                                 );
                      }

                      return bResult;
                  };

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef TEST_CORRECTNESS_INITIALIZATION
TEST_CASE("Test_Correctness_Initialization", "Test_Correctness_Initialization")
{
    ScopeTimer          timer("Test_Correctness_Initialization");
    StdVT<Vec3<Real>>   v3_data;
    StdVT<Mat3x3<Real>> m3_data;
    ////////////////////////////////////////////////////////////////////////////////
    {
        init_glmVec3(v3_data);
        for(UInt i = 0; i < DATA_SIZE; ++i) {
            auto& v3 = v3_data[i];
            {
                FastVec3<Real> fv3 = v3;
                REQUIRE(check_vec3(v3, fv3));

                Vec3<Real> v3_copy = fv3;
                REQUIRE(check_vec3(v3_copy, fv3));
                REQUIRE(       glm::length2(v3 - v3_copy) < PRECISION);
            }
            {
                FastVec3<Real> fv3(v3.x, v3.y, v3.z);
                REQUIRE(check_vec3(v3, fv3));
            }
            {
                FastVec3<Real> fv3     = v3;
                Vec3i          vi      = Vec3i(v3);
                Vec3i          fvi     = fv3.toVec3i();
                bool           bResult = (vi.x == fvi.x && vi.y == fvi.y && vi.z == fvi.z);
                REQUIRE(bResult);
            }
        }
    }
    ////////////////////////////////////////////////////////////////////////////////
    {
        init_glmMat3(m3_data);
        for(UInt i = 0; i < DATA_SIZE; ++i) {
            auto& m3 = m3_data[i];
            {
                FastMat3<Real> fm3 = m3;
                REQUIRE(check_mat3(m3, fm3));

                Mat3x3<Real> m3_copy1 = fm3;
                REQUIRE(check_mat3(m3_copy1, fm3));
                REQUIRE(       glm::length2(m3[0] - m3_copy1[0]) + glm::length2(m3[1] - m3_copy1[1]) + glm::length2(m3[2] - m3_copy1[2]) < PRECISION);

                Mat3x3<Real> m3_copy2;
                fm3.copyToMat3x3r(m3_copy2);
                REQUIRE(check_mat3(m3_copy2, fm3));
                REQUIRE(       glm::length2(m3[0] - m3_copy2[0]) + glm::length2(m3[1] - m3_copy2[1]) + glm::length2(m3[2] - m3_copy2[2]) < PRECISION);
            }
            {
                FastMat3<Real> fm3(m3[0], m3[1], m3[2]);
                REQUIRE(check_mat3(m3, fm3));
            }
        }
    }
}
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef TEST_PERFORMANCE_INITIALIZATION
TEST_CASE("Test_Performance_Initialization", "Test_Performance_Initialization")
{
    ScopeTimer            timer("Test_Performance_Initialization");
    StdVT<Vec3<Real>>     v3_data;
    StdVT<FastVec3<Real>> fv3_data;

    StdVT<Mat3x3<Real>>   m3_data;
    StdVT<FastMat3<Real>> fm3_data;
    ////////////////////////////////////////////////////////////////////////////////
    {
        ScopeTimer timer("GLM vec3 initialization");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            init_glmVec3(v3_data);
        }
    }
    {
        ScopeTimer timer("FastVec3 initialization");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            init_FastVec3(fv3_data);
        }
    }
    {
        ScopeTimer timer("GLM mat3 initialization");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            init_glmMat3(m3_data);
        }
    }
    {
        ScopeTimer timer("FastMat3 initialization");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            init_FastMat3(fm3_data);
        }
    }
    {
        ScopeTimer  timer("Overwrite GLM mat3 to GLM mat3");
        const auto& m3 = m3_data.front();
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            overwrite_glmMat3(m3, m3_data);
        }
    }
    {
        ScopeTimer  timer("Overwrite FastMat3 to GLM mat3");
        const auto& fm3 = fm3_data.front();
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            overwrite_FastMat3ToMat3(fm3, m3_data);
        }
    }
}
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef TEST_FAST_VEC3_OPS
TEST_CASE("Test_FastVec3_Ops", "Test_FastVec3_Ops")
{
    ScopeTimer timer("Test_FastVec3_Ops");
    for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
#ifdef PARALLEL_TEST
        Scheduler::parallel_for(DATA_SIZE,
                                [&](int i) {
#else
        for(int i = 0; i < DATA_SIZE; ++i) {
#endif
                                    __NT_UNUSED(i);
                                    Vec3<Real> v0      = NumberHelpers::fRand11<Real>::vrnd<Vec3<Real>>();
                                    Vec3<Real> v1      = NumberHelpers::fRand11<Real>::vrnd<Vec3<Real>>();
                                    FastVec3<Real> fv0 = v0;
                                    FastVec3<Real> fv1 = v1;
                                    {
                                        auto v2_add = v0 + v1;
                                        auto v2_sub = v0 - v1;
                                        auto v2_mul = v0 * v1;
                                        auto v2_div = v0 / v1;

                                        auto fv2_add = fv0 + fv1;
                                        auto fv2_sub = fv0 - fv1;
                                        auto fv2_mul = fv0 * fv1;
                                        auto fv2_div = fv0 / fv1;

                                        REQUIRE(check_vec3(v2_add, fv2_add));
                                        REQUIRE(check_vec3(v2_sub, fv2_sub));
                                        REQUIRE(check_vec3(v2_mul, fv2_mul));
                                        REQUIRE(check_vec3(v2_div, fv2_div));
                                        ////////////////////////////////////////////////////////////////////////////////
                                        v2_add += v0;
                                        v2_sub -= v0;
                                        v2_mul *= v0;
                                        v2_div /= v0;

                                        fv2_add += fv0;
                                        fv2_sub -= fv0;
                                        fv2_mul *= fv0;
                                        fv2_div /= fv0;

                                        REQUIRE(check_vec3(v2_add, fv2_add));
                                        REQUIRE(check_vec3(v2_sub, fv2_sub));
                                        REQUIRE(check_vec3(v2_mul, fv2_mul));
                                        REQUIRE(check_vec3(v2_div, fv2_div));
                                        ////////////////////////////////////////////////////////////////////////////////
                                        auto r  = NumberHelpers::fRand11<Real>::rnd();
                                        v2_add += r;
                                        v2_sub -= r;
                                        v2_mul *= r;
                                        v2_div /= r;

                                        fv2_add += r;
                                        fv2_sub -= r;
                                        fv2_mul *= r;
                                        fv2_div /= r;

                                        REQUIRE(check_vec3(v2_add, fv2_add));
                                        REQUIRE(check_vec3(v2_sub, fv2_sub));
                                        REQUIRE(check_vec3(v2_mul, fv2_mul));
                                        REQUIRE(check_vec3(v2_div, fv2_div));
                                    }
                                    ////////////////////////////////////////////////////////////////////////////////
                                    {
                                        auto v0_norm  = glm::normalize(v0);
                                        auto fv0_norm = fv0.normalized();
                                        REQUIRE(check_vec3(v0_norm, fv0_norm));
                                    }
                                    {
                                        auto v_dot  = glm::dot(v0, v1);
                                        auto fv_dot = fv0.dot(fv1);
                                        REQUIRE(std::abs(v_dot - fv_dot) < PRECISION);
                                    }
                                    {
                                        auto r   = glm::length(v0);
                                        auto r2  = glm::length2(v0);
                                        auto fr  = fv0.norm();
                                        auto fr2 = fv0.norm2();
                                        REQUIRE(std::abs(r - fr) < PRECISION);
                                        REQUIRE(std::abs(r2 - fr2) < PRECISION);
                                    }
                                    {
                                        auto v2  = glm::cross(v0, v1);
                                        auto fv2 = fv0.cross(fv1);
                                        REQUIRE(check_vec3(v2, fv2));
                                    }
                                    {
                                        auto fv2 = fv0;
                                        REQUIRE(check_vec3(fv2, fv0));
                                    }
                                    {
                                        FastVec3<Real> fv2(0);
                                        FastVec3<Real> fv3 = -fv0;
                                        FastVec3<Real> fv4 = fv0 + fv3;
                                        REQUIRE(check_vec3(fv2, fv4));
                                    }
#ifdef PARALLEL_TEST
                                });
#else
                                }
#endif
    }
}

#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef TEST_FAST_MAT3_OPS
TEST_CASE("Test_FastMat3_Ops", "Test_FastMat3_Ops")
{
    ScopeTimer timer("Test_FastMat3_Ops");
    for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
#ifdef PARALLEL_TEST
        Scheduler::parallel_for(DATA_SIZE,
                                [&](int i) {
#else
        for(int i = 0; i < DATA_SIZE; ++i) {
#endif
                                    __NT_UNUSED(i);
                                    Vec3<Real> v0      = NumberHelpers::fRand11<Real>::vrnd<Vec3<Real>>();
                                    Mat3x3<Real> m0    = NumberHelpers::fRand11<Real>::vrnd<Mat3x3<Real>>();
                                    Mat3x3<Real> m1    = NumberHelpers::fRand11<Real>::vrnd<Mat3x3<Real>>();
                                    FastVec3<Real> fv0 = v0;
                                    FastMat3<Real> fm0 = m0;
                                    FastMat3<Real> fm1 = m1;
                                    {
                                        auto m2_add = m0 + m1;
                                        auto m2_sub = m0 - m1;
                                        auto m2_mul = m0 * m1;
                                        auto v2_mul = m0 * v0;
                                        auto m2_div = m2_mul;

                                        auto fm2_add = fm0 + fm1;
                                        auto fm2_sub = fm0 - fm1;
                                        auto fm2_mul = fm0 * fm1;
                                        auto fv2_mul = fm0 * fv0;
                                        auto fm2_div = fm2_mul;

                                        REQUIRE(check_mat3(m2_add, fm2_add));
                                        REQUIRE(check_mat3(m2_sub, fm2_sub));
                                        REQUIRE(check_mat3(m2_mul, fm2_mul));
                                        REQUIRE(check_vec3(v2_mul, fv2_mul));
                                        ////////////////////////////////////////////////////////////////////////////////
                                        m2_add += m0;
                                        m2_sub -= m0;
                                        m2_mul *= m0;

                                        fm2_add += fm0;
                                        fm2_sub -= fm0;
                                        fm2_mul *= fm0;

                                        REQUIRE(check_mat3(m2_add, fm2_add));
                                        REQUIRE(check_mat3(m2_sub, fm2_sub));
                                        REQUIRE(check_mat3(m2_mul, fm2_mul));
                                        ////////////////////////////////////////////////////////////////////////////////
                                        auto r  = NumberHelpers::fRand11<Real>::rnd();
                                        m2_add += r;
                                        m2_sub -= r;
                                        m2_mul *= r;
                                        m2_div /= r;

                                        fm2_add += r;
                                        fm2_sub -= r;
                                        fm2_mul *= r;
                                        fm2_div /= r;

                                        REQUIRE(check_mat3(m2_add, fm2_add));
                                        REQUIRE(check_mat3(m2_sub, fm2_sub));
                                        REQUIRE(check_mat3(m2_mul, fm2_mul));
                                        REQUIRE(check_mat3(m2_div, fm2_div));
                                    }
                                    ////////////////////////////////////////////////////////////////////////////////
                                    {
                                        auto m0_T  = glm::transpose(m0);
                                        auto fm0_T = fm0.transposed();
                                        REQUIRE(check_mat3(m0_T, fm0_T));
                                    }
                                    {
                                        auto fm2 = fm0;
                                        REQUIRE(check_mat3(fm2, fm0));
                                    }
                                    {
                                        FastMat3<Real> fm2(0);
                                        FastMat3<Real> fv3 = -fm0;
                                        FastMat3<Real> fv4 = fm0 + fv3;
                                        REQUIRE(check_mat3(fm2, fv4));
                                    }
#ifdef PARALLEL_TEST
                                });
#else
                                }
#endif
    }
}

#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef TEST_PERFORMANCE_FAST_VEC3_FAST_MAT3
TEST_CASE("Compare_FastVec3_Performance", "Compare_FastVec3_Performance")
{
    ScopeTimer timer("Compare_FastVec3_Performance");
    StdVT<Real> mass(DATA_SIZE);
    StdVT<Vec3<Real>> positions(DATA_SIZE);
    StdVT<Vec3<Real>> velocities(DATA_SIZE);
    StdVT<Vec3<int>> kernelCornerNode(DATA_SIZE);
    StdVT<Mat3x3<Real>> C(DATA_SIZE);;

    Grid<3, Real> grid;
    FastGrid3<Real> fastGrid;
    Array<3, Vec3<Real>> gridData;
    const auto cellSize = Real(2.0 / 128.0);
    grid.setGrid(Vec3<Real>(0 - 2.0 * cellSize), Vec3<Real>(1.0 + 2.0 * cellSize), cellSize, cellSize);
    fastGrid.setGrid(Vec3<Real>(0 - 2.0 * cellSize), Vec3<Real>(1.0 + 2.0 * cellSize), cellSize, cellSize);

    {
        ScopeTimer timer("Init performance test data");
        gridData.resize(grid.getNNodes());
        Scheduler::parallel_for(DATA_SIZE,
                                [&](int i) {
                                    mass[i]       = NumberHelpers::fRand<Real>::rnd();
                                    positions[i]  = NumberHelpers::fRand01<Real>::vrnd<Vec3<Real>>();
                                    velocities[i] = NumberHelpers::fRand<Real>::vrnd<Vec3<Real>>();
                                    C[i]          = NumberHelpers::fRand<Real>::mrnd<Mat3x3<Real>>();

                                    const auto pg       = grid.getGridCoordinate(positions[i]);
                                    kernelCornerNode[i] = Vec3i(pg) - Vec3i(1);
                                });
    }

    {
        ScopeTimer timer("Test velocity transfer using GLM");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto mp       = mass[i];
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto& ppos    = positions[i];
                                        const auto& pvel    = velocities[i];
                                        const auto& pC      = C[i];
                                        const auto& lcorner = kernelCornerNode[i];
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xixp = grid.getWorldCoordinate(grid_x, grid_y, grid_z) - ppos;
                                                    const auto vp   = (pvel + pC * xixp) * (w * mp);
                                                    AtomicOperations::add(gridData(grid_x, grid_y, grid_z), vp);
                                                }
                                            }
                                        }
                                    });
        }
    }

    {
        ScopeTimer timer("Test velocity transfer using FastVec3 + FastMat3");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto mp       = mass[i];
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto ppos     = FastVec3<Real>(positions[i]);
                                        const auto pvel     = FastVec3<Real>(velocities[i]);
                                        const auto pC       = FastMat3<Real>(C[i]);
                                        const auto& lcorner = kernelCornerNode[i];
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xixp = FastVec3<Real>(grid.getWorldCoordinate(grid_x, grid_y, grid_z)) - ppos;
                                                    const auto vp   = (pvel + pC * xixp) * (w * mp);
                                                    AtomicOperations::add(gridData(grid_x, grid_y, grid_z), vp.v3);
                                                }
                                            }
                                        }
                                    });
        }
    }

    {
        ScopeTimer timer("Test velocity transfer using FastVec3 + FastMat3 + FastGrid");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto mp       = mass[i];
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto ppos     = FastVec3<Real>(positions[i]);
                                        const auto pvel     = FastVec3<Real>(velocities[i]);
                                        const auto pC       = FastMat3<Real>(C[i]);
                                        const auto& lcorner = kernelCornerNode[i];
                                        for(int z = 0; z < 4; ++z) {
                                            for(int y = 0; y < 4; ++y) {
                                                for(int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xixp = fastGrid.getWorldCoordinate(grid_x, grid_y, grid_z) - ppos;
                                                    const auto vp   = (pvel + pC * xixp) * (w * mp);
                                                    AtomicOperations::add(gridData(grid_x, grid_y, grid_z), vp.v3);
                                                }
                                            }
                                        }
                                    });
        }
    }

    {
        ScopeTimer timer("Test position interpolation using GLM");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto& lcorner = kernelCornerNode[i];
                                        auto ppos           = Vec3<Real>(0);
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xi = grid.getWorldCoordinate(grid_x, grid_y, grid_z);
                                                    const auto vi = gridData(grid_x, grid_y, grid_z);
                                                    ppos         += w * (xi + Real(1e-5) * vi);
                                                }
                                            }
                                        }
                                        positions[i] = ppos;
                                    });
        }
    }

    {
        ScopeTimer timer("Test position interpolation using FastVec3 + FastMat3");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto& lcorner = kernelCornerNode[i];
                                        auto ppos           = FastVec3<Real>(0);
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xi = FastVec3<Real>(grid.getWorldCoordinate(grid_x, grid_y, grid_z));
                                                    const auto vi = FastVec3<Real>(gridData(grid_x, grid_y, grid_z));
                                                    ppos         += w * (xi + Real(1e-5) * vi);
                                                }
                                            }
                                        }
                                        positions[i] = ppos;
                                    });
        }
    }

    {
        ScopeTimer timer("Test position interpolation using FastVec3 + FastMat3 + FastGrid");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto& lcorner = kernelCornerNode[i];
                                        auto ppos           = FastVec3<Real>(0);
                                        auto gvel           = FastVec3<Real>();
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xi = fastGrid.getWorldCoordinate(grid_x, grid_y, grid_z);
                                                    gvel.v3       = gridData(grid_x, grid_y, grid_z);
                                                    ppos         += w * (xi + Real(1e-5) * gvel);
                                                }
                                            }
                                        }
                                        positions[i] = ppos.toVec3r();
                                    });
        }
    }

    {
        ScopeTimer timer("Test position interpolation using FastVec3 + FastMat3 + FastGrid");
        for(int i = 0; i < PERFORMANCE_TEST_NUM; ++i) {
            Scheduler::parallel_for(DATA_SIZE,
                                    [&](int i) {
                                        const auto w        = NumberHelpers::fRand<Real>::rnd();
                                        const auto& lcorner = kernelCornerNode[i];
                                        auto ppos           = Vec3<Real>(0);
                                        for(Int z = 0; z < 4; ++z) {
                                            for(Int y = 0; y < 4; ++y) {
                                                for(Int x = 0; x < 4; ++x) {
                                                    const auto grid_x = lcorner.x + x;
                                                    const auto grid_y = lcorner.y + y;
                                                    const auto grid_z = lcorner.z + z;

                                                    const auto xi = fastGrid.getWorldCoordinate(grid_x, grid_y, grid_z).v3;
                                                    const auto vi = gridData(grid_x, grid_y, grid_z);
                                                    ppos         += w * (xi + Real(1e-5) * vi);
                                                }
                                            }
                                        }
                                        positions[i] = ppos;
                                    });
        }
    }
}

#endif
