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

#include <LibCommon/Math/FastVec3.h>
#include <LibCommon/Math/FastMat3.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
using Real = float;

//#define FAST_TEST
#define PERFORMANCE_TEST_NUM 10000

#ifdef FAST_TEST
#  define DATA_SIZE          1000
#else
#  define DATA_SIZE          100000000
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
auto init_glmVec3 = [] (auto& v3_data) {
                        v3_data.reserve(DATA_SIZE);
                        for(int i = 0; i < DATA_SIZE; ++i) {
                            v3_data.push_back(NumberHelpers::fRand01<Real>::vrnd<Vec3<Real>>());
                        }
                    };

auto init_FastVec3 = [](auto& fv3_data) {
                         fv3_data.reserve(DATA_SIZE);
                         for(int i = 0; i < DATA_SIZE; ++i) {
                             fv3_data.push_back(NumberHelpers::fRand01<Real>::vrnd<Vec3<Real>>());
                         }
                     };

auto init_glmMat3 = [] (auto& m3_data) {
                        m3_data.reserve(DATA_SIZE);
                        for(int i = 0; i < DATA_SIZE; ++i) {
                            m3_data.push_back(NumberHelpers::fRand01<Real>::vrnd<Mat3x3<Real>>());
                        }
                    };

auto init_FastMat3 = [](auto& fm3_data) {
                         fm3_data.reserve(DATA_SIZE);
                         for(int i = 0; i < DATA_SIZE; ++i) {
                             fm3_data.push_back(NumberHelpers::fRand01<Real>::vrnd<Mat3x3<Real>>());
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
                          ((std::abs(v3[0] - fv3[0]) < 1e-10) &&
                           (std::abs(v3[1] - fv3[1]) < 1e-10) &&
                           (std::abs(v3[2] - fv3[2]) < 1e-10) &&
                           (std::abs(fv3[3]) < 1e-10) &&
                           (std::abs(fv3.dummy) < 1e-10) &&
                           //
                           (std::abs(v3.x - fv3.x) < 1e-10) &&
                           (std::abs(v3.y - fv3.y) < 1e-10) &&
                           (std::abs(v3.z - fv3.z) < 1e-10) &&
                           (std::abs(fv3.w) < 1e-10)
                          );

                      if(!bResult) {
                          printf("Not correct: glm_v3 = %f, %f, %f  | FastVec3: %f, %f, %f, %f - %f, %f, %f, %f\n",
                                 v3.x, v3.y, v3.z,
                                 fv3.x, fv3.y, fv3.z, fv3.w,
                                 fv3[0], fv3[1], fv3[2], fv3[3]
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
                                 "FastVec3: \n  %f, %f, %f, %f\n  %f, %f, %f, %f\n  %f, %f, %f, %f\n",
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
TEST_CASE("Test_Correctness_Initialization", "Test_Correctness_Initialization")
{
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
                REQUIRE(glm::length2(v3 - v3_copy) < 1e-10);
            }
            {
                FastVec3<Real> fv3(v3.x, v3.y, v3.z);
                REQUIRE(check_vec3(v3, fv3));
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
                REQUIRE(glm::length2(m3[0] - m3_copy1[0]) + glm::length2(m3[1] - m3_copy1[1]) + glm::length2(m3[2] - m3_copy1[2]) < 1e-10);

                Mat3x3<Real> m3_copy2;
                fm3.copyToMat3x3r(m3_copy2);
                REQUIRE(check_mat3(m3_copy2, fm3));
                REQUIRE(glm::length2(m3[0] - m3_copy2[0]) + glm::length2(m3[1] - m3_copy2[1]) + glm::length2(m3[2] - m3_copy2[2]) < 1e-10);
            }
            {
                FastMat3<Real> fm3(m3[0], m3[1], m3[2]);
                REQUIRE(check_mat3(m3, fm3));
            }
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
TEST_CASE("Test_Performance_Initialization", "Test_Performance_Initialization")
{
    StdVT<Vec3<Real>>     v3_data;
    StdVT<FastVec3<Real>> fv3_data;

    StdVT<Mat3x3<Real>>   m3_data;
    StdVT<FastMat3<Real>> fm3_data;
    ////////////////////////////////////////////////////////////////////////////////
    {
        ScopeTimer timer("GLM vec3 initialization");
        init_glmVec3(v3_data);
    }
    {
        ScopeTimer timer("FastVec3 initialization");
        init_FastVec3(fv3_data);
    }
    {
        ScopeTimer timer("GLM mat3 initialization");
        init_glmMat3(m3_data);
    }
    {
        ScopeTimer timer("FastMat3 initialization");
        init_FastMat3(fm3_data);
    }
    {
        ScopeTimer  timer("Overwrite GLM mat3 to GLM mat3");
        const auto& m3 = m3_data.front();
        overwrite_glmMat3(m3, m3_data);
    }
    {
        ScopeTimer  timer("Overwrite FastMat3 to GLM mat3");
        const auto& fm3 = fm3_data.front();
        overwrite_FastMat3ToMat3(fm3, m3_data);
    }
}
