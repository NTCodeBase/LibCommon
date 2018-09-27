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

#define PRECISION            1e-7

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
                REQUIRE(glm::length2(v3 - v3_copy) < PRECISION);
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
                REQUIRE(glm::length2(m3[0] - m3_copy1[0]) + glm::length2(m3[1] - m3_copy1[1]) + glm::length2(m3[2] - m3_copy1[2]) < PRECISION);

                Mat3x3<Real> m3_copy2;
                fm3.copyToMat3x3r(m3_copy2);
                REQUIRE(check_mat3(m3_copy2, fm3));
                REQUIRE(glm::length2(m3[0] - m3_copy2[0]) + glm::length2(m3[1] - m3_copy2[1]) + glm::length2(m3[2] - m3_copy2[2]) < PRECISION);
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

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
TEST_CASE("Test_FastVec3_Ops", "Test_FastVec3_Ops")
{
    {
        for(Int i = 0; i < DATA_SIZE; ++i) {
            Vec3<Real>     v0  = NumberHelpers::fRand01<Real>::vrnd<Vec3<Real>>();
            Vec3<Real>     v1  = NumberHelpers::fRand01<Real>::vrnd<Vec3<Real>>();
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
                auto r = NumberHelpers::fRand01<Real>::rnd();
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
        }
    }
}
