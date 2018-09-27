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

#define PERFORMANCE_TEST_NUM 10000
#define DATA_SIZE            1000
//#define DATA_SIZE            100000000

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

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
TEST_CASE("Test_Correctness_Initialization", "Test_Correctness_Initialization")
{
    StdVT<Vec3<Real>>     v3_data;
    StdVT<FastVec3<Real>> fv3_data;

    StdVT<Mat3x3<Real>>   m3_data;
    StdVT<FastMat3<Real>> fm3_data;
    ////////////////////////////////////////////////////////////////////////////////
    {
        auto check_vec3 = [](const auto& v3, const auto& fv3) {
                              volatile bool bResult =
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
        ////////////////////////////////////////////////////////////////////////////////
        init_glmVec3(v3_data);
        for(UInt i = 0; i < DATA_SIZE; ++i) {
            auto& v3 = v3_data[i];
            {
                FastVec3<Real> fv3 = v3;
                REQUIRE(check_vec3(v3, fv3));
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
        init_FastMat3(fm3_data);
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
}
