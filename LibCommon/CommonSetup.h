//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTGraphics                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <limits>
#include <vector>
#include <map>
#include <set>
#include <cstdint>
#include <string>
#include <memory>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define GLM_FORCE_INLINE
#define GLM_ENABLE_EXPERIMENTAL
#define GLM_FORCE_CTOR_INIT
#define GLM_SWIZZLE

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/component_wise.hpp>

#include <json.hpp>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Basic types
using Int8  = int8_t;
using Int16 = int16_t;
using Int   = int32_t;
using Int32 = int32_t;
using Int64 = int64_t;

////////////////////////////////////////////////////////////////////////////////
using UChar  = unsigned char;
using UInt8  = uint8_t;
using UInt16 = uint16_t;
using UInt   = uint32_t;
using UInt32 = uint32_t;
using UInt64 = uint64_t;

////////////////////////////////////////////////////////////////////////////////
using String  = std::string;
using JParams = nlohmann::json;

////////////////////////////////////////////////////////////////////////////////
using PairInt8  = std::pair<Int8, Int8>;
using PairInt16 = std::pair<Int16, Int16>;
using PairInt32 = std::pair<Int, Int>;
using PairInt   = std::pair<Int, Int>;
using PairInt64 = std::pair<Int64, Int64>;

////////////////////////////////////////////////////////////////////////////////
using PairUInt8  = std::pair<UInt8, UInt8>;
using PairUInt16 = std::pair<UInt16, UInt16>;
using PairUInt32 = std::pair<UInt, UInt>;
using PairUInt   = std::pair<UInt, UInt>;
using PairUInt64 = std::pair<UInt64, UInt64>;

////////////////////////////////////////////////////////////////////////////////
using Pairf = std::pair<float, float>;
using Paird = std::pair<double, double>;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// generic types
template<class Type> using SharedPtr = std::shared_ptr<Type>;
template<class Type> using UniquePtr = std::unique_ptr<Type>;

template<class K, class V> using Map = std::map<K, V>;
template<class Type> using StdVT     = std::vector<Type>;

template<class Type> using Set    = std::set<Type>;
template<class Type> using Quat   = glm::tquat<Type>;
template<class Type> using Vec2   = glm::vec<2, Type>;
template<class Type> using Vec3   = glm::vec<3, Type>;
template<class Type> using Vec4   = glm::vec<4, Type>;
template<class Type> using Mat2x2 = glm::mat<2, 2, Type>;
template<class Type> using Mat3x3 = glm::mat<3, 3, Type>;
template<class Type> using Mat4x4 = glm::mat<4, 4, Type>;

template<int N, class Type> using VecX          = glm::vec<N, Type>;
template<int N, class Type> using MatXxX        = glm::mat<N, N, Type>;
template<int M, int N, class Type> using MatMxN = glm::mat<N, M, Type>; // col major: N = cols, M = rows

template<int N, class Type> using StdVT_VecX          = StdVT<VecX<N, Type>>;
template<int N, class Type> using StdVT_MatXxX        = StdVT<MatXxX<N, Type>>;
template<int M, int N, class Type> using StdVT_MatMxN = StdVT<MatMxN<M, N, Type>>;

////////////////////////////////////////////////////////////////////////////////
template<class Type> using StdVT_Vec    = StdVT<StdVT<Type>>;
template<class Type> using StdVT_Vec2   = StdVT<Vec2<Type>>;
template<class Type> using StdVT_Vec3   = StdVT<Vec3<Type>>;
template<class Type> using StdVT_Vec4   = StdVT<Vec4<Type>>;
template<class Type> using StdVT_Mat2x2 = StdVT<Mat2x2<Type>>;
template<class Type> using StdVT_Mat3x3 = StdVT<Mat3x3<Type>>;
template<class Type> using StdVT_Mat4x4 = StdVT<Mat4x4<Type>>;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
using Quatf = Quat<float>;
using Quatd = Quat<double>;

using Vec2i = Vec2<Int>;
using Vec3i = Vec3<Int>;
using Vec4i = Vec4<Int>;

using Vec2ui = Vec2<UInt>;
using Vec3ui = Vec3<UInt>;
using Vec4ui = Vec4<UInt>;

using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Vec4f = Vec4<float>;

using Vec2d = Vec2<double>;
using Vec3d = Vec3<double>;
using Vec4d = Vec4<double>;

using Mat2x2f = Mat2x2<float>;
using Mat3x3f = Mat3x3<float>;
using Mat4x4f = Mat4x4<float>;

using Mat2x2d = Mat2x2<double>;
using Mat3x3d = Mat3x3<double>;
using Mat4x4d = Mat4x4<double>;

////////////////////////////////////////////////////////////////////////////////
using StdVT_Int8  = StdVT<Int8>;
using StdVT_Int16 = StdVT<Int16>;
using StdVT_Int   = StdVT<Int>;
using StdVT_Int32 = StdVT<Int>;
using StdVT_Int64 = StdVT<Int64>;

using StdVT_UInt8  = StdVT<UInt8>;
using StdVT_UInt16 = StdVT<UInt16>;
using StdVT_UInt   = StdVT<UInt>;
using StdVT_UInt32 = StdVT<UInt>;
using StdVT_UInt64 = StdVT<UInt64>;

using StdVT_Char   = StdVT<char>;
using StdVT_UChar  = StdVT<unsigned char>;
using StdVT_Float  = StdVT<float>;
using StdVT_Double = StdVT<double>;
using StdVT_String = StdVT<String>;

using StdVT_Vec2i = StdVT<Vec2i>;
using StdVT_Vec3i = StdVT<Vec3i>;
using StdVT_Vec4i = StdVT<Vec4i>;

using StdVT_Vec2ui = StdVT<Vec2ui>;
using StdVT_Vec3ui = StdVT<Vec3ui>;
using StdVT_Vec4ui = StdVT<Vec4ui>;

using StdVT_Vec2f = StdVT<Vec2f>;
using StdVT_Vec3f = StdVT<Vec3f>;
using StdVT_Vec4f = StdVT<Vec4f>;

using StdVT_Vec2d = StdVT<Vec2d>;
using StdVT_Vec3d = StdVT<Vec3d>;
using StdVT_Vec4d = StdVT<Vec4d>;

using StdVT_Mat2x2f = StdVT<Mat2x2f>;
using StdVT_Mat3x3f = StdVT<Mat3x3f>;
using StdVT_Mat4x4f = StdVT<Mat4x4f>;

using StdVT_Mat2x2d = StdVT<Mat2x2d>;
using StdVT_Mat3x3d = StdVT<Mat3x3d>;
using StdVT_Mat4x4d = StdVT<Mat4x4d>;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Conversion operators
float constexpr operator"" _f(long double x)
{
    return static_cast<float>(x);
}

float constexpr operator"" _f(unsigned long long int x)
{
    return static_cast<float>(x);
}

double constexpr operator"" _d(long double x)
{
    return static_cast<double>(x);
}

double constexpr operator"" _d(unsigned long long int x)
{
    return static_cast<double>(x);
}

Int32 constexpr operator"" _int(unsigned long long int x)
{
    return static_cast<Int32>(x);
}

UInt32 constexpr operator"" _uint(unsigned long long int x)
{
    return static_cast<UInt32>(x);
}

UInt64 constexpr operator"" _uint64(unsigned long long int x)
{
    return static_cast<UInt64>(x);
}

std::size_t constexpr operator"" _sz(unsigned long long int x)
{
    return static_cast<std::size_t>(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T> constexpr auto MEpsilon() { return std::numeric_limits<T>::epsilon(); }
template<class T> constexpr auto Tiny() { return std::numeric_limits<T>::min(); }
template<class T> constexpr auto Huge() { return std::numeric_limits<T>::max(); }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#include <CommonMacros.h>
#include <CommonForward.h>
