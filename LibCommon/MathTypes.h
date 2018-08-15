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
#include <BasicTypes.h>

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

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
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
