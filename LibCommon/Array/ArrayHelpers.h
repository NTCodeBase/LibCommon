//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
//    .--------------------------------------------------.
//    |  This file is part of NTCodeBase                 |
//    |  Created 2018 by NT (https://ttnghia.github.io)  |
//    '--------------------------------------------------'
//                            \o/
//                             |
//                            / |
//
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#pragma once

#include <array>
#include <LibCommon/Array/Array.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace ArrayHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t> void getCoordinatesAndWeights(const Vec2<Real_t>& point, const Vec2ui& size, std::array<Vec2i, 8>& indices, std::array<Real_t, 8>& weights);
template<class Real_t> void getCoordinatesAndWeights(const Vec3<Real_t>& point, const Vec3ui& size, std::array<Vec3i, 8>& indices, std::array<Real_t, 8>& weights);

template<class Real_t> Real_t interpolateValueLinear(const Vec2<Real_t>& point, const Array2<Real_t>& grid);
template<class Real_t> Real_t interpolateValueLinear(const Vec3<Real_t>& point, const Array3<Real_t>& grid);

template<class Real_t> Vec2<Real_t> grad_bilerp(Real_t v00, Real_t v10, Real_t v01, Real_t v11, Real_t fi, Real_t fj);

template<class Real_t> Vec2<Real_t> interpolateValueAffine(const Vec2<Real_t>& point, const Array2<Real_t>& grid);

template<class Real_t> Real_t interpolateValueCubicBSpline(const Vec2<Real_t>& point, const Array2<Real_t>& grid);
template<class Real_t> Real_t interpolateValueCubicBSpline(const Vec3<Real_t>& point, const Array3<Real_t>& grid);

template<class Real_t> Vec2<Real_t> interpolateGradient(const Vec2<Real_t>& point, const Array2<Real_t>& grid);
template<class Real_t> Vec3<Real_t> interpolateGradient(const Vec3<Real_t>& point, const Array3<Real_t>& grid);

template<class Real_t> Vec2<Real_t> interpolateGradientValue(const Vec2<Real_t>& point, const Array2<Real_t>& grid, Real_t cellSize);
template<class Real_t> Vec3<Real_t> interpolateGradientValue(const Vec3<Real_t>& point, const Array3<Real_t>& grid, Real_t cellSize);

template<class Real_t> Real_t interpolateValueAndGradient(Vec2<Real_t>& gradient, const Vec2<Real_t>& point, const Array2<Real_t>& grid);
template<class Real_t> Real_t interpolateValueAndGradient(Vec3<Real_t>& gradient, const Vec3<Real_t>& point, const Array3<Real_t>& grid);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace ArrayHelpers
