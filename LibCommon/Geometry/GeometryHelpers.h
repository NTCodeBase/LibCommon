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

#include <LibCommon/Array/Array.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace GeometryHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t point_line_distance(const VecX<N, Real_t>& p, const VecX<N, Real_t>& x0, const VecX<N, Real_t>& x1);

template<class Real_t>
Real_t point_segment_distance(const Vec3<Real_t>& x0, const Vec3<Real_t>& x1, const Vec3<Real_t>& x2);

template<class Real_t>
Real_t point_triangle_distance(const Vec3<Real_t>& x0, const Vec3<Real_t>& x1, const Vec3<Real_t>& x2, const Vec3<Real_t>& x3);

template<class Real_t>
void check_neighbour(const StdVT<Vec3ui>& tri, const StdVT_Vec3<Real_t>& x, Array<3, Real_t>& phi, Array3ui& closest_tri,
                     const Vec3<Real_t>& gx,
                     Int i0, Int j0, Int k0,
                     Int i1, Int j1, Int k1);

template<class Real_t>
void sweep(const StdVT<Vec3ui>& tri, const StdVT_Vec3<Real_t>& x,
           Array<3, Real_t>& phi, Array3ui& closest_tri, const Vec3<Real_t>& origin, Real_t dx,
           Int di, Int dj, Int dk);

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
template<class Real_t>
Int orientation(Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t& twice_signed_area);

// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
template<class Real_t>
bool point_in_triangle_2d(Real_t x0, Real_t y0,
                          Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t x3, Real_t y3,
                          Real_t& a, Real_t& b, Real_t& c);

template<class Real_t>
bool point_in_triangle_2d(Real_t x0, Real_t y0, Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t x3, Real_t y3);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace GeometryHelpers
