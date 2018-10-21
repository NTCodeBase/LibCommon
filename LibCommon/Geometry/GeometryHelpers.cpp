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

#include <LibCommon/Geometry/GeometryHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace GeometryHelpers
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
Real_t point_line_distance(const VecX<N, Real_t>& p, const VecX<N, Real_t>& x0, const VecX<N, Real_t>& x1)
{
    auto x01 = x1 - x0;
    auto x0p = p - x0;
    auto prj = glm::dot(x01, x0p) * glm::normalize(x01);
    return glm::length(x0p - prj);
}

template float  point_line_distance<2, float>(const VecX<2, float>& p, const VecX<2, float>& x0, const VecX<2, float>& x1);
template double point_line_distance<2, double>(const VecX<2, double>& p, const VecX<2, double>& x0, const VecX<2, double>& x1);

template float  point_line_distance<3, float>(const VecX<3, float>& p, const VecX<3, float>& x0, const VecX<3, float>& x1);
template double point_line_distance<3, double>(const VecX<3, double>& p, const VecX<3, double>& x0, const VecX<3, double>& x1);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// find distance x0 is from segment x1-x2
template<class Real_t>
Real_t point_segment_distance(const Vec3<Real_t>& x0, const Vec3<Real_t>& x1, const Vec3<Real_t>& x2)
{
    Vec3<Real_t> dx(x2 - x1);

    Real_t m2 = glm::length2(dx);
    // find parameter value of closest point on segment
    Real_t s12 = glm::dot(dx, x2 - x0) / m2;

    if(s12 < 0) {
        s12 = 0;
    } else if(s12 > 1) {
        s12 = 1;
    }

    return glm::length(x0 - s12 * x1 + (1 - s12) * x2);
}

template float  point_segment_distance<float>(const Vec3<float>& x0, const Vec3<float>& x1, const Vec3<float>& x2);
template double point_segment_distance<double>(const Vec3<double>& x0, const Vec3<double>& x1, const Vec3<double>& x2);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// find distance x0 is from triangle x1-x2-x3
template<class Real_t>
Real_t point_triangle_distance(const Vec3<Real_t>& x0, const Vec3<Real_t>& x1, const Vec3<Real_t>& x2, const Vec3<Real_t>& x3)
{
    // first find barycentric coordinates of closest point on infinite plane
    Vec3<Real_t> x13(x1 - x3), x23(x2 - x3), x03(x0 - x3);
    Real_t       m13 = glm::length2(x13), m23 = glm::length2(x23), d = glm::dot(x13, x23);

    Real_t invdet = 1.f / fmax(m13 * m23 - d * d, 1e-30f);
    Real_t a = glm::dot(x13, x03), b = glm::dot(x23, x03);

    // the barycentric coordinates themselves
    Real_t w23 = invdet * (m23 * a - d * b);
    Real_t w31 = invdet * (m13 * b - d * a);
    Real_t w12 = 1 - w23 - w31;

    if(w23 >= 0 && w31 >= 0 && w12 >= 0) {     // if we're inside the triangle
        return glm::length(x0 - w23 * x1 + w31 * x2 + w12 * x3);
    } else {                                   // we have to clamp to one of the edges
        if(w23 > 0) {                          // this rules out edge 2-3 for us
            return std::min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x1, x3));
        } else if(w31 > 0) {                   // this rules out edge 1-3
            return std::min(point_segment_distance(x0, x1, x2), point_segment_distance(x0, x2, x3));
        } else {                               // w12 must be >0, ruling out edge 1-2
            return std::min(point_segment_distance(x0, x1, x3), point_segment_distance(x0, x2, x3));
        }
    }
}

template float  point_triangle_distance<float>(const Vec3<float>& x0, const Vec3<float>& x1, const Vec3<float>& x2, const Vec3<float>& x3);
template double point_triangle_distance<double>(const Vec3<double>& x0, const Vec3<double>& x1, const Vec3<double>& x2, const Vec3<double>& x3);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void check_neighbour(const StdVT<Vec3ui>& tri, const StdVT_Vec3<Real_t>& x, Array<3, Real_t>& phi, Array3ui& closest_tri,
                     const Vec3<Real_t>& gx,
                     Int i0, Int j0, Int k0,
                     Int i1, Int j1, Int k1)
{
    if(closest_tri(i1, j1, k1) != 0xffffffff) {
        UInt p = tri[closest_tri(i1, j1, k1)][0];
        UInt q = tri[closest_tri(i1, j1, k1)][1];
        UInt r = tri[closest_tri(i1, j1, k1)][2];

        Real_t d = point_triangle_distance(gx, x[p], x[q], x[r]);

        if(d < phi(i0, j0, k0)) {
            phi(i0, j0, k0)         = Real_t(d);
            closest_tri(i0, j0, k0) = closest_tri(i1, j1, k1);
        }
    }
}

template void check_neighbour<float>(const StdVT<Vec3ui>& tri, const StdVT_Vec3<float>& x, Array<3, float>& phi, Array3ui& closest_tri,
                                     const Vec3<float>& gx,
                                     Int i0, Int j0, Int k0,
                                     Int i1, Int j1, Int k1);
template void check_neighbour<double>(const StdVT<Vec3ui>& tri, const StdVT_Vec3<double>& x, Array<3, double>& phi, Array3ui& closest_tri,
                                      const Vec3<double>& gx,
                                      Int i0, Int j0, Int k0,
                                      Int i1, Int j1, Int k1);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

template<class Real_t>
void sweep(const StdVT<Vec3ui>& tri, const StdVT_Vec3<Real_t>& x,
           Array<3, Real_t>& phi, Array3ui& closest_tri, const Vec3<Real_t>& origin, Real_t dx,
           Int di, Int dj, Int dk)
{
    Int i0, i1;
    Int j0, j1;
    Int k0, k1;

    if(di > 0) {
        i0 = 1;
        i1 = static_cast<Int>(phi.vsize()[0]);
    } else {
        i0 = static_cast<Int>(phi.vsize()[0]) - 2;
        i1 = -1;
    }

    if(dj > 0) {
        j0 = 1;
        j1 = static_cast<Int>(phi.vsize()[1]);
    } else {
        j0 = static_cast<Int>(phi.vsize()[1]) - 2;
        j1 = -1;
    }

    if(dk > 0) {
        k0 = 1;
        k1 = static_cast<Int>(phi.vsize()[2]);
    } else {
        k0 = static_cast<Int>(phi.vsize()[2]) - 2;
        k1 = -1;
    }

    //    Scheduler::parallel_for<Int>(i0, i1 + 1, j0, j1 + 1, k0, k1 + 1,
    //                                       [&](Int i, Int j, Int k)

    for(Int k = k0; k != k1; k += dk) {
        for(Int j = j0; j != j1; j += dj) {
            for(Int i = i0; i != i1; i += di) {
                Vec3<Real_t> gx = Vec3<Real_t>(i, j, k) * dx + origin;

                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j,      k);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,      j - dj, k);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj, k);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,      j,      k - dk);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j,      k - dk);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i,      j - dj, k - dk);
                check_neighbour(tri, x, phi, closest_tri, gx, i, j, k, i - di, j - dj, k - dk);
            }
        }
    }
}

template void sweep<float>(const StdVT<Vec3ui>& tri,
                           const StdVT_Vec3<float>& x,
                           Array<3, float>& phi, Array3ui& closest_tri, const Vec3<float>& origin, float dx,
                           Int di, Int dj, Int dk);
template void sweep<double>(const StdVT<Vec3ui>& tri,
                            const StdVT_Vec3<double>& x,
                            Array<3, double>& phi, Array3ui& closest_tri, const Vec3<double>& origin, double dx,
                            Int di, Int dj, Int dk);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

// calculate twice signed area of triangle (0,0)-(x1,y1)-(x2,y2)
// return an SOS-determined sign (-1, +1, or 0 only if it's a truly degenerate triangle)
template<class Real_t>
Int orientation(Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t& twice_signed_area)
{
    twice_signed_area = y1 * x2 - x1 * y2;

    if(twice_signed_area > 0) {
        return 1;
    } else if(twice_signed_area < 0) {
        return -1;
    } else if(y2 > y1) {
        return 1;
    } else if(y2 < y1) {
        return -1;
    } else if(x1 > x2) {
        return 1;
    } else if(x1 < x2) {
        return -1;
    } else {
        return 0;                    // only true when x1==x2 and y1==y2
    }
}

template Int orientation<float>(float x1, float y1, float x2, float y2, float& twice_signed_area);
template Int orientation<double>(double x1, double y1, double x2, double y2, double& twice_signed_area);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// robust test of (x0,y0) in the triangle (x1,y1)-(x2,y2)-(x3,y3)
// if true is returned, the barycentric coordinates are set in a,b,c.
template<class Real_t>
bool point_in_triangle_2d(Real_t x0, Real_t y0,
                          Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t x3, Real_t y3,
                          Real_t& a, Real_t& b, Real_t& c)
{
    x1 -= x0;
    x2 -= x0;
    x3 -= x0;
    y1 -= y0;
    y2 -= y0;
    y3 -= y0;
    Int signa = orientation(x2, y2, x3, y3, a);

    if(signa == 0) {
        return false;
    }

    Int signb = orientation(x3, y3, x1, y1, b);

    if(signb != signa) {
        return false;
    }

    Int signc = orientation(x1, y1, x2, y2, c);

    if(signc != signa) {
        return false;
    }

    Real_t sum = a + b + c;
    __NT_REQUIRE(sum != 0);                                 // if the SOS signs match and are nonkero, there's no way all of a, b, and c are zero.
    a /= sum;
    b /= sum;
    c /= sum;
    return true;
}

template<class Real_t>
bool point_in_triangle_2d(Real_t x0, Real_t y0, Real_t x1, Real_t y1, Real_t x2, Real_t y2, Real_t x3, Real_t y3)
{
    x1 -= x0;
    x2 -= x0;
    x3 -= x0;
    y1 -= y0;
    y2 -= y0;
    y3 -= y0;

    Real_t a;
    Int      signa = orientation(x2, y2, x3, y3, a);
    if(signa == 0) {
        return false;
    }

    Real_t b;
    Int      signb = orientation(x3, y3, x1, y1, b);
    if(signb != signa) {
        return false;
    }

    Real_t c;
    Int      signc = orientation(x1, y1, x2, y2, c);
    if(signc != signa) {
        return false;
    }

    return true;
}

template bool point_in_triangle_2d<float>(float x0, float y0,
                                          float x1, float y1, float x2, float y2, float x3, float y3,
                                          float& a, float& b, float& c);
template bool point_in_triangle_2d<double>(double x0, double y0,
                                           double x1, double y1, double x2, double y2, double x3, double y3,
                                           double& a, double& b, double& c);

template bool point_in_triangle_2d<float>(float x0, float y0, float x1, float y1, float x2, float y2, float x3, float y3);
template bool point_in_triangle_2d<double>(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace GeometryHelpers
