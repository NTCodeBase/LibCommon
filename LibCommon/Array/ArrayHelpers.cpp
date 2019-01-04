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

#include <LibCommon/Array/ArrayHelpers.h>
#include <LibCommon/Math/MathHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::ArrayHelpers {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void getCoordinatesAndWeights(const Vec2<Real_t>& point, const Vec2ui& size, std::array<Vec2i, 8>& indices, std::array<Real_t, 8>& weights) {
    Int    i, j;
    Real_t fi, fj;

    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(size[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(size[1]));

    indices[0] = Vec2i(i, j);
    indices[1] = Vec2i(i + 1, j);
    indices[2] = Vec2i(i, j + 1);
    indices[3] = Vec2i(i + 1, j + 1);

    weights[0] = (Real_t(1.0) - fi) * (Real_t(1.0) - fj);
    weights[1] = fi * (Real_t(1.0) - fj);
    weights[2] = (Real_t(1.0) - fi) * fj;
    weights[3] = fi * fj;
}

template<class Real_t>
void getCoordinatesAndWeights(const Vec3<Real_t>& point, const Vec3ui& size, std::array<Vec3i, 8>& indices, std::array<Real_t, 8>& weights) {
    Int    i, j, k;
    Real_t fi, fj, fk;

    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(size[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(size[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(size[2]));

    indices[0] = Vec3i(i, j, k);
    indices[1] = Vec3i(i + 1, j, k);
    indices[2] = Vec3i(i, j + 1, k);
    indices[3] = Vec3i(i + 1, j + 1, k);
    indices[4] = Vec3i(i, j, k + 1);
    indices[5] = Vec3i(i + 1, j, k + 1);
    indices[6] = Vec3i(i, j + 1, k + 1);
    indices[7] = Vec3i(i + 1, j + 1, k + 1);

    weights[0] = (Real_t(1.0) - fi) * (Real_t(1.0) - fj) * (Real_t(1.0) - fk);
    weights[1] = fi * (Real_t(1.0) - fj) * (Real_t(1.0) - fk);
    weights[2] = (Real_t(1.0) - fi) * fj * (Real_t(1.0) - fk);
    weights[3] = fi * fj * (Real_t(1.0) - fk);
    weights[4] = (Real_t(1.0) - fi) * (Real_t(1.0) - fj) * fk;
    weights[5] = fi * (Real_t(1.0) - fj) * fk;
    weights[6] = (Real_t(1.0) - fi) * fj * fk;
    weights[7] = fi * fj * fk;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Real_t interpolateValueLinear(const Vec2<Real_t>& point, const Array2<Real_t>& grid) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    return MathHelpers::bilerp(grid(i, j), grid(i + 1, j), grid(i, j + 1), grid(i + 1, j + 1), fi, fj);
}

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
Real_t interpolateValueLinear(const Vec3<Real_t>& point, const Array3<Real_t>& grid) {
    Int    i, j, k;
    Real_t fi, fj, fk;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(grid.vsize()[2]));
    return MathHelpers::trilerp(
        grid(i, j, k), grid(i + 1, j, k), grid(i, j + 1, k), grid(i + 1, j + 1, k),
        grid(i, j, k + 1), grid(i + 1, j, k + 1), grid(i, j + 1, k + 1), grid(i + 1, j + 1, k + 1),
        fi, fj, fk);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Vec2<Real_t> grad_bilerp(Real_t v00, Real_t v10, Real_t v01, Real_t v11, Real_t fi, Real_t fj) {
    return Vec2<Real_t>(fj - Real_t(1.0), fi - Real_t(1.0)) * v00 +
           Vec2<Real_t>(Real_t(1.0) - fj, -fi) * v10 +
           Vec2<Real_t>(-fj,              Real_t(1.0) - fi) * v01 +
           Vec2<Real_t>(fj,               fi) * v11;
}

template<class Real_t>
Vec2<Real_t> interpolateValueAffine(const Vec2<Real_t>& point, const Array2<Real_t>& grid) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    return grad_bilerp(
        grid(i, j), grid(i + 1, j),
        grid(i, j + 1), grid(i + 1, j + 1),
        fi, fj);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Real_t interpolateValueCubicBSpline(const Vec2<Real_t>& point, const Array2<Real_t>& grid) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));

    Real_t sumW   = 0;
    Real_t sumVal = 0;
    for(Int lj = -1; lj <= 2; ++lj) {
        for(Int li = -1; li <= 2; ++li) {
            const Vec2i ind = Vec2i(i + li, j + lj);
            if(grid.isValidIndex(ind)) {
                const Real_t weight = MathHelpers::cubic_bspline_2d(fi - static_cast<Real_t>(li), fj - static_cast<Real_t>(lj));
                sumW   += weight;
                sumVal += weight * grid(ind);
            }
        }
    }
    if(sumW > 0) {
        return sumVal / sumW;
    } else {
        return 0;
    }
}

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
Real_t interpolateValueCubicBSpline(const Vec3<Real_t>& point, const Array3<Real_t>& grid) {
    Int    i, j, k;
    Real_t fi, fj, fk;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(grid.vsize()[2]));

    Real_t sumW   = 0;
    Real_t sumVal = 0;
    for(Int lk = -1; lk <= 2; ++lk) {
        for(Int lj = -1; lj <= 2; ++lj) {
            for(Int li = -1; li <= 2; ++li) {
                const Vec3i ind = Vec3i(i + li, j + lj, k + lk);
                if(grid.isValidIndex(ind)) {
                    const Real_t weight = MathHelpers::cubic_bspline_3d(fi - static_cast<Real_t>(li),
                                                                        fj - static_cast<Real_t>(lj),
                                                                        fk - static_cast<Real_t>(lk));
                    sumW   += weight;
                    sumVal += weight * grid(ind);
                }
            }
        }
    }
    if(sumW > 0) {
        return sumVal / sumW;
    } else {
        return 0;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Vec2<Real_t> interpolateGradient(const Vec2<Real_t>& point, const Array2<Real_t>& grid) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));

    Real_t v00 = grid(i, j);
    Real_t v01 = grid(i, j + 1);
    Real_t v10 = grid(i + 1, j);
    Real_t v11 = grid(i + 1, j + 1);

    Real_t ddy0 = (v01 - v00);
    Real_t ddy1 = (v11 - v10);

    Real_t ddx0 = (v10 - v00);
    Real_t ddx1 = (v11 - v01);

    return Vec2<Real_t>(MathHelpers::lerp(ddx0, ddx1, fj),
                        MathHelpers::lerp(ddy0, ddy1, fi));
}

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
Vec3<Real_t> interpolateGradient(const Vec3<Real_t>& point, const Array3<Real_t>& grid) {
    Int    i, j, k;
    Real_t fi, fj, fk;

    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(grid.vsize()[2]));

    Real_t v000 = grid(i, j, k);
    Real_t v001 = grid(i, j, k + 1);
    Real_t v010 = grid(i, j + 1, k);
    Real_t v011 = grid(i, j + 1, k + 1);
    Real_t v100 = grid(i + 1, j, k);
    Real_t v101 = grid(i + 1, j, k + 1);
    Real_t v110 = grid(i + 1, j + 1, k);
    Real_t v111 = grid(i + 1, j + 1, k + 1);

    Real_t ddx00 = (v100 - v000);
    Real_t ddx10 = (v110 - v010);
    Real_t ddx01 = (v101 - v001);
    Real_t ddx11 = (v111 - v011);
    Real_t dv_dx = MathHelpers::bilerp(ddx00, ddx10, ddx01, ddx11, fj, fk);

    Real_t ddy00 = (v010 - v000);
    Real_t ddy10 = (v110 - v100);
    Real_t ddy01 = (v011 - v001);
    Real_t ddy11 = (v111 - v101);
    Real_t dv_dy = MathHelpers::bilerp(ddy00, ddy10, ddy01, ddy11, fi, fk);

    Real_t ddz00 = (v001 - v000);
    Real_t ddz10 = (v101 - v100);
    Real_t ddz01 = (v011 - v010);
    Real_t ddz11 = (v111 - v110);
    Real_t dv_dz = MathHelpers::bilerp(ddz00, ddz10, ddz01, ddz11, fi, fj);

    return Vec3<Real_t>(dv_dx, dv_dy, dv_dz);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Vec2<Real_t> interpolateGradientValue(const Vec2<Real_t>& point, const Array2<Real_t>& grid, Real_t cellSize) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));

    Real_t v00 = grid(i, j);
    Real_t v01 = grid(i, j + 1);
    Real_t v10 = grid(i + 1, j);
    Real_t v11 = grid(i + 1, j + 1);

    Vec2<Real_t> grad = v00 * Vec2<Real_t>(fj - Real_t(1.0), fi - Real_t(1.0)) +
                        v10 * Vec2<Real_t>(Real_t(1.0) - fj, -fi) +
                        v01 * Vec2<Real_t>(-fj, Real_t(1.0) - fi) +
                        v11 * Vec2<Real_t>(fj, fi);
    return grad / cellSize;
}

template<class Real_t>
Vec3<Real_t> interpolateGradientValue(const Vec3<Real_t>& point, const Array3<Real_t>& grid, Real_t cellSize) {
    Int    i, j, k;
    Real_t fi, fj, fk;

    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(grid.vsize()[2]));

    Real_t v000 = grid(i, j, k);
    Real_t v001 = grid(i, j, k + 1);
    Real_t v010 = grid(i, j + 1, k);
    Real_t v011 = grid(i, j + 1, k + 1);
    Real_t v100 = grid(i + 1, j, k);
    Real_t v101 = grid(i + 1, j, k + 1);
    Real_t v110 = grid(i + 1, j + 1, k);
    Real_t v111 = grid(i + 1, j + 1, k + 1);

    Vec3<Real_t> grad = v000 * Vec3<Real_t>(-(Real_t(1.0) - fj) * (Real_t(1.0) - fk), -(Real_t(1.0) - fi) * (Real_t(1.0) - fk), -(Real_t(1.0) - fi) * (Real_t(1.0) - fj)) +
                        v001 * Vec3<Real_t>(-(Real_t(1.0) - fj) * fk, -(Real_t(1.0) - fi) * fk, (Real_t(1.0) - fi) * (Real_t(1.0) - fj)) +
                        v010 * Vec3<Real_t>(-fj * (Real_t(1.0) - fk), (Real_t(1.0) - fi) * (Real_t(1.0) - fk), -(Real_t(1.0) - fi) * fj) +
                        v011 * Vec3<Real_t>(-fj * fk, (Real_t(1.0) - fi) * fk, (Real_t(1.0) - fi) * fj) +
                        v100 * Vec3<Real_t>((Real_t(1.0) - fj) * (Real_t(1.0) - fk), -fi * (Real_t(1.0) - fk), -fi * (Real_t(1.0) - fj)) +
                        v101 * Vec3<Real_t>((Real_t(1.0) - fj) * fk, -fi * fk, fi * (Real_t(1.0) - fj)) +
                        v110 * Vec3<Real_t>(fj * (Real_t(1.0) - fk), fi * (Real_t(1.0) - fk), -fi * fj) +
                        v111 * Vec3<Real_t>(fj * fk, fi * fk, fi * fj);
    return grad / cellSize;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
Real_t interpolateValueAndGradient(Vec2<Real_t>& gradient, const Vec2<Real_t>& point, const Array2<Real_t>& grid) {
    Int    i, j;
    Real_t fi, fj;
    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));

    Real_t v00 = grid(i, j);
    Real_t v01 = grid(i, j + 1);
    Real_t v10 = grid(i + 1, j);
    Real_t v11 = grid(i + 1, j + 1);

    Real_t ddy0 = (v01 - v00);
    Real_t ddy1 = (v11 - v10);

    Real_t ddx0 = (v10 - v00);
    Real_t ddx1 = (v11 - v01);

    gradient[0] = MathHelpers::lerp(ddx0, ddx1, fj);
    gradient[1] = MathHelpers::lerp(ddy0, ddy1, fi);

    //may as well return value too
    return MathHelpers::bilerp(v00, v10, v01, v11, fi, fj);
}

template<class Real_t>
Real_t interpolateValueAndGradient(Vec3<Real_t>& gradient, const Vec3<Real_t>& point, const Array3<Real_t>& grid) {
    Int    i, j, k;
    Real_t fi, fj, fk;

    MathHelpers::get_barycentric(point[0], i, fi, 0, static_cast<Int>(grid.vsize()[0]));
    MathHelpers::get_barycentric(point[1], j, fj, 0, static_cast<Int>(grid.vsize()[1]));
    MathHelpers::get_barycentric(point[2], k, fk, 0, static_cast<Int>(grid.vsize()[2]));

    Real_t v000 = grid(i, j, k);
    Real_t v001 = grid(i, j, k + 1);
    Real_t v010 = grid(i, j + 1, k);
    Real_t v011 = grid(i, j + 1, k + 1);
    Real_t v100 = grid(i + 1, j, k);
    Real_t v101 = grid(i + 1, j, k + 1);
    Real_t v110 = grid(i + 1, j + 1, k);
    Real_t v111 = grid(i + 1, j + 1, k + 1);

    Real_t ddx00 = (v100 - v000);
    Real_t ddx10 = (v110 - v010);
    Real_t ddx01 = (v101 - v001);
    Real_t ddx11 = (v111 - v011);
    Real_t dv_dx = MathHelpers::bilerp(ddx00, ddx10, ddx01, ddx11, fj, fk);

    Real_t ddy00 = (v010 - v000);
    Real_t ddy10 = (v110 - v100);
    Real_t ddy01 = (v011 - v001);
    Real_t ddy11 = (v111 - v101);
    Real_t dv_dy = MathHelpers::bilerp(ddy00, ddy10, ddy01, ddy11, fi, fk);

    Real_t ddz00 = (v001 - v000);
    Real_t ddz10 = (v101 - v100);
    Real_t ddz01 = (v011 - v010);
    Real_t ddz11 = (v111 - v110);
    Real_t dv_dz = MathHelpers::bilerp(ddz00, ddz10, ddz01, ddz11, fi, fj);

    gradient[0] = dv_dx;
    gradient[1] = dv_dy;
    gradient[2] = dv_dz;

    //return value for good measure.
    return MathHelpers::trilerp(v000, v100,
                                v010, v110,
                                v001, v101,
                                v011, v111,
                                fi, fj, fk);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template void getCoordinatesAndWeights<float>(const Vec2<float>& point, const Vec2ui& size, std::array<Vec2i, 8>& indices, std::array<float, 8>& weights);
template void getCoordinatesAndWeights<float>(const Vec3<float>& point, const Vec3ui& size, std::array<Vec3i, 8>& indices, std::array<float, 8>& weights);

template float interpolateValueLinear<float>(const Vec2<float>& point, const Array2<float>& grid);
template float interpolateValueLinear<float>(const Vec3<float>& point, const Array3<float>& grid);

template Vec2<float> grad_bilerp<float>(float v00, float v10, float v01, float v11, float fi, float fj);

template Vec2<float> interpolateValueAffine<float>(const Vec2<float>& point, const Array2<float>& grid);

template float interpolateValueCubicBSpline<float>(const Vec2<float>& point, const Array2<float>& grid);
template float interpolateValueCubicBSpline<float>(const Vec3<float>& point, const Array3<float>& grid);

template Vec2<float> interpolateGradient<float>(const Vec2<float>& point, const Array2<float>& grid);
template Vec3<float> interpolateGradient<float>(const Vec3<float>& point, const Array3<float>& grid);

template Vec2<float> interpolateGradientValue<float>(const Vec2<float>& point, const Array2<float>& grid, float cellSize);
template Vec3<float> interpolateGradientValue<float>(const Vec3<float>& point, const Array3<float>& grid, float cellSize);

template float interpolateValueAndGradient<float>(Vec2<float>& gradient, const Vec2<float>& point, const Array2<float>& grid);
template float interpolateValueAndGradient<float>(Vec3<float>& gradient, const Vec3<float>& point, const Array3<float>& grid);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template void getCoordinatesAndWeights<double>(const Vec2<double>& point, const Vec2ui& size, std::array<Vec2i, 8>& indices, std::array<double, 8>& weights);
template void getCoordinatesAndWeights<double>(const Vec3<double>& point, const Vec3ui& size, std::array<Vec3i, 8>& indices, std::array<double, 8>& weights);

template double interpolateValueLinear<double>(const Vec2<double>& point, const Array2<double>& grid);
template double interpolateValueLinear<double>(const Vec3<double>& point, const Array3<double>& grid);

template Vec2<double> grad_bilerp<double>(double v00, double v10, double v01, double v11, double fi, double fj);

template Vec2<double> interpolateValueAffine<double>(const Vec2<double>& point, const Array2<double>& grid);

template double interpolateValueCubicBSpline<double>(const Vec2<double>& point, const Array2<double>& grid);
template double interpolateValueCubicBSpline<double>(const Vec3<double>& point, const Array3<double>& grid);

template Vec2<double> interpolateGradient<double>(const Vec2<double>& point, const Array2<double>& grid);
template Vec3<double> interpolateGradient<double>(const Vec3<double>& point, const Array3<double>& grid);

template Vec2<double> interpolateGradientValue<double>(const Vec2<double>& point, const Array2<double>& grid, double cellSize);
template Vec3<double> interpolateGradientValue<double>(const Vec3<double>& point, const Array3<double>& grid, double cellSize);

template double interpolateValueAndGradient<double>(Vec2<double>& gradient, const Vec2<double>& point, const Array2<double>& grid);
template double interpolateValueAndGradient<double>(Vec3<double>& gradient, const Vec3<double>& point, const Array3<double>& grid);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::ArrayHelpers
