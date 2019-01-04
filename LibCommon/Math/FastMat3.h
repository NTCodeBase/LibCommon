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

#include <LibCommon/Math/FastVec3.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
class FastMat3 {
public:
    inline FastMat3() { static_assert(sizeof(T) == sizeof(float) && alignof(FastMat3<T>) == 16, "Error: Size or alignment is not correct!"); }
    inline FastMat3(__m128 c0, __m128 c1, __m128 c2) : col{{ c0 }, { c1 }, { c2 }} {}
    inline FastMat3(T x) : col{FastVec3<T>(x, T(0), T(0)), FastVec3<T>(T(0), x, T(0)), FastVec3<T>(T(0), T(0), x)} {}
    inline FastMat3(const Mat3x3<T>& m) : col{{ m[0] }, { m[1] }, { m[2] }} {}
    inline FastMat3(const FastMat3<T>& other) : col{{ other.col[0] }, { other.col[1] }, { other.col[2] }} {}
    inline FastMat3(const FastVec3<T>& c0, const FastVec3<T>& c1, const FastVec3<T>& c2) : col{{ c0 }, { c1 }, { c2 }} {}
    ////////////////////////////////////////////////////////////////////////////////
    inline FastMat3<T>& operator=(const FastMat3<T>& other) { col[0] = other.col[0]; col[1] = other.col[1]; col[2] = other.col[2]; return *this; }
    inline FastMat3<T>& operator=(const Mat3x3<T>& other) { col[0] = other[0]; col[1] = other[1]; col[2] = other[2]; return *this; }
    inline operator Mat3x3<T>() const { return Mat3x3<T>(col[0], col[1], col[2]); }
    inline void copyToMat3x3r(Mat3x3<T>& m) const { m[0] = col[0].v3; m[1] = col[1].v3; m[2] = col[2].v3; }
    ////////////////////////////////////////////////////////////////////////////////
    /// Arithmetic operators with FastMat3
    inline FastMat3<T> operator+(const FastMat3<T>& B) const {
        return FastMat3<T>(_mm_add_ps(col[0].mmvalue, B.col[0].mmvalue),
                           _mm_add_ps(col[1].mmvalue, B.col[1].mmvalue),
                           _mm_add_ps(col[2].mmvalue, B.col[2].mmvalue));
    }

    inline FastMat3<T> operator-(const FastMat3<T>& B) const {
        return FastMat3<T>(_mm_sub_ps(col[0].mmvalue, B.col[0].mmvalue),
                           _mm_sub_ps(col[1].mmvalue, B.col[1].mmvalue),
                           _mm_sub_ps(col[2].mmvalue, B.col[2].mmvalue));
    }

    inline FastVec3<T> operator*(const FastVec3<T>& v) const {
        __m128 v0 = _mm_shuffle_ps(v.mmvalue, v.mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
        __m128 v1 = _mm_shuffle_ps(v.mmvalue, v.mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
        __m128 v2 = _mm_shuffle_ps(v.mmvalue, v.mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

        __m128 m0 = _mm_mul_ps(col[0].mmvalue, v0);
        __m128 m1 = _mm_mul_ps(col[1].mmvalue, v1);
        __m128 m2 = _mm_mul_ps(col[2].mmvalue, v2);

        return _mm_add_ps(_mm_add_ps(m0, m1), m2);
    }

    inline FastMat3<T> operator*(const FastMat3<T>& B) const {
        FastMat3<T> result;
        {
            __m128 e0 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            result.col[0].mmvalue = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        {
            __m128 e0 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            result.col[1].mmvalue = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        {
            __m128 e0 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            result.col[2].mmvalue = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        ////////////////////////////////////////////////////////////////////////////////
        return result;
    }

    ////////////////////////////////////////////////////////////////////////////////
    inline FastMat3<T>& operator+=(const FastMat3<T>& B) {
        col[0].mmvalue = _mm_add_ps(col[0].mmvalue, B.col[0].mmvalue);
        col[1].mmvalue = _mm_add_ps(col[1].mmvalue, B.col[1].mmvalue);
        col[2].mmvalue = _mm_add_ps(col[2].mmvalue, B.col[2].mmvalue);
        return *this;
    }

    inline FastMat3<T>& operator-=(const FastMat3<T>& B) {
        col[0].mmvalue = _mm_sub_ps(col[0].mmvalue, B.col[0].mmvalue);
        col[1].mmvalue = _mm_sub_ps(col[1].mmvalue, B.col[1].mmvalue);
        col[2].mmvalue = _mm_sub_ps(col[2].mmvalue, B.col[2].mmvalue);
        return *this;
    }

    inline FastMat3<T>& operator*=(const FastMat3<T>& B) {
        __m128 col0m, col1m, col2m;
        {
            __m128 e0 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[0].mmvalue, B.col[0].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            col0m = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        {
            __m128 e0 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[1].mmvalue, B.col[1].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            col1m = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        {
            __m128 e0 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(0, 0, 0, 0));
            __m128 e1 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(1, 1, 1, 1));
            __m128 e2 = _mm_shuffle_ps(B.col[2].mmvalue, B.col[2].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));

            __m128 m0 = _mm_mul_ps(col[0].mmvalue, e0);
            __m128 m1 = _mm_mul_ps(col[1].mmvalue, e1);
            __m128 m2 = _mm_mul_ps(col[2].mmvalue, e2);

            col2m = _mm_add_ps(_mm_add_ps(m0, m1), m2);
        }
        col[0].mmvalue = col0m;
        col[1].mmvalue = col1m;
        col[2].mmvalue = col2m;
        ////////////////////////////////////////////////////////////////////////////////
        return *this;
    }

    ////////////////////////////////////////////////////////////////////////////////
    inline FastMat3<T> operator+(T a) const {
        return FastMat3<T>(_mm_add_ps(col[0].mmvalue, _mm_set1_ps(a)),
                           _mm_add_ps(col[1].mmvalue, _mm_set1_ps(a)),
                           _mm_add_ps(col[2].mmvalue, _mm_set1_ps(a)));
    }

    inline FastMat3<T> operator-(T a) const {
        return FastMat3<T>(_mm_sub_ps(col[0].mmvalue, _mm_set1_ps(a)),
                           _mm_sub_ps(col[1].mmvalue, _mm_set1_ps(a)),
                           _mm_sub_ps(col[2].mmvalue, _mm_set1_ps(a)));
    }

    inline FastMat3<T> operator*(T a) const {
        return FastMat3<T>(_mm_mul_ps(col[0].mmvalue, _mm_set1_ps(a)),
                           _mm_mul_ps(col[1].mmvalue, _mm_set1_ps(a)),
                           _mm_mul_ps(col[2].mmvalue, _mm_set1_ps(a)));
    }

    inline FastMat3<T> operator/(T a) const {
        assert(a != 0);
        return FastMat3<T>(_mm_div_ps(col[0].mmvalue, _mm_set1_ps(a)),
                           _mm_div_ps(col[1].mmvalue, _mm_set1_ps(a)),
                           _mm_div_ps(col[2].mmvalue, _mm_set1_ps(a)));
    }

    ////////////////////////////////////////////////////////////////////////////////
    inline FastMat3<T>& operator+=(T a) {
        col[0].mmvalue = _mm_add_ps(col[0].mmvalue, _mm_set1_ps(a));
        col[1].mmvalue = _mm_add_ps(col[1].mmvalue, _mm_set1_ps(a));
        col[2].mmvalue = _mm_add_ps(col[2].mmvalue, _mm_set1_ps(a));
        return *this;
    }

    inline FastMat3<T>& operator-=(T a) {
        col[0].mmvalue = _mm_sub_ps(col[0].mmvalue, _mm_set1_ps(a));
        col[1].mmvalue = _mm_sub_ps(col[1].mmvalue, _mm_set1_ps(a));
        col[2].mmvalue = _mm_sub_ps(col[2].mmvalue, _mm_set1_ps(a));
        return *this;
    }

    inline FastMat3<T>& operator*=(T a) {
        col[0].mmvalue = _mm_mul_ps(col[0].mmvalue, _mm_set1_ps(a));
        col[1].mmvalue = _mm_mul_ps(col[1].mmvalue, _mm_set1_ps(a));
        col[2].mmvalue = _mm_mul_ps(col[2].mmvalue, _mm_set1_ps(a));
        return *this;
    }

    inline FastMat3<T>& operator/=(T a) {
        assert(a != 0);
        col[0].mmvalue = _mm_div_ps(col[0].mmvalue, _mm_set1_ps(a));
        col[1].mmvalue = _mm_div_ps(col[1].mmvalue, _mm_set1_ps(a));
        col[2].mmvalue = _mm_div_ps(col[2].mmvalue, _mm_set1_ps(a));
        return *this;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // geometric functions
    inline FastMat3<T> transposed() const {
        __m128 t0 = _mm_shuffle_ps(col[0].mmvalue, col[1].mmvalue, _MM_SHUFFLE(1, 0, 1, 0));
        __m128 t1 = _mm_shuffle_ps(col[0].mmvalue, col[1].mmvalue, _MM_SHUFFLE(2, 2, 2, 2));
        __m128 c0 = _mm_shuffle_ps(t0, col[2].mmvalue, _MM_SHUFFLE(0, 0, 2, 0));
        __m128 c1 = _mm_shuffle_ps(t0, col[2].mmvalue, _MM_SHUFFLE(0, 1, 3, 1));
        __m128 c2 = _mm_shuffle_ps(t1, col[2].mmvalue, _MM_SHUFFLE(0, 2, 2, 0));
        return FastMat3<T>(c0, c1, c2);
    }

    inline T norm2() const { return col[0].norm2() + col[1].norm2() + col[2].norm2(); }
    inline T norm() const { return std::sqrt(norm2()); }
    ////////////////////////////////////////////////////////////////////////////////
    inline bool operator==(const FastMat3<T>& B) const { return col[0] == B.col[0] && col[1] == B.col[1] && col[2] == B.col[2]; }
    inline bool operator!=(const FastMat3<T>& B) const { return !(*this == B); }
    inline FastMat3<T> operator-() const { return FastMat3<T>(-col[0], -col[1], -col[2]); }

    ////////////////////////////////////////////////////////////////////////////////
    // accessor
    template<class IndexType> inline const FastVec3<T>& operator[](const IndexType i) const { assert(i < 4); return col[i]; }
    template<class IndexType> inline FastVec3<T>& operator[](const IndexType i) { assert(i < 4); return col[i]; }
    ////////////////////////////////////////////////////////////////////////////////
    // member variable
    FastVec3<T> col[3];
};
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}     // end namespace NTCodeBase
