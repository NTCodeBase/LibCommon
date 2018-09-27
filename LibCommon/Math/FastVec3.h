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

#include <cassert>
#include <smmintrin.h>
#include <LibCommon/MathTypes.h>
#include <LibCommon/CommonMacros.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
class __NT_ALIGN16 FastVec3
{
public:
    inline FastVec3() { static_assert(sizeof(T) == sizeof(float) && alignof(FastVec3<T>) == 16, "Error: Size or alignment is not correct!"); }
    inline FastVec3(__m128 mval) : mmvalue(mval) {}
    inline FastVec3(T x) : mmvalue(_mm_set1_ps(x)) {}
    inline FastVec3(T x, T y, T z) : mmvalue(_mm_setr_ps(x, y, z, 1)) {}
    inline FastVec3(const Vec3<T>& v) : mmvalue(_mm_setr_ps(v.x, v.y, v.z, 1)) {}
    inline FastVec3(const FastVec3<T>& other) : mmvalue(other.mmvalue) {}
    ////////////////////////////////////////////////////////////////////////////////
    inline FastVec3<T>& operator=(const FastVec3<T>& other) { mmvalue = other.mmvalue; return *this; }
    inline FastVec3<T>& operator=(const Vec3<T>& v) { mmvalue = _mm_setr_ps(v.x, v.y, v.z, 0); return *this; }
    inline operator Vec3<T>() const { return v3; }
    ////////////////////////////////////////////////////////////////////////////////
    inline FastVec3<T> operator+(const FastVec3<T>& b) const { return _mm_add_ps(mmvalue, b.mmvalue); }
    inline FastVec3<T> operator-(const FastVec3<T>& b) const { return _mm_sub_ps(mmvalue, b.mmvalue); }
    inline FastVec3<T> operator*(const FastVec3<T>& b) const { return _mm_mul_ps(mmvalue, b.mmvalue); }
    inline FastVec3<T> operator/(const FastVec3<T>& b) const { return _mm_div_ps(mmvalue, b.mmvalue); }

    inline FastVec3<T>& operator+=(const FastVec3<T>& b) { mmvalue = _mm_add_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3<T>& operator-=(const FastVec3<T>& b) { mmvalue = _mm_sub_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3<T>& operator*=(const FastVec3<T>& b) { mmvalue = _mm_mul_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3<T>& operator/=(const FastVec3<T>& b) { mmvalue = _mm_div_ps(mmvalue, b.mmvalue); return *this; }

    inline FastVec3<T> operator+(T a) const { return _mm_add_ps(mmvalue, _mm_set1_ps(a)); }
    inline FastVec3<T> operator-(T a) const { return _mm_sub_ps(mmvalue, _mm_set1_ps(a)); }
    inline FastVec3<T> operator*(T a) const { return _mm_mul_ps(mmvalue, _mm_set1_ps(a)); }
    inline FastVec3<T> operator/(T a) const { assert(a != 0); return _mm_div_ps(mmvalue, _mm_set1_ps(a)); }

    inline FastVec3<T>& operator+=(T a) { mmvalue = _mm_add_ps(mmvalue, _mm_set1_ps(a)); return *this; }
    inline FastVec3<T>& operator-=(T a) { mmvalue = _mm_sub_ps(mmvalue, _mm_set1_ps(a)); return *this; }
    inline FastVec3<T>& operator*=(T a) { mmvalue = _mm_mul_ps(mmvalue, _mm_set1_ps(a)); return *this; }
    inline FastVec3<T>& operator/=(T a) { assert(a != 0); mmvalue = _mm_div_ps(mmvalue, _mm_set1_ps(a)); return *this; }
    ////////////////////////////////////////////////////////////////////////////////
    // geometric functions
    inline FastVec3<T> cross(const FastVec3<T>& b) const
    {
        return _mm_sub_ps(
            _mm_mul_ps(_mm_shuffle_ps(mmvalue, mmvalue, _MM_SHUFFLE(3, 0, 2, 1)),
                       _mm_shuffle_ps(b.mmvalue, b.mmvalue, _MM_SHUFFLE(3, 1, 0, 2))),
            _mm_mul_ps(_mm_shuffle_ps(mmvalue, mmvalue, _MM_SHUFFLE(3, 1, 0, 2)),
                       _mm_shuffle_ps(b.mmvalue, b.mmvalue, _MM_SHUFFLE(3, 0, 2, 1)))
            );
    }

    inline FastVec3<T> normalized() const { return _mm_div_ps(mmvalue, _mm_sqrt_ps(_mm_dp_ps(mmvalue, mmvalue, 0x7F))); }
    inline T dot(const FastVec3<T>& b) const { return _mm_cvtss_f32(_mm_dp_ps(mmvalue, b.mmvalue, 0x71)); }
    inline T norm2() const { return _mm_cvtss_f32(_mm_dp_ps(mmvalue, mmvalue, 0x71)); }
    inline T norm() const { return std::sqrt(norm2()); }
    ////////////////////////////////////////////////////////////////////////////////
    inline bool operator==(const FastVec3<T>& b) const { return (((_mm_movemask_ps(_mm_cmpeq_ps(mmvalue, b.mmvalue))) & 0x7) == 0x7); }
    inline bool operator!=(const FastVec3<T>& b) const { return !(*this == b); }
    inline FastVec3<T> operator-() const { return _mm_xor_ps(mmvalue, s_SignMask); }
    ////////////////////////////////////////////////////////////////////////////////
    // accessor
    template<class IndexType> inline const T& operator[](const IndexType i) const { assert(i < 4); return val[i]; }
    template<class IndexType> inline T& operator[](const IndexType i) { assert(i < 4); return val[i]; }
    ////////////////////////////////////////////////////////////////////////////////
    // member variable
    union
    {
        struct { Vec3<T> v3; T dummy; };
        struct { T x, y, z, w; };
        struct { T val[4]; };
        __m128 mmvalue;
    };
private:
    // __m128 bits mask to target the Ting point sign bit.
    inline static const __m128 s_SignMask = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T> inline FastVec3<T> operator+(T a, const FastVec3<T>& b) { return b + a; }
template<class T> inline FastVec3<T> operator-(T a, const FastVec3<T>& b) { return FastVec3(_mm_set1_ps(a)) - b; }
template<class T> inline FastVec3<T> operator*(T a, const FastVec3<T>& b) { return b * a; }
template<class T> inline FastVec3<T> operator/(T a, const FastVec3<T>& b) { return FastVec3(_mm_set1_ps(a)) / b; }
