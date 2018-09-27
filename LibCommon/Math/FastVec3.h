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

#include <smmintrin.h>
#include <LibCommon/MathTypes.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __GNUC__
class __attribute__((aligned(16))) FastVec3
#else
_MM_ALIGN16 class FastVec3
#endif
{
public:
    /// Constructors
    inline FastVec3() : mmvalue(_mm_setzero_ps()) {}
    inline FastVec3(float x, float y, float z) : mmvalue(_mm_set_ps(0, z, y, x)) {}
    inline FastVec3(__m128 m) : mmvalue(m) {}

    /// Arithmetic operators with FastVec3
    inline FastVec3 operator+(const FastVec3& b) const { return _mm_add_ps(mmvalue, b.mmvalue); }
    inline FastVec3 operator-(const FastVec3& b) const { return _mm_sub_ps(mmvalue, b.mmvalue); }
    inline FastVec3 operator*(const FastVec3& b) const { return _mm_mul_ps(mmvalue, b.mmvalue); }
    inline FastVec3 operator/(const FastVec3& b) const { return _mm_div_ps(mmvalue, b.mmvalue); }

    /// Assignation and arithmetic operators with FastVec3
    inline FastVec3& operator+=(const FastVec3& b) { mmvalue = _mm_add_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3& operator-=(const FastVec3& b) { mmvalue = _mm_sub_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3& operator*=(const FastVec3& b) { mmvalue = _mm_mul_ps(mmvalue, b.mmvalue); return *this; }
    inline FastVec3& operator/=(const FastVec3& b) { mmvalue = _mm_div_ps(mmvalue, b.mmvalue); return *this; }

    /// Arithmetic operators with floats
    inline FastVec3 operator+(float b) const { return _mm_add_ps(mmvalue, _mm_set1_ps(b)); }
    inline FastVec3 operator-(float b) const { return _mm_sub_ps(mmvalue, _mm_set1_ps(b)); }
    inline FastVec3 operator*(float b) const { return _mm_mul_ps(mmvalue, _mm_set1_ps(b)); }
    inline FastVec3 operator/(float b) const { return _mm_div_ps(mmvalue, _mm_set1_ps(b)); }

    /// Assignation and arithmetic operators with float
    inline FastVec3& operator+=(float b) { mmvalue = _mm_add_ps(mmvalue, _mm_set1_ps(b)); return *this; }
    inline FastVec3& operator-=(float b) { mmvalue = _mm_sub_ps(mmvalue, _mm_set1_ps(b)); return *this; }
    inline FastVec3& operator*=(float b) { mmvalue = _mm_mul_ps(mmvalue, _mm_set1_ps(b)); return *this; }
    inline FastVec3& operator/=(float b) { mmvalue = _mm_div_ps(mmvalue, _mm_set1_ps(b)); return *this; }

    /// Equality operators
    inline bool operator==(const FastVec3& b) const
    {
        return (((_mm_movemask_ps(_mm_cmpeq_ps(mmvalue, b.mmvalue))) & 0x7) == 0x7);
    }

    inline bool operator!=(const FastVec3& b) const { return !(*this == b); }

    /// Unary minus operator
    inline FastVec3 operator-() const { return _mm_xor_ps(mmvalue, SIGN_MASK); }

    /// Subscript operator
    // Note: there is not bound checking here.
    inline const float& operator[](const int i) const
    {
        return i == 0 ? this->x : (i == 1 ? this->y : this->z);
    }

    inline float& operator[](const int i)
    {
        return i == 0 ? this->x : (i == 1 ? this->y : this->z);
    }

    /// Cross product
    inline FastVec3 cross(const FastVec3& b) const
    {
        return _mm_sub_ps(
            _mm_mul_ps(_mm_shuffle_ps(mmvalue, mmvalue, _MM_SHUFFLE(3, 0, 2, 1)), _mm_shuffle_ps(b.mmvalue, b.mmvalue, _MM_SHUFFLE(3, 1, 0, 2))),
            _mm_mul_ps(_mm_shuffle_ps(mmvalue, mmvalue, _MM_SHUFFLE(3, 1, 0, 2)), _mm_shuffle_ps(b.mmvalue, b.mmvalue, _MM_SHUFFLE(3, 0, 2, 1)))
            );
    }

    /// Dot product
    inline float dot(const FastVec3& b) const { return _mm_cvtss_f32(_mm_dp_ps(mmvalue, b.mmvalue, 0x71)); }
    /// Length of the vector
    inline float length() const { return _mm_cvtss_f32(_mm_sqrt_ss(_mm_dp_ps(mmvalue, mmvalue, 0x71))); }
    /// Returns the normalized vector
    inline FastVec3 normalize() const
    {
        // multiplying by rsqrt does not yield an accurate enough result, so we
        // divide by sqrt instead.
        return _mm_div_ps(mmvalue, _mm_sqrt_ps(_mm_dp_ps(mmvalue, mmvalue, 0xFF)));
    }

    ////////////////////////////////////////////////////////////////////////////////
    /// Direct access member variables.
    union
    {
        struct { float x, y, z; };
        __m128 mmvalue;
    };
private:
    // __m128 bits mask to target the floating point sign bit.
    static const __m128 SIGN_MASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));
};

/**
** External operators that maps to the actual FastVec3 method.
*/
inline FastVec3 operator+(float a, const FastVec3& b) { return b + a; }
inline FastVec3 operator-(float a, const FastVec3& b) { return FastVec3(_mm_set1_ps(a)) - b; }
inline FastVec3 operator*(float a, const FastVec3& b) { return b * a; }
inline FastVec3 operator/(float a, const FastVec3& b) { return FastVec3(_mm_set1_ps(a)) / b; }
