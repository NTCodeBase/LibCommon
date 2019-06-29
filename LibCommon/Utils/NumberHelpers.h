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

#include <cstdint>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <locale>
#include <random>

#include <LibCommon/CommonSetup.h>
#include <LibCommon/ParallelHelpers/ParallelExec.h>
#include <LibCommon/ParallelHelpers/ParallelObjects.h>
#include <LibCommon/Math/MathHelpers.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::NumberHelpers {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
constexpr bool isFloat() { return std::is_same_v<T, float>; }

template<class T>
String nameRealT() { return isFloat<T>() ? String("float") : String("double"); }
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
bool isValidNumber(T x) {
    return !std::isnan(x) && !std::isinf(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
T isPowerOfTwo(T n) {
    return !(n == 0) && !(n & (n - 1));
}

template<class T>
T nextPowerOfTwo(T n) {
    if(n > 0) {
        --n;

        n |= n >> 1;
        n |= n >> 2;
        n |= n >> 4;
        n |= n >> 8;
        n |= n >> 16;

        return n + static_cast<T>(1);
    }
    return 0;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IndexType, Int N, class Real_t>
auto createGrid(const VecX<N, Real_t>& bmin, const VecX<N, Real_t>& bmax, const VecX<N, Real_t>& spacing) {
    VecX<N, Real_t>    fgrid = (bmax - bmin) / spacing;
    VecX<N, IndexType> result;
    for(Int d = 0; d < N; ++d) {
        result[d] = static_cast<IndexType>(ceil(fgrid[d]));
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IndexType, Int N, class Real_t>
auto createGrid(const VecX<N, Real_t>& bmin, const VecX<N, Real_t>& bmax, Real_t spacing) {
    VecX<N, Real_t>    fgrid = (bmax - bmin) / spacing;
    VecX<N, IndexType> result;
    for(Int d = 0; d < N; ++d) {
        result[d] = static_cast<IndexType>(ceil(fgrid[d]));
    }
    return result;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
class MT_iRandom {
    using Distribution = std::uniform_int_distribution<T>;
    constexpr static UInt s_CacheSize = 1048576u;
public:
    MT_iRandom(T start = T(0), T end = std::numeric_limits<T>::max()) : m_Dist(Distribution(start, end)) { generateCache(); }

    T rnd() {
        m_Lock.lock();
        if(m_CacheIdx >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        T tmp = m_Cache[m_CacheIdx++];
        m_Lock.unlock();
        return tmp;
    }

    template<class VT>
    VT vrnd() {
        static_assert(sizeof(typename VT::value_type) == sizeof(T));
        UInt N = static_cast<UInt>(VT::length());
        m_Lock.lock();
        if(m_CacheIdx + N >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        auto oldIdx = m_CacheIdx;
        m_CacheIdx += N;
        m_Lock.unlock();
        VT result;
        std::memcpy(&result, &m_Cache[oldIdx], sizeof(VT));
        return result;
    }

    template<class Matrix>
    Matrix mrnd() {
        static_assert(sizeof(typename Matrix::value_type) == sizeof(T));
        UInt MxN = static_cast<UInt>(Matrix::length() * Matrix::col_type::length());
        m_Lock.lock();
        if(m_CacheIdx + MxN >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        auto oldIdx = m_CacheIdx;
        m_CacheIdx += MxN;
        m_Lock.unlock();
        Matrix result;
        std::memcpy(&result, &m_Cache[oldIdx], sizeof(Matrix));
        return result;
    }

    void setSeed(UInt seed) {
        m_Lock.lock();
        m_CacheIdx  = 0;
        m_Generator = std::mt19937(seed);
        m_Lock.unlock();
    }

private:
    void generateCache() {
        for(UInt i = 0; i < s_CacheSize; ++i) {
            m_Cache[i] = m_Dist(m_Generator);
        }
    }

    ParallelObjects::SpinLock m_Lock {};
    std::mt19937              m_Generator { (std::random_device())() };
    Distribution              m_Dist;

    T    m_Cache[s_CacheSize];
    UInt m_CacheIdx = 0;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
class MT_fRandom {
    using Distribution = std::uniform_real_distribution<T>;
    constexpr static UInt s_CacheSize = 1048576u;
public:
    MT_fRandom(T start = T(0), T end = std::numeric_limits<T>::max()) : m_Dist(Distribution(start, end)) { generateCache(); }

    T rnd() {
        m_Lock.lock();
        if(m_CacheIdx >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        T tmp = m_Cache[m_CacheIdx++];
        m_Lock.unlock();
        return tmp;
    }

    template<class VT>
    VT vrnd() {
        static_assert(sizeof(typename VT::value_type) == sizeof(T));
        UInt N = static_cast<UInt>(VT::length());
        m_Lock.lock();
        if(m_CacheIdx + N >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        auto oldIdx = m_CacheIdx;
        m_CacheIdx += N;
        m_Lock.unlock();
        VT result;
        std::memcpy(&result, &m_Cache[oldIdx], sizeof(VT));
        return result;
    }

    template<class Matrix>
    Matrix mrnd() {
        static_assert(sizeof(typename Matrix::value_type) == sizeof(T));
        UInt MxN = static_cast<UInt>(Matrix::length() * Matrix::col_type::length());
        m_Lock.lock();
        if(m_CacheIdx + MxN >= s_CacheSize) {
            m_CacheIdx = 0;
        }
        auto oldIdx = m_CacheIdx;
        m_CacheIdx += MxN;
        m_Lock.unlock();
        Matrix result;
        std::memcpy(&result, &m_Cache[oldIdx], sizeof(Matrix));
        return result;
    }

    void setSeed(UInt seed) {
        m_Lock.lock();
        m_CacheIdx  = 0;
        m_Generator = std::mt19937(seed);
        m_Lock.unlock();
    }

private:
    void generateCache() {
        for(UInt i = 0; i < s_CacheSize; ++i) {
            m_Cache[i] = m_Dist(m_Generator);
        }
    }

    ParallelObjects::SpinLock m_Lock {};
    std::mt19937              m_Generator { (std::random_device())() };
    Distribution              m_Dist;

    T    m_Cache[s_CacheSize];
    UInt m_CacheIdx = 0;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
class iRand {
public:
    static auto rnd() { return s_Rand.rnd(); }
    template<class VT> static auto vrnd() { return s_Rand.template vrnd<VT>(); }
    template<class Matrix> static auto mrnd() { return s_Rand.template mrnd<Matrix>(); }
    static void setSeed(UInt seed) { s_Rand.setSeed(seed); }
private:
    static inline MT_iRandom<T> s_Rand = MT_iRandom<T>();
};

template<class T>
class fRand {
public:
    static auto rnd() { return s_Rand.rnd(); }
    template<class VT> static auto vrnd() { return s_Rand.template vrnd<VT>(); }
    template<class Matrix> static auto mrnd() { return s_Rand.template mrnd<Matrix>(); }
    static void setSeed(UInt seed) { s_Rand.setSeed(seed); }
private:
    static inline MT_fRandom<T> s_Rand = MT_fRandom<T>();
};

template<class T>
class fRand01 {
public:
    static auto rnd() { return s_Rand.rnd(); }
    template<class VT> static auto vrnd() { return s_Rand.template vrnd<VT>(); }
    template<class Matrix> static auto mrnd() { return s_Rand.template mrnd<Matrix>(); }
    static void setSeed(UInt seed) { s_Rand.setSeed(seed); }
private:
    static inline MT_fRandom<T> s_Rand = MT_fRandom<T>(T(0), T(1));
};

template<class T>
class fRand11 {
public:
    static auto rnd() { return s_Rand.rnd(); }
    template<class VT> static auto vrnd() { return s_Rand.template vrnd<VT>(); }
    template<class Matrix> static auto mrnd() { return s_Rand.template mrnd<Matrix>(); }
    static void setSeed(UInt seed) { s_Rand.setSeed(seed); }
private:
    static inline MT_fRandom<T> s_Rand = MT_fRandom<T>(T(-1), T(1));
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IntType, class SizeType>
StdVT<IntType> generateRandomIntVector(SizeType size, IntType start = 0, IntType end = std::numeric_limits<IntType>::max()) {
    StdVT<IntType> v(size);
    for(SizeType i = 0; i < size; ++i) {
        v[i] = MT_iRandom<IntType>::rnd();
    }
    return v;
}

template<class Real_t, class SizeType>
StdVT<Real_t> generateRandomRealVector(SizeType size, Real_t start = Real_t(0), Real_t end = std::numeric_limits<Real_t>::max()) {
    StdVT<Real_t> v(size);
    for(SizeType i = 0; i < size; ++i) {
        v[i] = MT_fRandom<Real_t>::rnd();
    }
    return v;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Transforms even the sequence 0,1,2,3,... into reasonably good random numbers
// Challenge: improve on this in speed and "randomness"!
// This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
inline unsigned int randhash(unsigned int seed) {
    unsigned int i = (seed ^ 0xA3C59AC3u) * 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    return i;
}

// the inverse of randhash
inline unsigned int unhash(unsigned int h) {
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= 0xA3C59AC3u;
    return h;
}

// returns repeatable stateless pseudo-random number in [0,1]
template<class T>
inline T frandhash(unsigned int seed) {
    return T(randhash(seed)) / static_cast<T>(UINT_MAX);
}

// returns repeatable stateless pseudo-random number in [a,b]
template<class T>
inline T frandhash(T a, T b, unsigned int seed) {
    return (b - a) * (randhash(seed) / static_cast<T>(UINT_MAX)) + a;
}

template<class T>
inline T frandhash11(unsigned int seed) {
    return frandhash(T(-1.0), T(1.0), seed);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class FastRand {
public:
    void seed(UInt seed_) { s_Seed = seed_; }

    /// Compute a pseudo-random integer
    /// Output value in range [0, 32767]
    Int rand() {
        s_Seed = (214013u * s_Seed + 2531011u);
        return static_cast<Int>((s_Seed >> 16) & 0x7FFF);
    }

private:
    UInt s_Seed = 0u;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
bool isInside(const VecX<N, Real_t>& ppos, const VecX<N, Real_t>& bMin, const VecX<N, Real_t>& bMax) {
    for(Int d = 0; d < N; ++d) {
        if(ppos[d] < bMin[d] || ppos[d] > bMax[d]) {
            return false;
        }
    }
    return true;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class RealType1, class RealType2>
void jitter(VecX<N, RealType1>& ppos, RealType2 maxJitter) {
    for(Int j = 0; j < N; ++j) {
        ppos += fRand11<RealType1>::rnd() * static_cast<RealType1>(maxJitter);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void translate(StdVT_VecX<N, Real_t>& points, const VecX<N, Real_t>& translation) {
    ParallelExec::run(points.size(), [&](size_t i) {
                          points[i] = points[i] + translation;
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void scale(StdVT_VecX<N, Real_t>& points, const VecX<N, Real_t>& scale) {
    ParallelExec::run(points.size(), [&](size_t i) {
                          points[i] = points[i] * scale;
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void rotate(StdVT_VecX<N, Real_t>& points, const VecX<N, Real_t>& rotation) {
    Real_t azimuth = rotation[1];
    Real_t elevation = rotation[0];
    Real_t roll = rotation[2];
    Real_t sinA, cosA, sinE, cosE, sinR, cosR;

    Vec3<Real_t> R[3];

    sinA = static_cast<Real_t>(sin(azimuth));
    cosA = static_cast<Real_t>(cos(azimuth));
    sinE = static_cast<Real_t>(sin(elevation));
    cosE = static_cast<Real_t>(cos(elevation));
    sinR = static_cast<Real_t>(sin(roll));
    cosR = static_cast<Real_t>(cos(roll));

    R[0][0] = cosR * cosA - sinR * sinA * sinE;
    R[0][1] = sinR * cosA + cosR * sinA * sinE;
    R[0][2] = -sinA * cosE;

    R[1][0] = -sinR * cosE;
    R[1][1] = cosR * cosE;
    R[1][2] = sinE;

    R[2][0] = cosR * sinA + sinR * cosA * sinE;
    R[2][1] = sinR * sinA - cosR * cosA * sinE;
    R[2][2] = cosA * cosE;

    ParallelExec::run(points.size(), [&](size_t i) {
                          const auto& pi = points[i];
                          points[i]      = Vec3<Real_t>(glm::dot(R[0], pi),
                                                        glm::dot(R[1], pi),
                                                        glm::dot(R[2], pi));
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
void transform(StdVT_VecX<N, Real_t>& points, const VecX<N, Real_t>& translation, const VecX<N, Real_t>& scale) {
    ParallelExec::run(points.size(), [&](size_t i) {
                          points[i] = points[i] * scale + translation;
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
void transform(StdVT_Vec3<Real_t>& points, const Vec3<Real_t>& translation, const Vec3<Real_t>& scale, const Vec3<Real_t>& rotation) {
    Real_t azimuth = rotation[0];
    Real_t elevation = rotation[1];
    Real_t roll = rotation[2];
    Real_t sinA, cosA, sinE, cosE, sinR, cosR;

    Vec3<Real_t> R[3];

    sinA = static_cast<Real_t>(sin(azimuth));
    cosA = static_cast<Real_t>(cos(azimuth));
    sinE = static_cast<Real_t>(sin(elevation));
    cosE = static_cast<Real_t>(cos(elevation));
    sinR = static_cast<Real_t>(sin(roll));
    cosR = static_cast<Real_t>(cos(roll));

    R[0][0] = cosR * cosA - sinR * sinA * sinE;
    R[0][1] = sinR * cosA + cosR * sinA * sinE;
    R[0][2] = -sinA * cosE;

    R[1][0] = -sinR * cosE;
    R[1][1] = cosR * cosE;
    R[1][2] = sinE;

    R[2][0] = cosR * sinA + sinR * cosA * sinE;
    R[2][1] = sinR * sinA - cosR * cosA * sinE;
    R[2][2] = cosA * cosE;

    ParallelExec::run(points.size(), [&](size_t i) {
                          const auto& pi = points[i];
                          Vec3<Real_t> tmp(glm::dot(R[0], pi),
                                           glm::dot(R[1], pi),
                                           glm::dot(R[2], pi));
                          points[i] = tmp * scale + translation;
                      });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class IndexType, class Real_t>
auto generatePointsOnCircle(IndexType nPoints, Real_t sphereRadius = Real_t(1)) {
    StdVT_Vec2<Real_t> points;
    points.reserve(nPoints);
    ////////////////////////////////////////////////////////////////////////////////
    const auto b   = static_cast<IndexType>(std::round(2 /*=alpha*/ * std::sqrt(nPoints)));
    const auto phi = Real_t((std::sqrt(5) + 1.0) * 0.5);
    for(IndexType k = 1; k <= nPoints; ++k) {
        auto r = (k > nPoints - b) ? sphereRadius : sphereRadius* std::sqrt(k - 0.5) / std::sqrt(nPoints - 0.5 * b - 0.5);
        auto theta = Real_t(2.0 * M_PI * k / phi / phi);
        points.push_back(Vec2<Real_t>(r * std::cos(theta), r * std::sin(theta)));
    }
    return points;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class IndexType, class Real_t>
auto generatePointsOnSphere(IndexType nPoints, Real_t radius = Real_t(1), Real_t spanPolarAngle = Real_t(M_PI)) {
    if constexpr(N == 2) {
        NT_UNUSED(spanPolarAngle);
        StdVT_Vec2<Real_t> points;
        points.reserve(nPoints);
        ////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////
        const auto dphi = Real_t(2 * M_PI / nPoints);
        for(IndexType k = 0; k < nPoints; ++k) {
            auto phi = k * dphi;
            points.push_back(Vec2<Real_t>(radius * std::cos(phi), radius * std::sin(phi)));
        }
        return points;
    } else {
        StdVT_Vec3<Real_t> points;
        points.reserve(nPoints);
        auto a    = Real_t(4.0 * M_PI / nPoints);
        auto d    = std::sqrt(a);
        auto Mv   = static_cast<IndexType>(std::round(M_PI / d));
        auto dv   = Real_t(M_PI / Mv);
        auto dphi = a / dv;
        for(IndexType m = 0; m < Mv; ++m) {
            auto v      = Real_t(spanPolarAngle * (m + 0.5) / Mv);
            auto Mphi   = static_cast<IndexType>(std::round(2.0 * M_PI * std::sin(v) / dphi));
            auto jitter = M_PI * fRand01<Real_t>::rnd();
            for(IndexType n = 0; n < Mphi; ++n) {
                auto phi = Real_t(2 * M_PI * n / Mphi) + jitter;
                points.push_back(Vec3<Real_t>(radius * std::sin(v) * std::cos(phi),
                                              radius * std::cos(v),
                                              radius * std::sin(v) * std::sin(phi)));
            }
        }
        return points;
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::NumberHelpers
