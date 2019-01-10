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

#include <atomic>
#include <limits>
#include <tbb/tbb.h>

#include <LibCommon/CommonSetup.h>
#undef min
#undef max

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::ParallelObjects {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class SpinLock {
public:
    SpinLock() = default;
    SpinLock(const SpinLock&) {}
    SpinLock& operator=(const SpinLock&) { return *this; }

    void lock() {
        while(m_Lock.test_and_set(std::memory_order_acquire)) {}
    }

    void unlock() {
        m_Lock.clear(std::memory_order_release);
    }

private:
    std::atomic_flag m_Lock = ATOMIC_FLAG_INIT;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class FastRandMT {
public:
    void seed(UInt seed_) { s_Seed = seed_; }

    /// Compute a pseudo-random integer
    /// Output value in range [0, 32767]
    Int rand() {
        m_Lock.lock();
        s_Seed = 214013u * s_Seed + 2531011u;
        m_Lock.unlock();
        return static_cast<Int>((s_Seed >> 16) & 0x7FFF);
    }

private:
    UInt     s_Seed = 0u;
    SpinLock m_Lock;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class DotProduct {
public:
    DotProduct(const StdVT<VecX<N, Real_t>>& vec1, const StdVT<VecX<N, Real_t>>& vec2) : m_Vec1(vec1), m_Vec2(vec2) {}
    DotProduct(DotProduct<N, Real_t>& pObj, tbb::split) : m_Vec1(pObj.m_Vec1), m_Vec2(pObj.m_Vec2) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += glm::dot(m_Vec1[i], m_Vec2[i]);
        }
    }

    void     join(DotProduct<N, Real_t>& pObj) { m_Result += pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t m_Result = 0;
    const StdVT<VecX<N, Real_t>>& m_Vec1;
    const StdVT<VecX<N, Real_t>>& m_Vec2;
};

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
class DotProduct<1, Real_t> {
public:
    DotProduct(const StdVT<Real_t>& vec1, const StdVT<Real_t>& vec2) : m_Vec1(vec1), m_Vec2(vec2) {}
    DotProduct(DotProduct<1, Real_t>& pObj, tbb::split) : m_Vec1(pObj.m_Vec1), m_Vec2(pObj.m_Vec2) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += m_Vec1[i] * m_Vec2[i];
        }
    }

    void     join(DotProduct<1, Real_t>& pObj) { m_Result += pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t               m_Result = 0;
    const StdVT<Real_t>& m_Vec1;
    const StdVT<Real_t>& m_Vec2;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MinElement {
public:
    MinElement(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {}
    MinElement(MinElement<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t vmin = glm::compMin(m_Vector[i]);
            m_Result = m_Result < vmin ? m_Result : vmin;
        }
    }

    void     join(MinElement<N, Real_t>& pObj) { m_Result = m_Result < pObj.m_Result ? m_Result : pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t m_Result = std::numeric_limits<Real_t>::max();
    const StdVT<VecX<N, Real_t>>& m_Vector;
};

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
class MinElement<1, Real_t> {
public:
    MinElement(const StdVT<Real_t>& vec) : m_Vector(vec) {}
    MinElement(MinElement<1, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result = m_Result < m_Vector[i] ? m_Result : m_Vector[i];
        }
    }

    void     join(MinElement<1, Real_t>& pObj) { m_Result = m_Result < pObj.m_Result ? m_Result : pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t               m_Result = std::numeric_limits<Real_t>::max();
    const StdVT<Real_t>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MaxElement {
public:
    MaxElement(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {}
    MaxElement(MaxElement<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t vmax = glm::compMax(m_Vector[i]);
            m_Result = m_Result > vmax ? m_Result : vmax;
        }
    }

    void     join(MaxElement<N, Real_t>& pObj) { m_Result = m_Result > pObj.m_Result ? m_Result : pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t m_Result = -std::numeric_limits<Real_t>::max();
    const StdVT<VecX<N, Real_t>>& m_Vector;
};

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
class MaxElement<1, Real_t> {
public:
    MaxElement(const StdVT<Real_t>& vec) : m_Vector(vec) {}
    MaxElement(MaxElement<1, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result = m_Result > m_Vector[i] ? m_Result : m_Vector[i];
        }
    }

    void     join(MaxElement<1, Real_t>& pObj) { m_Result = m_Result > pObj.m_Result ? m_Result : pObj.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t               m_Result = -std::numeric_limits<Real_t>::max();
    const StdVT<Real_t>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MaxAbs {
public:
    MaxAbs(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {}
    MaxAbs(MaxAbs<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t tmp = glm::compMax(glm::abs(m_Vector[i]));
            m_Result = m_Result > tmp ? m_Result : tmp;
        }
    }

    void     join(MaxAbs<N, Real_t>& vma) { m_Result = m_Result > vma.m_Result ? m_Result : vma.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t m_Result = 0;
    const StdVT<VecX<N, Real_t>>& m_Vector;
};

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
class MaxAbs<1, Real_t> {
public:
    MaxAbs(const StdVT<Real_t>& vec) : m_Vector(vec) {}
    MaxAbs(MaxAbs<1, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t tmp = std::abs(m_Vector[i]);
            m_Result = m_Result > tmp ? m_Result : tmp;
        }
    }

    void     join(MaxAbs<1, Real_t>& vma) { m_Result = m_Result > vma.m_Result ? m_Result : vma.m_Result; }
    Real_t getResult() const noexcept { return m_Result; }

private:
    Real_t               m_Result = 0;
    const StdVT<Real_t>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MaxNorm2 {
public:
    MaxNorm2(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {}
    MaxNorm2(MaxNorm2<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t mag2 = glm::length2(m_Vector[i]);
            m_Result = m_Result > mag2 ? m_Result : mag2;
        }
    }

    void     join(MaxNorm2<N, Real_t>& pObj) { m_Result = m_Result > pObj.m_Result ? m_Result : pObj.m_Result; }
    Real_t getResult() const noexcept { return std::sqrt(m_Result); }

private:
    Real_t m_Result = 0;
    const StdVT<VecX<N, Real_t>>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MinMaxElements {
public:
    MinMaxElements(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {
        if(vec.size() > 0) { m_ResultMax = vec[0]; }
    }

    MinMaxElements(MinMaxElements<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            const auto& vec = m_Vector[i];
            for(int j = 0; j < N; ++j) {
                m_ResultMin[j] = (m_ResultMin[j] < vec[j]) ? m_ResultMin[j] : vec[j];
                m_ResultMax[j] = (m_ResultMax[j] > vec[j]) ? m_ResultMax[j] : vec[j];
            }
        }
    }

    void join(MinMaxElements<N, Real_t>& pObj) {
        for(int j = 0; j < N; ++j) {
            m_ResultMin[j] = (m_ResultMin[j] < pObj.m_ResultMin[j]) ? m_ResultMin[j] : pObj.m_ResultMin[j];
            m_ResultMax[j] = (m_ResultMax[j] > pObj.m_ResultMax[j]) ? m_ResultMax[j] : pObj.m_ResultMax[j];
        }
    }

    const VecX<N, Real_t>& getMin() const noexcept { return m_ResultMin; }
    const VecX<N, Real_t>& getMax() const noexcept { return m_ResultMax; }

private:
    VecX<N, Real_t>               m_ResultMin = VecX<N, Real_t>(std::numeric_limits<Real_t>::max());
    VecX<N, Real_t>               m_ResultMax = VecX<N, Real_t>(-std::numeric_limits<Real_t>::max());
    const StdVT<VecX<N, Real_t>>& m_Vector;
};

////////////////////////////////////////////////////////////////////////////////
// N = 0 : plain C array
template<class Real_t>
class MinMaxElements<0, Real_t> {
public:
    MinMaxElements(const Real_t* vec) : m_Vector(vec) {
        m_ResultMax = vec[0];
    }

    MinMaxElements(MinMaxElements<0, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_ResultMin = m_ResultMin < m_Vector[i] ? m_ResultMin : m_Vector[i];
            m_ResultMax = m_ResultMax > m_Vector[i] ? m_ResultMax : m_Vector[i];
        }
    }

    void join(MinMaxElements<0, Real_t>& pObj) {
        m_ResultMin = m_ResultMin < pObj.m_ResultMin ? m_ResultMin : pObj.m_ResultMin;
        m_ResultMax = m_ResultMax > pObj.m_ResultMax ? m_ResultMax : pObj.m_ResultMax;
    }

    Real_t getMin() const noexcept { return m_ResultMin; }
    Real_t getMax() const noexcept { return m_ResultMax; }

private:
    Real_t m_ResultMin = std::numeric_limits<Real_t>::max();
    Real_t m_ResultMax = -std::numeric_limits<Real_t>::max();

    const Real_t* m_Vector;
};

////////////////////////////////////////////////////////////////////////////////
template<class Real_t>
class MinMaxElements<1, Real_t> {
public:
    MinMaxElements(const StdVT<Real_t>& vec) : m_Vector(vec) {
        if(vec.size() > 0) { m_ResultMax = vec[0]; }
    }

    MinMaxElements(MinMaxElements<1, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_ResultMin = m_ResultMin < m_Vector[i] ? m_ResultMin : m_Vector[i];
            m_ResultMax = m_ResultMax > m_Vector[i] ? m_ResultMax : m_Vector[i];
        }
    }

    void join(MinMaxElements<1, Real_t>& pObj) {
        m_ResultMin = m_ResultMin < pObj.m_ResultMin ? m_ResultMin : pObj.m_ResultMin;
        m_ResultMax = m_ResultMax > pObj.m_ResultMax ? m_ResultMax : pObj.m_ResultMax;
    }

    Real_t getMin() const noexcept { return m_ResultMin; }
    Real_t getMax() const noexcept { return m_ResultMax; }

private:
    Real_t m_ResultMin = std::numeric_limits<Real_t>::max();
    Real_t m_ResultMax = -std::numeric_limits<Real_t>::max();

    const StdVT<Real_t>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class Real_t>
class MinMaxNorm2 {
public:
    MinMaxNorm2(const StdVT<VecX<N, Real_t>>& vec) : m_Vector(vec) {}
    MinMaxNorm2(MinMaxNorm2<N, Real_t>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            Real_t mag2 = glm::length2(m_Vector[i]);
            m_ResultMin = m_ResultMin < mag2 ? m_ResultMin : mag2;
            m_ResultMax = m_ResultMax > mag2 ? m_ResultMax : mag2;
        }
    }

    void join(MinMaxNorm2<N, Real_t>& pObj) {
        m_ResultMin = m_ResultMin < pObj.m_ResultMin ? m_ResultMin : pObj.m_ResultMin;
        m_ResultMax = m_ResultMax > pObj.m_ResultMax ? m_ResultMax : pObj.m_ResultMax;
    }

    Real_t getMin() const noexcept { return std::sqrt(m_ResultMin); }
    Real_t getMax() const noexcept { return std::sqrt(m_ResultMax); }

private:
    Real_t m_ResultMin = std::numeric_limits<Real_t>::max();
    Real_t m_ResultMax = 0;

    const StdVT<VecX<N, Real_t>>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
class VectorSum {
public:
    VectorSum(const StdVT<VecX<N, T>>& vec) : m_Vector(vec) {}
    VectorSum(VectorSum<N, T>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += m_Vector[i];
        }
    }

    void              join(VectorSum<N, T>& pObj) { m_Result += pObj.m_Result; }
    const VecX<N, T>& getResult() const noexcept { return m_Result; }
private:
    VecX<N, T>               m_Result = VecX<N, T>(0);
    const StdVT<VecX<N, T>>& m_Vector;
};

template<class T>
class VectorSum<1, T> {
public:
    VectorSum(const StdVT<T>& vec) : m_Vector(vec) {}
    VectorSum(VectorSum<1, T>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += m_Vector[i];
        }
    }

    void join(VectorSum<1, T>& pObj) { m_Result += pObj.m_Result; }
    T    getResult() const noexcept { return m_Result; }

private:
    T               m_Result = T(0);
    const StdVT<T>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<Int N, class T>
class VectorSumSqr {
public:
    VectorSumSqr(const StdVT<VecX<N, T>>& vec) : m_Vector(vec) {}
    VectorSumSqr(VectorSumSqr<N, T>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += glm::length2(m_Vector[i]);
        }
    }

    void join(VectorSumSqr<N, T>& pObj) { m_Result += pObj.m_Result; }
    T    getResult() const noexcept { return m_Result; }
private:
    T m_Result = T(0);
    const StdVT<VecX<N, T>>& m_Vector;
};

template<class T>
class VectorSumSqr<1, T> {
public:
    VectorSumSqr(const StdVT<T>& vec) : m_Vector(vec) {}
    VectorSumSqr(VectorSumSqr<1, T>& pObj, tbb::split) : m_Vector(pObj.m_Vector) {}

    void operator()(const tbb::blocked_range<size_t>& r) {
        for(size_t i = r.begin(); i != r.end(); ++i) {
            m_Result += m_Vector[i] * m_Vector[i];
        }
    }

    void join(VectorSumSqr<1, T>& pObj) { m_Result += pObj.m_Result; }
    T    getResult() const noexcept { return m_Result; }

private:
    T               m_Result = T(0);
    const StdVT<T>& m_Vector;
};

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::ParallelObjects
