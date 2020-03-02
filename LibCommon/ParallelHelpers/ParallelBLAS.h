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

#include <cassert>
#include <vector>

#if !defined(Q_MOC_RUN)
#include <tbb/tbb.h>
#endif

#include <LibCommon/ParallelHelpers/ParallelObjects.h>
#include <LibCommon/ParallelHelpers/ParallelExec.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase::ParallelBLAS {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// dot products
//
template<class Real_t>
inline Real_t dotProduct(const StdVT<Real_t>& x, const StdVT<Real_t>& y) {
    NT_REQUIRE(x.size() == y.size());
    ParallelObjects::DotProduct<1, Real_t> pObj(x, y);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

template<Int N, class Real_t>
inline Real_t dotProduct(const StdVT<VecX<N, Real_t>>& x, const StdVT<VecX<N, Real_t>>& y) {
    NT_REQUIRE(x.size() == y.size());
    ParallelObjects::DotProduct<N, Real_t> pObj(x, y);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Real_t>
inline Real_t norm2(const StdVT<Real_t>& x) {
    ParallelObjects::VectorSumSqr<1, Real_t> pObj(x);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

template<Int N, class Real_t>
inline Real_t norm2(const StdVT<VecX<N, Real_t>>& x) {
    ParallelObjects::VectorSumSqr<N, Real_t> pObj(x);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class VectorType>
inline StdVT<VectorType> add(const StdVT<VectorType>& x, const StdVT<VectorType>& y) {
    NT_REQUIRE(x.size() == y.size());
    StdVT<VectorType> z(x.size());
    ParallelExec::run(z.size(), [&](size_t i) { z[i] = x[i] + y[i]; });
    return z;
}

template<class VectorType>
inline StdVT<VectorType> minus(const StdVT<VectorType>& x, const StdVT<VectorType>& y) {
    NT_REQUIRE(x.size() == y.size());
    StdVT<VectorType> z(x.size());
    ParallelExec::run(z.size(), [&](size_t i) { z[i] = x[i] - y[i]; });
    return z;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// saxpy (y=alpha*x+y)
//
template<class Real_t, class VectorType>
inline void addScaled(Real_t alpha, const StdVT<VectorType>& x, StdVT<VectorType>& y) {
    ParallelExec::run(x.size(), [&, alpha](size_t i) { y[i] += alpha * x[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// y=x+beta*y
//
template<class Real_t, class VectorType>
inline void scaledAdd(Real_t beta, const StdVT<VectorType>& x, StdVT<VectorType>& y) {
    ParallelExec::run(x.size(), [&, beta](size_t i) { y[i] = beta * y[i] + x[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// x *= alpha
//
template<class Real_t, class VectorType>
inline void scale(Real_t alpha, StdVT<VectorType>& x) {
    ParallelExec::run(x.size(), [&, alpha](size_t i) { x[i] *= alpha; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// y = x * alpha
template<class Real_t, class VectorType>
inline StdVT<VectorType> multiply(Real_t alpha, const StdVT<VectorType>& x) {
    StdVT<VectorType> y(x.size());
    ParallelExec::run(x.size(), [&, alpha](size_t i) { y[i] = x[i] * alpha; });
    return y;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase::ParallelBLAS
