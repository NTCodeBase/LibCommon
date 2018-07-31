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
#include <vector>

#include <tbb/tbb.h>
#include <ParallelHelpers/ParallelObjects.h>
#include <ParallelHelpers/Scheduler.h>


//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace ParallelBLAS
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// dot products
//
template<class RealType>
inline RealType dotProduct(const Vector<RealType>& x, const Vector<RealType>& y)
{
    __NT_REQUIRE(x.size() == y.size());
    ParallelObjects::DotProduct<1, RealType> pObj(x, y);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

template<Int N, class RealType>
inline RealType dotProduct(const Vector<VecX<N, RealType> >& x, const Vector<VecX<N, RealType> >& y)
{
    __NT_REQUIRE(x.size() == y.size());
    ParallelObjects::DotProduct<N, RealType> pObj(x, y);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class RealType>
inline RealType norm2(const Vector<RealType>& x)
{
    ParallelObjects::VectorSumSqr<1, RealType> pObj(x);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

template<Int N, class RealType>
inline RealType norm2(const Vector<VecX<N, RealType> >& x)
{
    ParallelObjects::VectorSumSqr<N, RealType> pObj(x);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, x.size()), pObj);

    return pObj.getResult();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class VectorType>
inline Vector<VectorType> add(const Vector<VectorType>& x, const Vector<VectorType>& y)
{
    __NT_REQUIRE(x.size() == y.size());
    Vector<VectorType> z(x.size());
    Scheduler::parallel_for(z.size(), [&](size_t i) { z[i] = x[i] + y[i]; });
    return z;
}

template<class VectorType>
inline Vector<VectorType> minus(const Vector<VectorType>& x, const Vector<VectorType>& y)
{
    __NT_REQUIRE(x.size() == y.size());
    Vector<VectorType> z(x.size());
    Scheduler::parallel_for(z.size(), [&](size_t i) { z[i] = x[i] - y[i]; });
    return z;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// saxpy (y=alpha*x+y)
//
template<class RealType, class VectorType>
inline void addScaled(RealType alpha, const Vector<VectorType>& x, Vector<VectorType>& y)
{
    Scheduler::parallel_for(x.size(), [&, alpha](size_t i) { y[i] += alpha * x[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// y=x+beta*y
//
template<class RealType, class VectorType>
inline void scaledAdd(RealType beta, const Vector<VectorType>& x, Vector<VectorType>& y)
{
    Scheduler::parallel_for(x.size(), [&, beta](size_t i) { y[i] = beta * y[i] + x[i]; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// x *= alpha
//
template<class RealType, class VectorType>
inline void scale(RealType alpha, Vector<VectorType>& x)
{
    Scheduler::parallel_for(x.size(), [&, alpha](size_t i) { x[i] *= alpha; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// y = x * alpha
template<class RealType, class VectorType>
inline Vector<VectorType> multiply(RealType alpha, const Vector<VectorType>& x)
{
    Vector<VectorType> y(x.size());
    Scheduler::parallel_for(x.size(), [&, alpha](size_t i) { y[i] = x[i] * alpha; });
    return y;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
}   // end namespace ParallelBLAS

