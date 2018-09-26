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

#include <atomic>
#include <LibCommon/CommonSetup.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace AtomicOperations
{
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T, class Function>
__NT_FORCE_INLINE void atomicOp(T& target, const T& operand, Function&& f)
{
    std::atomic<T>& tgt = *((std::atomic<T>*)&target);

    T cur_val = target;
    T new_val;
    do {
        new_val = f(cur_val, operand);
    } while(!tgt.compare_exchange_weak(cur_val, new_val));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
__NT_FORCE_INLINE void add(T& target, const T& operand)
{
    atomicOp(target, operand, [](T a, T b) { return a + b; });
}

template<Int N, class T>
__NT_FORCE_INLINE void add(VecX<N, T>& target, const VecX<N, T>& operand)
{
    if constexpr (N == 4) {
        atomicOp(target, operand, [](T a, T b) { return a + b; });
    } else {
        for(Int d = 0; d < N; ++d) {
            atomicOp(target[d], operand[d], [](T a, T b) { return a + b; });
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
__NT_FORCE_INLINE void subtract(T& target, const T& operand)
{
    atomicOp(target, operand, [](T a, T b) { return a - b; });
}

template<Int N, class T>
__NT_FORCE_INLINE void subtract(VecX<N, T>& target, const VecX<N, T>& operand)
{
    if constexpr (N == 4) {
        atomicOp(target, operand, [](T a, T b) { return a + b; });
    } else {
        for(Int d = 0; d < N; ++d) {
            atomicOp(target[d], operand[d], [](T a, T b) { return a - b; });
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
__NT_FORCE_INLINE void multiply(T& target, const T& operand)
{
    atomicOp(target, operand, [](T a, T b) { return a * b; });
}

template<Int N, class T>
__NT_FORCE_INLINE void multiply(VecX<N, T>& target, const VecX<N, T>& operand)
{
    if constexpr (N == 4) {
        atomicOp(target, operand, [](T a, T b) { return a + b; });
    } else {
        for(Int d = 0; d < N; ++d) {
            atomicOp(target[d], operand[d], [](T a, T b) { return a * b; });
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
__NT_FORCE_INLINE void divide(T& target, const T& operand)
{
    atomicOp(target, operand, [](T a, T b) { return a / b; });
}

template<Int N, class T>
__NT_FORCE_INLINE void divide(VecX<N, T>& target, const VecX<N, T>& operand)
{
    if constexpr (N == 4) {
        atomicOp(target, operand, [](T a, T b) { return a + b; });
    } else {
        for(Int d = 0; d < N; ++d) {
            atomicOp(target[d], operand[d], [](T a, T b) { return a / b; });
        }
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace ParallelOperations
