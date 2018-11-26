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
#include <LibCommon/CommonSetup.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace AtomicOperations {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T, class Function>
inline void atomicOp(T& target, T operand, Function&& f) {
    std::atomic<T>& tgt = *((std::atomic<T>*)&target);

    T cur_val = target;
    T new_val;
    do {
        new_val = f(cur_val, operand);
    } while(!tgt.compare_exchange_weak(cur_val, new_val));
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
inline void add(T& target, T operand) {
    atomicOp(target, operand, [](T a, T b) { return a + b; });
}

template<class T>
inline void add(Vec2<T>& target, const Vec2<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a + b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a + b; });
}

template<class T>
inline void add(Vec3<T>& target, const Vec3<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a + b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a + b; });
    atomicOp(target[2], operand[2], [](T a, T b) { return a + b; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
inline void subtract(T& target, T operand) {
    atomicOp(target, operand, [](T a, T b) { return a - b; });
}

template<class T>
inline void subtract(Vec2<T>& target, const Vec2<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a - b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a - b; });
}

template<class T>
inline void subtract(Vec3<T>& target, const Vec3<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a - b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a - b; });
    atomicOp(target[2], operand[2], [](T a, T b) { return a - b; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
inline void multiply(T& target, T operand) {
    atomicOp(target, operand, [](T a, T b) { return a * b; });
}

template<class T>
inline void multiply(Vec2<T>& target, const Vec2<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a * b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a * b; });
}

template<class T>
inline void multiply(Vec3<T>& target, const Vec3<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a * b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a * b; });
    atomicOp(target[2], operand[2], [](T a, T b) { return a * b; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class T>
inline void divide(T& target, T operand) {
    atomicOp(target, operand, [](T a, T b) { return a / b; });
}

template<class T>
inline void divide(Vec2<T>& target, const Vec2<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a / b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a / b; });
}

template<class T>
inline void divide(Vec3<T>& target, const Vec3<T>& operand) {
    atomicOp(target[0], operand[0], [](T a, T b) { return a / b; });
    atomicOp(target[1], operand[1], [](T a, T b) { return a / b; });
    atomicOp(target[2], operand[2], [](T a, T b) { return a / b; });
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace ParallelOperations
