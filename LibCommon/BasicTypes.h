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

#include <limits>
#include <vector>
#include <map>
#include <set>
#include <cstdint>
#include <string>
#include <memory>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Basic types
using Int8  = int8_t;
using Int16 = int16_t;
using Int   = int32_t;
using Int32 = int32_t;
using Int64 = int64_t;

////////////////////////////////////////////////////////////////////////////////
using UChar  = unsigned char;
using UInt8  = uint8_t;
using UInt16 = uint16_t;
using UInt   = uint32_t;
using UInt32 = uint32_t;
using UInt64 = uint64_t;

////////////////////////////////////////////////////////////////////////////////
using String = std::string;

////////////////////////////////////////////////////////////////////////////////
using PairInt8  = std::pair<Int8, Int8>;
using PairInt16 = std::pair<Int16, Int16>;
using PairInt32 = std::pair<Int, Int>;
using PairInt   = std::pair<Int, Int>;
using PairInt64 = std::pair<Int64, Int64>;

////////////////////////////////////////////////////////////////////////////////
using PairUInt8  = std::pair<UInt8, UInt8>;
using PairUInt16 = std::pair<UInt16, UInt16>;
using PairUInt32 = std::pair<UInt, UInt>;
using PairUInt   = std::pair<UInt, UInt>;
using PairUInt64 = std::pair<UInt64, UInt64>;

using FuncPtr = const void*;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<class Type> using SharedPtr = std::shared_ptr<Type>;
template<class Type> using UniquePtr = std::unique_ptr<Type>;

template<class Type> using Set          = std::set<Type>;
template<class Type> using StdVT        = std::vector<Type>;
template<class Type> using StdVT_StdVec = StdVT<StdVT<Type>>;
template<class K, class V> using Map    = std::map<K, V>;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
using StdVT_Int8  = StdVT<Int8>;
using StdVT_Int16 = StdVT<Int16>;
using StdVT_Int   = StdVT<Int>;
using StdVT_Int32 = StdVT<Int>;
using StdVT_Int64 = StdVT<Int64>;

using StdVT_UInt8  = StdVT<UInt8>;
using StdVT_UInt16 = StdVT<UInt16>;
using StdVT_UInt   = StdVT<UInt>;
using StdVT_UInt32 = StdVT<UInt>;
using StdVT_UInt64 = StdVT<UInt64>;

using StdVT_Char   = StdVT<char>;
using StdVT_UChar  = StdVT<unsigned char>;
using StdVT_Float  = StdVT<float>;
using StdVT_Double = StdVT<double>;
using StdVT_String = StdVT<String>;

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Conversion operators
float constexpr operator"" _f(long double x)
{
    return static_cast<float>(x);
}

float constexpr operator"" _f(unsigned long long int x)
{
    return static_cast<float>(x);
}

double constexpr operator"" _d(long double x)
{
    return static_cast<double>(x);
}

double constexpr operator"" _d(unsigned long long int x)
{
    return static_cast<double>(x);
}

Int32 constexpr operator"" _int(unsigned long long int x)
{
    return static_cast<Int32>(x);
}

UInt32 constexpr operator"" _uint(unsigned long long int x)
{
    return static_cast<UInt32>(x);
}

UInt64 constexpr operator"" _uint64(unsigned long long int x)
{
    return static_cast<UInt64>(x);
}

std::size_t constexpr operator"" _sz(unsigned long long int x)
{
    return static_cast<std::size_t>(x);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#include <limits>
#undef min
#undef max
template<class T> constexpr auto MEpsilon() { return std::numeric_limits<T>::epsilon(); }
template<class T> constexpr auto Tiny() { return std::numeric_limits<T>::min(); }
template<class T> constexpr auto Huge() { return std::numeric_limits<T>::max(); }
