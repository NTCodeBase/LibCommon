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

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//#define __NT_SUPPORT_DOUBLE_NUMBER

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#if defined(_WIN32) || defined(_WIN64)
#   define __NT_WINDOWS_OS__
#elif defined(__APPLE__)
#   define __NT_OSX__
#elif defined(linux) || defined(__linux__)
#   define __NT_LINUX_OS__
#endif

#ifndef __NT_UNUSED
#   define __NT_UNUSED(x) ((void)x)
#endif

#ifdef __NT_WINDOWS_OS__
//  Exclude rarely-used stuff from Windows headers
#   define WIN32_LEAN_AND_MEAN
#   ifndef NOMINMAX
#     define NOMINMAX
#   endif
#   define __func__ __FUNCTION__
#   include <Windows.h>
#   include <exception>
#   include <string>

inline void throwIfFailed(HRESULT hr)
{
    if(FAILED(hr)) {
        throw std::exception(std::to_string(hr).c_str());
    }
}

#endif  // __NT_WINDOWS_OS__

#ifdef __NT_WINDOWS_OS__
#   define __NT_SPRINT sprintf_s
#   define __NT_SSCAN  sscanf_s
#else
#   define __NT_SPRINT sprintf
#   define __NT_SSCAN  sscanf
#endif

#if defined(DEBUG) || defined(_DEBUG)
#  define __NT_DEBUG__
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __NT_WINDOWS_OS__
// Enable memory leak detection
#  ifdef _DEBUG
#    define _CRTDBG_MAP_ALLOC
#    include <stdlib.h>
#    include <crtdbg.h>
#    define DEBUG_NEW           new(_NORMAL_BLOCK, __FILE__, __LINE__)
#    define REPORT_MEMORY_LEAKS _CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#  else
#    define REPORT_MEMORY_LEAKS
#  endif
#else
#    define REPORT_MEMORY_LEAKS
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#include <csignal>
#define __NT_RAISE_TERMINATION_SIGNAL \
    {                                 \
        std::raise(SIGTERM);          \
    }

#define __NT_RAISE_ABORT_SIGNAL \
    {                           \
        std::raise(SIGABRT);    \
    }
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Constants
#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795028841971694
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Utility
#define _STR(x)           #x
#define __NT_TO_STRING(x) _STR(x)

#ifdef __NT_WINDOWS_OS__
#define __NT_FORCE_INLINE __forceinline
#else
#define __NT_FORCE_INLINE __attribute__((always_inline))
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __NT_WINDOWS_OS__
#ifdef QT_CORE_LIB
#  define __NT_COMPILER_MESSAGE(msg) \
    __pragma(message("\033[38;5;214m+++>" msg "\033[0m"))

#  define __NT_TODO \
    __pragma(message("\033[38;5;214m+++>TODO: => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") \033[0m"))

#  define __NT_TODO_MSG(msg) \
    __pragma(message("\033[38;5;214m+++>TODO: " msg " => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") \033[0m"))
#else
#  define __NT_COMPILER_MESSAGE(msg) \
    __pragma(message("+++>" msg))

#  define __NT_TODO \
    __pragma(message("+++>TODO: => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") "))

#  define __NT_TODO_MSG(msg) \
    __pragma(message("+++>TODO: " msg " => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") "))
#endif
#else // not __NT_WINDOWS_OS__
#  ifdef QT_CORE_LIB
#    define PRAGMA_MESSAGE(x) _Pragma(#x)
#    define __NT_COMPILER_MESSAGE(msg) \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>" msg "\033[0m")

#    define __NT_TODO \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>TODO: => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") \033[0m")

#    define __NT_TODO_MSG(msg) \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>TODO: " msg " => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") \033[0m")
#  else
#    define __NT_COMPILER_MESSAGE(msg) \
    PRAGMA_MESSAGE(message "+++>" msg)

#    define __NT_TODO \
    PRAGMA_MESSAGE(message "+++>TODO: => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") ")

#    define __NT_TODO_MSG(msg) \
    PRAGMA_MESSAGE(message "+++>TODO: " msg " => " __FILE__ "(" __NT_TO_STRING(__LINE__) ") ")
#  endif
#endif
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_PRINT_LINE                         \
    {                                           \
        printf("%d: %s\n", __LINE__, __FILE__); \
        fflush(stdout);                         \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_PRINT_EXP(x)                                                                     \
    {                                                                                         \
        std::stringstream ss;                                                                 \
        ss << "Printing at line: " << __LINE__ << ", file: " << __FILE__ << ":" << std::endl; \
        ss << "    " << #x << ": " << x;                                                      \
        printf("%s\n", ss.str().c_str());                                                     \
        fflush(stdout);                                                                       \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_PRINT_LOCATION                                   \
    {                                                         \
        std::stringstream ss;                                 \
        ss << "Function: " << __func__ << std::endl;          \
        ss << "Line: " << __LINE__ << ", file: " << __FILE__; \
        printf("%s\n", ss.str().c_str());                     \
        fflush(stdout);                                       \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifndef __NT_INFO
#   define __NT_INFO(info)     \
    {                          \
        fprintf(stderr, info); \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifndef __NT_ERROR
#   define __NT_ERROR(err)    \
    {                         \
        fprintf(stderr, err); \
        __NT_PRINT_LOCATION   \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_DIE(err)                 \
    {                                 \
        __NT_ERROR(err)               \
        __NT_RAISE_TERMINATION_SIGNAL \
    }
//exit(EXIT_FAILURE);
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __NT_DEBUG__
#  define __NT_REQUIRE(condition)                                             \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            printf("%s\n", erMsg.c_str());                                    \
            __NT_PRINT_LOCATION                                               \
            assert(false);                                                    \
        }                                                                     \
    }
#else
#  define __NT_REQUIRE(condition)                                             \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            printf("%s\n", erMsg.c_str());                                    \
            __NT_PRINT_LOCATION                                               \
            exit(EXIT_FAILURE);                                               \
        }                                                                     \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef __NT_DEBUG__
#  define __NT_REQUIRE_MSG(condition, msg)                                    \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            String rsMsg = String("Reason: ") + String(msg);                  \
            printf("%s\n%s\n", erMsg.c_str(), rsMsg.c_str());                 \
            __NT_PRINT_LOCATION                                               \
            assert(false);                                                    \
        }                                                                     \
    }
#else
#  define __NT_REQUIRE_MSG(condition, msg)                                    \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            String rsMsg = String("Reason: ") + String(msg);                  \
            printf("%s\n%s\n", erMsg.c_str(), rsMsg.c_str());                 \
            __NT_PRINT_LOCATION                                               \
            exit(EXIT_FAILURE);                                               \
        }                                                                     \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_CHECK_ERROR(condition, msg)                                      \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            String rsMsg = String("Reason: ") + String(msg);                  \
            printf("%s\n%s\n", erMsg.c_str(), rsMsg.c_str());                 \
            __NT_PRINT_LOCATION                                               \
        }                                                                     \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define __NT_UNIMPLEMENTED_FUNC          __NT_DIE("Called function is unimplemented.")
#define __NT_CALLED_TO_WRONG_PLACE       __NT_DIE("This function should not be reached.")
#define __NT_DENIED_CALL_TO_BASE_CLASS   __NT_DIE("This function should not be called from base class.")
#define __NT_DENIED_SWITCH_DEFAULT_VALUE __NT_DIE("Invalid default value in switch statement.")
#define __NT_DIE_UNKNOWN_ERROR           __NT_DIE("An unknown error has occured...")

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// test
#define __NT_PERORMANCE_TEST_BEGIN \
    Timer testTimer;               \
    testTimer.tick();

#define __NT_PERORMANCE_TEST_END(funcName) \
    testTimer.tock();                      \
    printf("Test %s finished. Time: %s\n", funcName, testTimer.getRunTime().c_str());

#define __NT_PERORMANCE_TEST_END_WITH_RUN_TIME(funcName, runTime)                     \
    runTime = testTimer.tock();                                                       \
    printf("Test %s finished. Time: %s\n", funcName, testTimer.getRunTime().c_str()); \

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// data macros
#define __NT_REQUIRE_EQUAL(a, b)                     __NT_REQUIRE_MSG(a == b, "Numbers are not equal.");
#define __NT_REQUIRE_APPROX_NUMBERS(a, b, threshold) __NT_REQUIRE_MSG(threshold > fabs(a - b), "Numbers are not equal.");
#define __NT_TO_CSTRING(x)                           Formatters::toString(x).c_str()
#define __NT_TO_CSTRING_5(x)                         Formatters::toString5(x).c_str()
#define __NT_TO_CSTRING_7(x)                         Formatters::toString7(x).c_str()

#define __NT_TYPE_ALIASING                                                            \
    using VecN             = VecX<N, RealType>;                                       \
    using VecNp1           = VecX<N + 1, RealType>;                                   \
    using MatNxN           = MatXxX<N, RealType>;                                     \
    using MatNp1xNp1       = MatXxX<N + 1, RealType>;                                 \
    using StdVT_VecN       = StdVT_VecX<N, RealType>;                                 \
    using StdVT_VecNp1     = StdVT_VecX<N + 1, RealType>;                             \
    using StdVT_MatNxN     = StdVT_MatXxX<N, RealType>;                               \
    using StdVT_MatNp1xNp1 = StdVT_MatXxX<N + 1, RealType>;                           \
    using StdVT_RealType   = StdVT<RealType>;                                         \
    static constexpr auto TinyReal() { return std::numeric_limits<RealType>::min(); } \
    static constexpr auto HugeReal() { return std::numeric_limits<RealType>::max(); }

#ifdef __NT_SUPPORT_DOUBLE_NUMBER
#  define __NT_INSTANTIATE_CLASS_COMMON_TYPES(ClassName) \
    template class ClassName<float>;                     \
    template class ClassName<double>;

#  define __NT_INSTANTIATE_STRUCT_COMMON_TYPES(StructName) \
    template struct StructName<float>;                     \
    template struct StructName<double>;

#  define __NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ClassName) \
    template class ClassName<2, float>;                                 \
    template class ClassName<3, float>;                                 \
    template class ClassName<2, double>;                                \
    template class ClassName<3, double>;

#  define __NT_INSTANTIATE_STRUCT_COMMON_DIMENSIONS_AND_TYPES(StructName) \
    template struct StructName<2, float>;                                 \
    template struct StructName<3, float>;                                 \
    template struct StructName<2, double>;                                \
    template struct StructName<3, double>;
#else
#  define __NT_INSTANTIATE_CLASS_COMMON_TYPES(ClassName) \
    template class ClassName<float>;

#  define __NT_INSTANTIATE_STRUCT_COMMON_TYPES(StructName) \
    template struct StructName<float>;

#  define __NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ClassName) \
    template class ClassName<2, float>;                                 \
    template class ClassName<3, float>;

#  define __NT_INSTANTIATE_STRUCT_COMMON_DIMENSIONS_AND_TYPES(StructName) \
    template struct StructName<2, float>;                                 \
    template struct StructName<3, float>;
#endif
