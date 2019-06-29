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

#include <LibCommon/Utils/Formatters.h>
#include <debugbreak.h>

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <csignal>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
namespace NTCodeBase {
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_SUPPORT_DOUBLE_NUMBER

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#if defined(_WIN32) || defined(_WIN64)
#   define NT_IN_WINDOWS_OS
#elif defined(__APPLE__)
#   define NT_IN_OSX
#elif defined(linux) || defined(__linux__)
#   define NT_IN_LINUX_OS
#endif

#ifndef NT_UNUSED
#   define NT_UNUSED(x) ((void)x)
#endif

#ifdef NT_IN_WINDOWS_OS
//  Exclude rarely-used stuff from Windows headers
#   define WIN32_LEAN_AND_MEAN
#   ifndef NOMINMAX
#     define NOMINMAX
#   endif
#   define __func__ __FUNCTION__
#endif // NT_IN_WINDOWS_OS

#ifdef NT_IN_WINDOWS_OS
#   define NT_SPRINT sprintf_s
#   define NT_SSCAN  sscanf_s
#else
#   define NT_SPRINT sprintf
#   define NT_SSCAN  sscanf
#endif

#if defined(DEBUG) || defined(_DEBUG)
#  define NT_DEBUG
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef NT_IN_WINDOWS_OS
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
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// signals
#define NT_RAISE_TERMINATION_SIGNAL std::raise(SIGTERM);

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Constants
#ifndef M_PI
# define M_PI 3.1415926535897932384626433832795028841971694
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// Utility
#define STR(x)            #x
#define NT_TOTO_STRING(x) STR(x)

#ifdef NT_IN_WINDOWS_OS
#define NT_ALIGN16 _MM_ALIGN16
#else
#define NT_ALIGN16 __attribute__((aligned(16)))
#endif

#ifdef NT_DEBUG
#  define NT_DEBUG_BREAK_OR_TERMINATE debug_break();
#else
#  define NT_DEBUG_BREAK_OR_TERMINATE NT_RAISE_TERMINATION_SIGNAL
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifdef NT_IN_WINDOWS_OS
#ifdef QT_CORE_LIB
#  define NT_COMPILER_MESSAGE(msg) \
    __pragma(message("\033[38;5;214m+++>" msg "\033[0m"))

#  define NT_TODO \
    __pragma(message("\033[38;5;214m+++>TODO: => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") \033[0m"))

#  define NT_TODO_MSG(msg) \
    __pragma(message("\033[38;5;214m+++>TODO: " msg " => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") \033[0m"))
#else
#  define NT_COMPILER_MESSAGE(msg) \
    __pragma(message("+++>" msg))

#  define NT_TODO \
    __pragma(message("+++>TODO: => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") "))

#  define NT_TODO_MSG(msg) \
    __pragma(message("+++>TODO: " msg " => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") "))
#endif
#else // not NT_IN_WINDOWS_OS
#  ifdef QT_CORE_LIB
#    define PRAGMA_MESSAGE(x) _Pragma(#x)
#    define NT_COMPILER_MESSAGE(msg) \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>" msg "\033[0m")

#    define NT_TODO \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>TODO: => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") \033[0m")

#    define NT_TODO_MSG(msg) \
    PRAGMA_MESSAGE(message "\033[38;5;214m+++>TODO: " msg " => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") \033[0m")
#  else
#    define NT_COMPILER_MESSAGE(msg) \
    PRAGMA_MESSAGE(message "+++>" msg)

#    define NT_TODO \
    PRAGMA_MESSAGE(message "+++>TODO: => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") ")

#    define NT_TODO_MSG(msg) \
    PRAGMA_MESSAGE(message "+++>TODO: " msg " => " __FILE__ "(" NT_TOTO_STRING(__LINE__) ") ")
#  endif
#endif
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_PRINT_LINE                           \
    {                                           \
        printf("%d: %s\n", __LINE__, __FILE__); \
        fflush(stdout);                         \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_PRINT_EXP(x)                                                                       \
    {                                                                                         \
        std::stringstream ss;                                                                 \
        ss << "Printing at line: " << __LINE__ << ", file: " << __FILE__ << ":" << std::endl; \
        ss << "    " << #x << ": " << x;                                                      \
        printf("%s\n", ss.str().c_str());                                                     \
        fflush(stdout);                                                                       \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_PRINT_LOCATION                                     \
    {                                                         \
        std::stringstream ss;                                 \
        ss << "Function: " << __func__ << std::endl;          \
        ss << "Line: " << __LINE__ << ", file: " << __FILE__; \
        printf("%s\n", ss.str().c_str());                     \
        fflush(stdout);                                       \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifndef NT_INFO
#   define NT_INFO(info)       \
    {                          \
        fprintf(stderr, info); \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#ifndef NT_ERROR
#   define NT_ERROR(err)       \
    {                          \
        fprintf(stderr, err);  \
        fprintf(stderr, "\n"); \
        NT_PRINT_LOCATION      \
    }
#endif

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_DIE(err)                 \
    {                               \
        NT_ERROR(err)               \
        NT_RAISE_TERMINATION_SIGNAL \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#define NT_REQUIRE(condition)                                                 \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            printf("%s\n", erMsg.c_str());                                    \
            NT_PRINT_LOCATION                                                 \
                NT_DEBUG_BREAK_OR_TERMINATE                                   \
        }                                                                     \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_REQUIRE_MSG(condition, msg)                                        \
    {                                                                         \
        if(!(condition))                                                      \
        {                                                                     \
            String erMsg = String("Assertion failed: ") + String(#condition); \
            String rsMsg = String("Reason: ") + String(msg);                  \
            printf("%s\n%s\n", erMsg.c_str(), rsMsg.c_str());                 \
            NT_PRINT_LOCATION                                                 \
                NT_DEBUG_BREAK_OR_TERMINATE                                   \
        }                                                                     \
    }

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_UNIMPLEMENTED_FUNC          NT_DIE("Called function is unimplemented.")
#define NT_CALLED_TO_WRONG_PLACE       NT_DIE("This function should not be reached.")
#define NT_DENIED_CALL_TO_BASE_CLASS   NT_DIE("This function should not be called from base class.")
#define NT_DENIED_SWITCH_DEFAULT_VALUE NT_DIE("Invalid default value in switch statement.")
#define NT_DIE_UNKNOWN_ERROR           NT_DIE("An unknown error has occured...")

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
// test
#define NT_PERORMANCE_TEST_BEGIN \
    Timer testTimer;             \
    testTimer.tick();

#define NT_PERORMANCE_TEST_END(funcName) \
    testTimer.tock();                    \
    printf("Test %s finished. Time: %s\n", funcName, testTimer.getRunTime().c_str());

#define NT_PERORMANCE_TEST_END_WITH_RUN_TIME(funcName, runTime)                       \
    runTime = testTimer.tock();                                                       \
    printf("Test %s finished. Time: %s\n", funcName, testTimer.getRunTime().c_str()); \

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_REQUIRE_EQUAL(a, b)                     NT_REQUIRE_MSG(a == b, "Numbers are not equal.");
#define NT_REQUIRE_APPROX_NUMBERS(a, b, threshold) NT_REQUIRE_MSG(threshold > std::abs(a - b), "Numbers are not equal.");

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_TO_CSTRING(x)                           Formatters::toString7(x).c_str()
#define NT_TO_CSTRING_2(x)                         Formatters::toString2(x).c_str()
#define NT_TO_CSTRING_5(x)                         Formatters::toString5(x).c_str()
#define NT_TO_CSTRING_7(x)                         Formatters::toString7(x).c_str()

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define NT_TYPE_ALIAS                                                               \
    using Vec1r            = Vec1<Real_t>;                                          \
    using Vec2r            = Vec2<Real_t>;                                          \
    using Vec3r            = Vec3<Real_t>;                                          \
    using Vec4r            = Vec4<Real_t>;                                          \
    using VecNi            = VecX<N, Int>;                                          \
    using VecNui           = VecX<N, UInt>;                                         \
    using VecN             = VecX<N, Real_t>;                                       \
    using VecNp1           = VecX<N + 1, Real_t>;                                   \
    using Mat2x2r          = Mat2x2<Real_t>;                                        \
    using Mat2x3r          = Mat2x3<Real_t>;                                        \
    using Mat3x2r          = Mat3x2<Real_t>;                                        \
    using Mat3x3r          = Mat3x3<Real_t>;                                        \
    using Mat4x4r          = Mat4x4<Real_t>;                                        \
    using MatNxN           = MatXxX<N, Real_t>;                                     \
    using MatNp1xNp1       = MatXxX<N + 1, Real_t>;                                 \
    using StdVT_Vec2r      = StdVT_Vec2<Real_t>;                                    \
    using StdVT_Vec3r      = StdVT_Vec3<Real_t>;                                    \
    using StdVT_Vec4r      = StdVT_Vec4<Real_t>;                                    \
    using StdVT_VecN       = StdVT_VecX<N, Real_t>;                                 \
    using StdVT_VecNp1     = StdVT_VecX<N + 1, Real_t>;                             \
    using StdVT_MatNxN     = StdVT_MatXxX<N, Real_t>;                               \
    using StdVT_MatNp1xNp1 = StdVT_MatXxX<N + 1, Real_t>;                           \
    using StdVT_Realt      = StdVT<Real_t>;                                         \
    static constexpr auto TinyReal() { return std::numeric_limits<Real_t>::min(); } \
    static constexpr auto HugeReal() { return std::numeric_limits<Real_t>::max(); }

#ifdef NT_SUPPORT_DOUBLE_NUMBER
#  define NT_INSTANTIATE_CLASS_COMMON_TYPES(ClassName) \
    template class ClassName<float>;                   \
    template class ClassName<double>;

#  define NT_INSTANTIATE_STRUCT_COMMON_TYPES(StructName) \
    template struct StructName<float>;                   \
    template struct StructName<double>;

#  define NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ClassName) \
    template class ClassName<2, float>;                               \
    template class ClassName<3, float>;                               \
    template class ClassName<2, double>;                              \
    template class ClassName<3, double>;

#  define NT_INSTANTIATE_STRUCT_COMMON_DIMENSIONS_AND_TYPES(StructName) \
    template struct StructName<2, float>;                               \
    template struct StructName<3, float>;                               \
    template struct StructName<2, double>;                              \
    template struct StructName<3, double>;
#else
#  define NT_INSTANTIATE_CLASS_COMMON_TYPES(ClassName) \
    template class ClassName<float>;

#  define NT_INSTANTIATE_STRUCT_COMMON_TYPES(StructName) \
    template struct StructName<float>;

#  define NT_INSTANTIATE_CLASS_COMMON_DIMENSIONS_AND_TYPES(ClassName) \
    template class ClassName<2, float>;                               \
    template class ClassName<3, float>;

#  define NT_INSTANTIATE_STRUCT_COMMON_DIMENSIONS_AND_TYPES(StructName) \
    template struct StructName<2, float>;                               \
    template struct StructName<3, float>;
#endif
//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
} // end namespace NTCodeBase
