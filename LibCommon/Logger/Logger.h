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

#include <LibCommon/CommonSetup.h>
#include <LibCommon/Timer/Timer.h>

#include <spdlog/spdlog.h>

#include <chrono>
#include <map>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class Logger
{
    using Clock = std::chrono::system_clock;
public:
    using LogLevel = spdlog::level::level_enum;
    ////////////////////////////////////////////////////////////////////////////////
    Logger(const String& loggerName, const String& rootPath, bool bLog2Console = true, bool bLog2File = false, LogLevel logLevel = LogLevel::trace);
    ~Logger();
    ////////////////////////////////////////////////////////////////////////////////
    void printLog(const String& s, LogLevel logLevel = LogLevel::info)
    {
        if(m_bLog2Console) { m_ConsoleLogger->log(logLevel, s); }
        if(m_bLog2File) { m_FileLogger->log(logLevel, s); }
    }

    void printLogPadding(const String& s, LogLevel logLevel       = LogLevel::info);
    void printLogIndent(const String& s, UInt indentLevel         = 1, LogLevel logLevel = LogLevel::info);
    void printLogPaddingIndent(const String& s,  UInt indentLevel = 1, LogLevel logLevel = LogLevel::info);

    void printCenterAligned(const String& s, char paddingChar      = Logger::s_PrefixPadding, LogLevel logLevel = LogLevel::info);
    void printTextBox(const String& s, LogLevel logLevel           = LogLevel::info);
    void printTextBox(const StdVT<String>& strs, LogLevel logLevel = LogLevel::info);

    void newLine(LogLevel logLevel       = LogLevel::info) { printLog("", logLevel); }
    void separatorLine(LogLevel logLevel = LogLevel::info) { printCenterAligned("", '=', logLevel); }
    ////////////////////////////////////////////////////////////////////////////////
    template<typename... Args> void printLogIf(bool bCondition, Args&& ... args) { if(bCondition) { printLog(std::forward<Args>(args)...); } }
    template<typename... Args> void printLogIndentIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogIndent(std::forward<Args>(args)...); } }
    template<typename... Args> void printLogPaddingIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogPadding(std::forward<Args>(args)...); } }
    template<typename... Args> void printLogPaddingIndentIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogPaddingIndent(std::forward<Args>(args)...); } }
    template<typename... Args> void printCenterAlignedIf(bool bCondition, Args&& ... args) { if(bCondition) { printCenterAligned(std::forward<Args>(args)...); } }
    template<typename... Args> void printTextBoxIf(bool bCondition, Args&& ... args) { if(bCondition) { printTextBox(std::forward<Args>(args)...); } }
    template<typename... Args> void newLineIf(bool bCondition, Args&& ... args) { if(bCondition) { newLine(std::forward<Args>(args)...); } }
    template<typename... Args> void separatorLineIf(bool bCondition, Args&& ... args) { if(bCondition) { separatorLine(std::forward<Args>(args)...); } }
    ////////////////////////////////////////////////////////////////////////////////
    void printDebug(const String& s) { printLog(s, LogLevel::debug); }
    void printDebugIndent(const String& s, UInt indentLevel = 1) { printLogIndent(s, indentLevel, LogLevel::debug); }
    void printWarning(const String& s) { printLog(s, LogLevel::warn); }
    void printWarningIndent(const String& s, UInt indentLevel = 1) { printLogIndent(s, indentLevel, LogLevel::warn); }
    void printError(const String& s) { printLog(s, LogLevel::err); }
    void printErrorIndent(const String& s, UInt indentLevel = 1) { printLogIndent(s, indentLevel, LogLevel::err); }
    ////////////////////////////////////////////////////////////////////////////////
    template<class Function> void printRunTime(const char* caption, const Function& function, LogLevel logLevel = LogLevel::info)
    {
        printLog(Timer::getRunTime<Function>(caption, function), logLevel);
    }

    template<class Function> void printRunTimeIndent(const char* caption, const Function& function, UInt indentLevel = 1, LogLevel logLevel = LogLevel::info)
    {
        printLogIndent(Timer::getRunTime<Function>(caption, function), indentLevel, logLevel);
    }

    template<class Function> void printRunTimeIf(const char* caption, const Function& function, LogLevel logLevel = LogLevel::info)
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIf(bResult, timer.getRunTime(caption), logLevel);
    }

    template<class Function> void printRunTimeIndentIf(const char* caption, const Function& function, UInt indentLevel = 1, LogLevel logLevel = LogLevel::info)
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIndentIf(bResult, timer.getRunTime(caption), indentLevel, logLevel);
    }

    ////////////////////////////////////////////////////////////////////////////////
    void flush();
    void printMemoryUsage(LogLevel logLevel  = LogLevel::info);
    void printTotalRunTime(LogLevel logLevel = LogLevel::info);
    ////////////////////////////////////////////////////////////////////////////////
    void        cleanup(int signal = 0);
    static void signalHandler(int signal);
    ////////////////////////////////////////////////////////////////////////////////
    inline const static String   s_Wrapper { "||" };
    inline const static char     s_PrefixPadding { ' ' };
    inline const static char     s_SuffixPadding { '*' };
    inline const static UInt     s_IndentSize { 4u };
    inline const static size_t   s_PaddingMaxSize { 100u };
    inline const static UInt     s_BufferLength { 128u };
    inline static StdVT<Logger*> s_Instances {};
private:
    String getTotalRunTime();
    ////////////////////////////////////////////////////////////////////////////////
    Clock::time_point         m_StartupTime {};
    SharedPtr<spdlog::logger> m_ConsoleLogger;
    SharedPtr<spdlog::logger> m_FileLogger;
    bool                      m_bLog2Console;
    bool                      m_bLog2File;
};
