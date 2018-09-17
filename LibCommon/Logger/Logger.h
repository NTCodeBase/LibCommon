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
    Logger(const String& loggerName, const String& rootPath, bool bLog2Console = true, bool bLog2File = false,
           LogLevel consoleLogLevel = LogLevel::trace, LogLevel fileLogLevel = LogLevel::trace);
    ~Logger();
    ////////////////////////////////////////////////////////////////////////////////
    void printLog(const String& s, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info)
    {
        if(m_bLog2Console) { m_ConsoleLogger->log(consoleLogLevel, s); }
        if(m_bLog2File) { m_FileLogger->log(fileLogLevel, s); }
    }

    void printLogPadding(const String& s, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    void printLogIndent(const String& s, UInt indentLevel          = 1, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    void printLogPaddingIndent(const String& s,  UInt indentLevel  = 1, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);

    void printCenterAligned(const String& s, char paddingChar             = Logger::s_PrefixPadding, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    void printTextBox(const String& s, LogLevel consoleLogLevel           = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    void printTextBox(const StdVT<String>& strs, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);

    void newLine(LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info) { printLog("", consoleLogLevel, fileLogLevel); }
    void separatorLine(int stype          = 0, LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    ////////////////////////////////////////////////////////////////////////////////
    template<class... Args> void printLogIf(bool bCondition, Args&& ... args) { if(bCondition) { printLog(std::forward<Args>(args)...); } }
    template<class... Args> void printLogIndentIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogIndent(std::forward<Args>(args)...); } }
    template<class... Args> void printLogPaddingIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogPadding(std::forward<Args>(args)...); } }
    template<class... Args> void printLogPaddingIndentIf(bool bCondition, Args&& ... args) { if(bCondition) { printLogPaddingIndent(std::forward<Args>(args)...); } }
    template<class... Args> void printCenterAlignedIf(bool bCondition, Args&& ... args) { if(bCondition) { printCenterAligned(std::forward<Args>(args)...); } }
    template<class... Args> void printTextBoxIf(bool bCondition, Args&& ... args) { if(bCondition) { printTextBox(std::forward<Args>(args)...); } }
    template<class... Args> void newLineIf(bool bCondition, Args&& ... args) { if(bCondition) { newLine(std::forward<Args>(args)...); } }
    template<class... Args> void separatorLineIf(bool bCondition, Args&& ... args) { if(bCondition) { separatorLine(std::forward<Args>(args)...); } }
    ////////////////////////////////////////////////////////////////////////////////
    void printDebug(const String& s) { printLog(s, LogLevel::debug, LogLevel::debug); }
    void printDebugIndent(const String& s, UInt indentLevel = 1) { printLogIndent(s, indentLevel, LogLevel::debug, LogLevel::debug); }
    void printWarning(const String& s) { printLogPadding(s, LogLevel::warn, LogLevel::warn); }
    void printWarningIndent(const String& s, UInt indentLevel = 1) { printLogPaddingIndent(s, indentLevel, LogLevel::warn, LogLevel::warn); }
    void printError(const String& s) { printLogPadding(s, LogLevel::err, LogLevel::err); }
    void printErrorIndent(const String& s, UInt indentLevel = 1) { printLogPaddingIndent(s, indentLevel, LogLevel::err, LogLevel::err); }
    ////////////////////////////////////////////////////////////////////////////////
    template<class Function, class... Args>
    void printRunTime(const char* caption, Function&& function, Args&& ... args)
    {
        printLog(Timer::getRunTime<Function>(caption, std::forward<Function>(function)), std::forward<Args>(args)...);
    }

    template<class Function, class... Args>
    void printRunTimeIndent(const char* caption, Function&& function, UInt indentLevel = 1, Args&& ... args)
    {
        printLogIndent(Timer::getRunTime<Function>(caption, std::forward<Function>(function)), indentLevel, std::forward<Args>(args)...);
    }

    template<class Function, class... Args>
    void printRunTimeIf(const char* caption, Function&& function, Args&& ... args)
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIf(bResult, timer.getRunTime(caption), std::forward<Args>(args)...);
    }

    template<class Function, class... Args>
    void printRunTimeIndentIf(const char* caption, Function&& function, UInt indentLevel = 1, Args&& ... args)
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIndentIf(bResult, timer.getRunTime(caption), indentLevel, std::forward<Args>(args)...);
    }

    ////////////////////////////////////////////////////////////////////////////////
    template<class... Args> void printRunTimeDebug(Args&& ... args) { printRunTime(std::forward<Args>(args)..., LogLevel::debug, LogLevel::debug); }
    template<class... Args> void printRunTimeIfDebug(Args&& ... args) { printRunTimeIf(std::forward<Args>(args)..., LogLevel::debug, LogLevel::debug); }

    template<class Function>
    void printRunTimeIndentDebug(const char* caption, Function&& function, UInt indentLevel = 1)
    {
        printRunTimeIndent(caption, std::forward<Function>(function), indentLevel, LogLevel::debug, LogLevel::debug);
    }

    template<class Function>
    void printRunTimeIndentIfDebug(const char* caption, Function&& function, UInt indentLevel = 1)
    {
        printRunTimeIndentIf(caption, std::forward<Function>(function), indentLevel, LogLevel::debug, LogLevel::debug);
    }

    ////////////////////////////////////////////////////////////////////////////////
    void flush();
    void printMemoryUsage(LogLevel consoleLogLevel  = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    void printTotalRunTime(LogLevel consoleLogLevel = LogLevel::info, LogLevel fileLogLevel = LogLevel::info);
    ////////////////////////////////////////////////////////////////////////////////
    void        cleanup(int signal = 0);
    static void signalHandler(int signal);
    ////////////////////////////////////////////////////////////////////////////////
    inline const static String   s_Wrapper { "||" };
    inline const static char     s_PrefixPadding { ' ' };
    inline const static char     s_SuffixPadding { '*' };
    inline const static UInt     s_IndentSize { 4u };
    inline const static size_t   s_PaddingMaxSize { 120u };
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
