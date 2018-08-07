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

#include <CommonSetup.h>
#include <Utils/Timer.h>

#include <spdlog/spdlog.h>

#include <chrono>
#include <map>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
#define WRAPPER           "||"
#define PADDING           ' '
#define INDENT_SIZE       4
#define MAX_BUFFER_LENGTH 128

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class Logger
{
    using Clock = std::chrono::system_clock;
public:
    Logger(const String& loggerName, const String& rootPath, bool bLog2Console = true, bool bLog2File = false,
           spdlog::level::level_enum level                                     = spdlog::level::level_enum::trace);
    ~Logger();
    ////////////////////////////////////////////////////////////////////////////////
    void newLine() { printLog(""); }
    void newLineIf(bool bCondition) { if(bCondition) { printLog(""); } }
    void printSeparator();
    void printAligned(const String& s, char padding = PADDING, const String& wrapper = WRAPPER, UInt maxSize = 100);
    void printTextBox(const String& s);
    void printTextBox(const Vector<String>& strs);
    void printWarning(const String& s, UInt maxSize            = 100);
    void printWarningIndent(const String& s,  UInt indentLevel = 1, char trailing = ' ', UInt maxSize = 100);
    void printError(const String& s, UInt maxSize              = 100);
    void printErrorIndent(const String& s, UInt indentLevel    = 1, char trailing = ' ', UInt maxSize = 100);

    ////////////////////////////////////////////////////////////////////////////////
    template<class Function> void printRunTime(const char* caption, const Function& function)
    {
        printLog(Timer::getRunTime<Function>(caption, function));
    }

    template<class Function> void printRunTimeIndent(const char* caption, const Function& function, UInt indentLevel = 1, char trailing = ' ')
    {
        printLogIndent(Timer::getRunTime<Function>(caption, function), indentLevel, trailing);
    }

    template<class Function> void printRunTimeIf(const char* caption, const Function& function)
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIf(bResult, timer.getRunTime(caption));
    }

    template<class Function> void printRunTimeIndentIf(const char* caption, const Function& function, UInt indentLevel = 1, char trailing = ' ')
    {
        Timer timer;
        timer.tick();
        bool bResult = function();
        timer.tock();
        printLogIndentIf(bResult, timer.getRunTime(caption), indentLevel, trailing);
    }

    ////////////////////////////////////////////////////////////////////////////////
    void printLog(const String& s);
    void printLog(const String& s, spdlog::level::level_enum level);
    void printLogIndent(const String& s, UInt indentLevel = 1, char trailing = ' ');

    void printLogPadding(const String& s, UInt maxSize                                            = 100);
    void printLogPadding(const String& s, spdlog::level::level_enum level, UInt maxSize           = 100);
    void printLogPaddingIndent(const String& s, UInt indentLevel                                  = 1, char trailing = ' ', UInt maxSize = 100);
    void printLogPaddingIndent(const String& s, spdlog::level::level_enum level, UInt indentLevel = 1, char trailing = ' ', UInt maxSize = 100);

    void printLogIf(bool bCondition, const String& s);
    void printLogIf(bool bCondition, const String& s, spdlog::level::level_enum level);
    void printLogIndentIf(bool bCondition, const String& s, UInt indentLevel = 1, char trailing = ' ');

    ////////////////////////////////////////////////////////////////////////////////
    void flush();
    void printMemoryUsage();
    void printTotalRunTime();
    ////////////////////////////////////////////////////////////////////////////////
    void        cleanup(int signal = 0);
    static void signalHandler(int signal);

private:
    String getTotalRunTime();
    ////////////////////////////////////////////////////////////////////////////////
    static inline Vector<Logger*> s_Instances {};
    Clock::time_point             m_StartupTime {};
    SharedPtr<spdlog::logger>     m_ConsoleLogger = nullptr;
    SharedPtr<spdlog::logger>     m_FileLogger    = nullptr;
    bool                          m_bLog2Console;
    bool                          m_bLog2File;
};
