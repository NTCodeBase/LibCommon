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

#include <Utils/Logger.h>
#include <Utils/FileHelpers.h>
#include <Utils/Formatters.h>
#include <Utils/MemoryUsage.h>

#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <csignal>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
Logger::Logger(const String& loggerName, const String& rootPath, bool bLog2Console /*= true*/, bool bLog2File /*= false*/,
               spdlog::level::level_enum level /*= spdlog::level::level_enum::trace*/) :
    m_bLog2Console(bLog2Console), m_bLog2File(bLog2File)
{
    if(!m_bLog2Console && !m_bLog2File) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(m_bLog2Console) {
        m_ConsoleLogger = spdlog::stdout_color_mt(loggerName + String("[console_logger]"));
        m_ConsoleLogger->set_pattern("[%Y-%m-%d] [%H:%M:%S.%e] [%^%l%$] %v");
        m_ConsoleLogger->set_level(level);
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(m_bLog2File) {
        FileHelpers::createFolder(String(rootPath + "/Log"));
        auto prefix = rootPath + String("/Log/log_");
        UInt i      = 1;
        do {
            if(String file = prefix + std::to_string(i) + String(".txt");
               !FileHelpers::fileExisted(file)) {
                m_FileLogger = spdlog::create_async<spdlog::sinks::basic_file_sink_mt>(loggerName + String("[file_logger]"), file);
                m_FileLogger->set_pattern("[%Y-%m-%d] [%H:%M:%S.%e] [%^%l%$] %v");
                m_FileLogger->set_level(level);
                break;
            }
            ++i;
        } while(true);
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_StartupTime = Clock::now();
    ////////////////////////////////////////////////////////////////////////////////
    s_Instances.push_back(this);
    signal(SIGINT,   Logger::signalHandler);
    signal(SIGFPE,   Logger::signalHandler);
    signal(SIGSEGV,  Logger::signalHandler);
    signal(SIGTERM,  Logger::signalHandler);
    signal(SIGBREAK, Logger::signalHandler);
    signal(SIGABRT,  Logger::signalHandler);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printSeparator()
{
    printAligned("", '=');
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printAligned(const String& s, char padding, const String& wrapper, UInt maxSize /*= 100*/)
{
    size_t       length = s.length();
    const String str    = length == 0 ? s : String(" " + s + " ");
    length = str.length();

    size_t paddingSize = ((size_t)maxSize - length - 2 * wrapper.length()) / 2;
    size_t finalLength = wrapper.length() * 2 + length + paddingSize * 2;

    String finalStr;
    finalStr.reserve(finalLength + 1);

    finalStr += wrapper;
    finalStr += String(paddingSize, padding);
    finalStr += str;
    finalStr += (finalLength == static_cast<size_t>(maxSize)) ? String(paddingSize, padding) : String(paddingSize + 1, padding);
    finalStr += wrapper;

    printLog(finalStr);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTextBox(const String& s)
{
    printSeparator();
    printAligned("", ' ');
    printAligned(s);
    printAligned("", ' ');
    printSeparator();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTextBox(const Vector<String>& strs)
{
    printSeparator();
    printAligned("", ' ');
    for(const auto& s: strs) {
        printAligned(s);
    }
    printAligned("", ' ');
    printSeparator();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printWarning(const String& s, UInt maxSize)
{
    printLogPadding(s, spdlog::level::warn, maxSize);
}

void Logger::printWarningIndent(const String& s, UInt indentLevel /* = 1 */, char trailing /* = ' ' */, UInt maxSize /* = 100 */)
{
    printLogPaddingIndent(s, spdlog::level::warn, indentLevel, trailing, maxSize);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printError(const String& s, UInt maxSize)
{
    printLogPadding(s, spdlog::level::err, maxSize);
}

void Logger::printErrorIndent(const String& s, UInt indentLevel /* = 1 */, char trailing /* = ' ' */, UInt maxSize /* = 100 */)
{
    printLogPaddingIndent(s, spdlog::level::err, indentLevel, trailing, maxSize);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLog(const String& s)
{
    if(m_bLog2Console) {
        m_ConsoleLogger->info(s);
    }

    if(m_bLog2File) {
        m_FileLogger->info(s);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLog(const String& s, spdlog::level::level_enum level)
{
    if(m_bLog2Console) {
        m_ConsoleLogger->log(level, s);
    }

    if(m_bLog2File) {
        m_FileLogger->log(level, s);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogIndent(const String& s, UInt indentLevel /*= 1*/, char trailing /*= ' '*/)
{
    String str;
    str.reserve(256);
    str += String(INDENT_SIZE * indentLevel, trailing);
    str += s;

    printLog(str);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogPadding(const String& s, UInt maxSize /*= 100*/)
{
    auto str = s;
    str += String(" ");
    size_t paddingSize = (static_cast<size_t>(maxSize) - str.length());
    str += String(paddingSize, '*');

    if(m_bLog2Console) {
        m_ConsoleLogger->info(str);
    }

    if(m_bLog2File) {
        m_FileLogger->info(str);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogPadding(const String& s, spdlog::level::level_enum level, UInt maxSize /*= 100*/)
{
    auto str = s;
    str += String(" ");
    size_t paddingSize = (static_cast<size_t>(maxSize) - str.length());
    str += String(paddingSize, '*');

    if(m_bLog2Console) {
        m_ConsoleLogger->log(level, str);
    }

    if(m_bLog2File) {
        m_FileLogger->log(level, str);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogPaddingIndent(const String& s, UInt indentLevel /*= 1*/, char trailing /*= ' '*/, UInt maxSize /* = 100 */)
{
    String str;
    str.reserve(256);
    str += String(INDENT_SIZE * indentLevel, trailing);
    str += s;
    str += String(" ");

    size_t paddingSize = (static_cast<size_t>(maxSize) - str.length());
    str += String(paddingSize, '*');

    printLog(str);
}

void Logger::printLogPaddingIndent(const String& s, spdlog::level::level_enum level, UInt indentLevel /*= 1*/, char trailing /*= ' '*/, UInt maxSize /* = 100 */)
{
    String str;
    str.reserve(256);
    str += String(INDENT_SIZE * indentLevel, trailing);
    str += s;
    str += String(" ");

    size_t paddingSize = (static_cast<size_t>(maxSize) - str.length());
    str += String(paddingSize, '*');

    if(m_bLog2Console) {
        m_ConsoleLogger->log(level, str);
    }

    if(m_bLog2File) {
        m_FileLogger->log(level, str);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogIf(bool bCondition, const String& s)
{
    if(bCondition) {
        printLog(s);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogIf(bool bCondition, const String& s, spdlog::level::level_enum level)
{
    if(bCondition) {
        printLog(s, level);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogIndentIf(bool bCondition, const String& s, UInt indentLevel /*= 1*/, char trailing /*= ' '*/)
{
    if(bCondition) {
        printLogIndent(s, indentLevel, trailing);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::flush()
{
    if(m_bLog2Console) {
        m_ConsoleLogger->flush();
    }
    if(m_bLog2File) {
        m_FileLogger->flush();
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printMemoryUsage()
{
    String str;
    str.reserve(256);
    str += String("Memory usage: ");
    str += Formatters::toString(static_cast<double>(getCurrentRSS()) / 1048576.0);
    str += String(" MB(s). Peak: ");
    str += Formatters::toString(static_cast<double>(getPeakRSS()) / 1048576.0) + " MB(s).";

    printLog(str);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTotalRunTime()
{
    printLogPadding(getTotalRunTime(), 120);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::cleanup(int signal /*= 0*/)
{
    if(!m_bLog2Console && !m_bLog2File) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(signal != EXIT_SUCCESS && signal != SIGABRT) {
        newLine();
        switch(signal) {
            case SIGINT:
                printWarning("Interrupt signal caught (Ctrl_C pressed)");
                break;
            case SIGFPE:
                printWarning("Floating-point exception signal caught");
                break;
            case SIGSEGV:
                printWarning("Segmentation violation signal caught");
                break;
            case SIGTERM:
                printWarning("Termination signal caught");
                break;
            case SIGBREAK:
                printWarning("Ctrl-Break pressed");
                break;
            case SIGABRT:
                // dont handle this case: this is normal termination, signal raised by main window
                break;
            default:
                printWarning("Unknown signal caught: " + std::to_string(signal));
        }
        printMemoryUsage();
        printTotalRunTime();
    }
    flush();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::signalHandler(int signal)
{
    for(auto& logger: s_Instances) {
        logger->cleanup(signal);
    }
    spdlog::drop_all();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<typename Rep, typename Period>
void getDuration(std::chrono::duration<Rep, Period> t, UInt& n_days, UInt& n_hours, UInt& n_mins, UInt& n_secs)
{
    assert(0 <= t.count());

    // approximate because a day doesn't have a fixed length
    typedef std::chrono::duration<int, std::ratio<60* 60* 24>> days_t;

    auto days  = std::chrono::duration_cast<days_t>(t);
    auto hours = std::chrono::duration_cast<std::chrono::hours>(t - days);
    auto mins  = std::chrono::duration_cast<std::chrono::minutes>(t - days - hours);
    auto secs  = std::chrono::duration_cast<std::chrono::seconds>(t - days - hours - mins);

    n_days  = static_cast<UInt>(days.count());
    n_hours = static_cast<UInt>(hours.count());
    n_mins  = static_cast<UInt>(mins.count());
    n_secs  = static_cast<UInt>(secs.count());
}

String Logger::getTotalRunTime()
{
    UInt days, hours, mins, secs;
    getDuration(Clock::now() - m_StartupTime, days, hours, mins, secs);

    std::stringstream strBuilder;
    strBuilder.str("");
    strBuilder << "Total time: "
               << days << "(days), "
               << hours << "(hours), "
               << mins << "(mins), "
               << secs << "(secs).";

    return strBuilder.str();
}
