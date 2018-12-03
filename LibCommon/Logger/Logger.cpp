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

#include <LibCommon/Logger/Logger.h>
#include <LibCommon/Utils/FileHelpers.h>
#include <LibCommon/Utils/Formatters.h>
#include <LibCommon/Utils/MemoryUsage.h>

#include <spdlog/async.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <csignal>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
SharedPtr<Logger> Logger::createLogger(const String& loggerName, const String& rootPath, bool bLog2Console /*= true*/, bool bLog2File /*= false*/,
                                       LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    auto logger = std::make_shared<MakeSharedEnabler>(loggerName, rootPath, bLog2Console, bLog2File, consoleLogLevel, fileLogLevel);
    s_Instances.push_back(logger);
    return logger;
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::removeLogger(const SharedPtr<Logger>& logger) {
    __NT_REQUIRE(std::find(std::begin(s_Instances), std::end(s_Instances), logger) != s_Instances.end());
    s_Instances.remove(logger);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
Logger::Logger(const String& loggerName, const String& rootPath, bool bLog2Console /*= true*/, bool bLog2File /*= false*/,
               LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) :
    m_bLog2Console(bLog2Console), m_bLog2File(bLog2File) {
    if(!m_bLog2Console && !m_bLog2File) {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////
    std::scoped_lock lock{ s_InstancingMutex };
    __NT_UNUSED(lock);
    m_InstanceIdx = s_NumCreatedInstances++;
    ////////////////////////////////////////////////////////////////////////////////
    if(m_bLog2Console) {
        m_ConsoleLogger = spdlog::stdout_color_mt(loggerName + String("[console_logger][idx=") + std::to_string(m_InstanceIdx) + String("]"));
        m_ConsoleLogger->set_pattern("[%Y-%m-%d] [%H:%M:%S.%e] [%^%l%$] %v");
        m_ConsoleLogger->set_level(consoleLogLevel);
    }
    ////////////////////////////////////////////////////////////////////////////////
    if(m_bLog2File) {
        FileHelpers::createFolder(String(rootPath + "/Log"));
        auto prefix = rootPath + String("/Log/log_");
        UInt i      = 1;
        do {
            if(String file = prefix + std::to_string(i) + String(".txt");
               !FileHelpers::fileExisted(file)) {
                m_FileLogger = spdlog::create_async<spdlog::sinks::basic_file_sink_mt>(loggerName + String("[file_logger][idx=") +
                                                                                       std::to_string(m_InstanceIdx) + String("]"), file);
                m_FileLogger->set_pattern("[%Y-%m-%d] [%H:%M:%S.%e] [%^%l%$] %v");
                m_FileLogger->set_level(fileLogLevel);
                break;
            }
            ++i;
        } while(true);
    }
    ////////////////////////////////////////////////////////////////////////////////
    m_StartupTime = Clock::now();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
Logger::~Logger() {
    if(m_bLog2Console) {
        spdlog::drop(m_ConsoleLogger->name());
    }
    if(m_bLog2File) {
        spdlog::drop(m_FileLogger->name());
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogIndent(const String& s, UInt indentLevel /*= 1*/,
                            LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    String s_formatted;
    s_formatted.reserve(256);
    s_formatted += String(s_IndentSize * indentLevel, s_PrefixPadding);
    s_formatted += s;
    ////////////////////////////////////////////////////////////////////////////////
    printLog(s_formatted, consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogPadding(const String& s, LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    auto s_formatted = s;
    s_formatted.reserve(256);
    s_formatted += String(" ");
    s_formatted += String(s_PaddingMaxSize - s_formatted.length(), s_SuffixPadding);
    printLog(s_formatted, consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printLogPaddingIndent(const String& s, UInt indentLevel /*= 1*/, LogLevel consoleLogLevel /*= LogLevel::trace*/,
                                   LogLevel fileLogLevel /*= LogLevel::trace*/) {
    String s_formatted;
    s_formatted.reserve(256);
    s_formatted += String(s_IndentSize * indentLevel, s_PrefixPadding);
    s_formatted += s;
    s_formatted += String(" ");
    s_formatted += String(s_PaddingMaxSize - s_formatted.length(), s_SuffixPadding);
    printLog(s_formatted, consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printCenterAligned(const String& s, char paddingChar /*= Logger::s_PrefixPadding*/,
                                LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    size_t       length = s.length();
    const String str    = length == 0 ? s : String(" " + s + " ");
    length = str.length();

    size_t paddingSize = (s_PaddingMaxSize - length - 2 * s_Wrapper.length()) / 2;
    size_t finalLength = s_Wrapper.length() * 2 + length + paddingSize * 2;

    String s_formatted;
    s_formatted.reserve(finalLength + 1);
    s_formatted += s_Wrapper;
    s_formatted += String(paddingSize, paddingChar);
    s_formatted += str;
    s_formatted += (finalLength == s_PaddingMaxSize) ? String(paddingSize, paddingChar) : String(paddingSize + 1, paddingChar);
    s_formatted += s_Wrapper;
    printLog(s_formatted, consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTextBox(const String& s, LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    separatorLine(consoleLogLevel, fileLogLevel);
    printCenterAligned("", ' ', consoleLogLevel, fileLogLevel);
    printCenterAligned(s,  ' ', consoleLogLevel, fileLogLevel);
    printCenterAligned("", ' ', consoleLogLevel, fileLogLevel);
    separatorLine(consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTextBox(const StdVT<String>& strs, LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    separatorLine(consoleLogLevel, fileLogLevel);
    printCenterAligned("", ' ', consoleLogLevel, fileLogLevel);
    for(const auto& s: strs) {
        printCenterAligned(s, ' ', consoleLogLevel, fileLogLevel);
    }
    printCenterAligned("", ' ', consoleLogLevel, fileLogLevel);
    separatorLine(consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::separatorLine(int stype /*= 0*/, LogLevel consoleLogLevel /*= LogLevel::info*/, LogLevel fileLogLevel /*= LogLevel::info*/) {
    if(stype == 0) {
        printCenterAligned("", '=', consoleLogLevel, fileLogLevel);
    } else {
        printLog(String(s_PaddingMaxSize, s_SuffixPadding), consoleLogLevel, fileLogLevel);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::flush() {
    if(m_bLog2Console) {
        m_ConsoleLogger->flush();
    }
    if(m_bLog2File) {
        m_FileLogger->flush();
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printMemoryUsage(LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    String str;
    str.reserve(256);
    str += String("Memory usage: ");
    str += Formatters::toString(static_cast<double>(getCurrentRSS()) / 1048576.0);
    str += String(" MB(s). Peak: ");
    str += Formatters::toString(static_cast<double>(getPeakRSS()) / 1048576.0) + " MB(s).";
    printLog(str, consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::printTotalRunTime(LogLevel consoleLogLevel /*= LogLevel::trace*/, LogLevel fileLogLevel /*= LogLevel::trace*/) {
    printLogPadding(getTotalRunTime(), consoleLogLevel, fileLogLevel);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::flushAll(int signal) {
    std::list<SharedPtr<Logger>> deletedLoggers;
    for(auto& logger: s_Instances) {
        if(logger.get()) {
            if(signal > 0) {
                logger->newLine();
                switch(signal) {
                    case SIGINT:
                        logger->printWarning("Interrupt signal caught (Ctrl_C pressed)");
                        break;
                    case SIGSEGV:
                        logger->printWarning("Segmentation violation signal caught");
                        break;
                    case SIGTERM:
                        logger->printWarning("Termination signal caught");
                        break;
                    case SIGABRT:
                        logger->printWarning("Abortion signal caught");
                        break;
                    default:
                        logger->printWarning("Unknown signal caught: " + std::to_string(signal));
                }
                logger->printMemoryUsage();
                logger->printTotalRunTime();
            } // end signal > 0
            logger->flush();
        } else {
            deletedLoggers.push_back(logger);
        }
    }
    for(auto& logger: deletedLoggers) {
        s_Instances.remove(logger);
    }
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::signalHandler(int signal) {
    flushAll(signal);
    shutdown();
    exit(signal);
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
void Logger::shutdown() {
    spdlog::shutdown();
}

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
template<typename Rep, typename Period>
void getDuration(std::chrono::duration<Rep, Period> t, UInt& n_days, UInt& n_hours, UInt& n_mins, UInt& n_secs) {
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

String Logger::getTotalRunTime() {
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
