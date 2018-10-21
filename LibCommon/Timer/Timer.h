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

#include <chrono>
#include <string>
#include <cassert>
#include <sstream>
#include <iomanip>
#include <functional>

#include <LibCommon/Utils/Formatters.h>

//-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
class Timer
{
protected:
    using Clock = std::chrono::high_resolution_clock;

public:
    Timer()          = default;
    virtual ~Timer() = default;

    void tick()
    {
        assert(!m_TimerTicked);
        m_StartTime   = Clock::now();
        m_TimerTicked = true;
    }

    double tock()
    {
        assert(m_TimerTicked);
        m_EndTime     = Clock::now();
        m_TimerTicked = false;
        m_ElapsedTime = std::chrono::duration<double, std::milli>(m_EndTime - m_StartTime).count();

        return m_ElapsedTime;
    }

    String getRunTime()
    {
        if(m_TimerTicked) {
            tock();
        }
        m_Str.clear();
        m_Str += Formatters::toString(m_ElapsedTime);
        m_Str += String("ms");
        return m_Str;
    }

    String getRunTime(const String& caption)
    {
        if(m_TimerTicked) {
            tock();
        }

        m_Str.clear();
        m_Str += caption;
        m_Str += String(": ");
        m_Str += Formatters::toString(m_ElapsedTime);
        m_Str += String("ms");

        return m_Str;
    }

    template<class Function>
    static String getRunTime(const String& caption, const Function& function)
    {
        Timer timer;
        ////////////////////////////////////////////////////////////////////////////////
        timer.tick();
        function();
        timer.tock();
        ////////////////////////////////////////////////////////////////////////////////
        return timer.getRunTime(caption);
    }

private:
    Clock::time_point m_StartTime;
    Clock::time_point m_EndTime;
    String            m_Str;

    double m_ElapsedTime { 0.0 };
    bool   m_TimerTicked { false };
};
