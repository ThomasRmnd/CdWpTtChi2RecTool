#ifndef CDWPTTCHI2RECTOOL_UTILS_TIMER_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_TIMER_HPP_

#include <chrono>
#include <iostream>
#include <string>
#include <vector>

#include "SniperKernel/SniperLog.h"

class timer {

public:

    timer(const std::string& msg) : m_msg(msg), m_stopped(true) {};

    ~timer() {
        if (m_durations.empty()) return;
        std::chrono::duration<double> total = m_durations[0];
        // std::ostringstream oss;
        // oss << "Detailed time consumption ==> 0: " << m_durations[0].count() << "s, "; 
        for (std::size_t i = 1; i < m_durations.size() - 1; ++i) {
            // oss << i << ": " << m_durations[i].count() << "s, ";
            total += m_durations[i];
        }
        // oss << m_durations.size() - 1 << ": " << m_durations.back().count() << "s\n";
        // LogInfo << oss.str();
        total += m_durations.back();
        LogInfo << "Total time consumption: " << total.count() << "s ==> Mean: " << total.count() / m_durations.size() << "s\n";
    };

    void start() {
        m_start = std::chrono::high_resolution_clock::now();
        m_stopped = false;
    };

    void stop() {
        if (!m_stopped) {
            m_durations.push_back(std::chrono::high_resolution_clock::now() - m_start);
            m_stopped = true;
        }
    }

private:

    std::string m_msg;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::vector<std::chrono::duration<double>> m_durations;
    bool m_stopped;

};

class timer_guard {

public:

    timer_guard(timer& t) : m_timer(t) {
        m_timer.start();
    };

    ~timer_guard() {
        m_timer.stop();
    };

private:

    timer& m_timer;

};

#endif // CDWPTTCHI2RECTOOL_UTILS_TIMER_HPP_