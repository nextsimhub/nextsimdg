

#ifndef __stopwatch_H
#define __stopwatch_H

/**
 * This class will ne Removed after merge with physis
 **/

#include <cassert>
#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

/*----------------------------------------------------*/

namespace Nextsim {
/**
 *
 * StopWatch that measure the 'real time'
 *
 **/
class StopWatch {
protected:
    std::chrono::high_resolution_clock::time_point last_time;
    std::chrono::duration<double> sum_time;

    std::clock_t last_cpu_time;
    double sum_cpu_time;

    size_t activations;

    bool running;

public:
    StopWatch()
        : sum_time(std::chrono::duration<double>::zero())
        , sum_cpu_time(0)
        , activations(0)
        , running(false)
    {
    }

    virtual void reset()
    {
        sum_time = std::chrono::duration<double>::zero();
        sum_cpu_time = 0;
        activations = 0;
        running = false;
    }

    virtual void start()
    {
        if (!running) {
            activations++;
            last_time = std::chrono::high_resolution_clock::now();
            last_cpu_time = std::clock();
            running = true;
        }
    }

    virtual void stop()
    {
        std::chrono::high_resolution_clock::time_point now
            = std::chrono::high_resolution_clock::now();
        sum_time += std::chrono::duration_cast<std::chrono::duration<double>>(now - last_time);
        last_time = std::chrono::high_resolution_clock::now();

        std::clock_t now_cpu = std::clock();
        sum_cpu_time += static_cast<double>(now_cpu - last_cpu_time) / CLOCKS_PER_SEC;
        last_cpu_time = std::clock();

        running = false;
    }

    virtual double read() const
    {
        assert(!running);
        return sum_time.count();
    }

    double readPerActivations() const
    {
        assert(!running);
        return sum_time.count() / activations;
    }

    virtual double readCPU() const
    {
        assert(!running);
        return sum_cpu_time;
    }

    virtual double read100() const
    {
        return static_cast<double>(static_cast<int>(read() * 100)) / 100.0;
    }

    double readCPU100() const
    {
        return static_cast<double>(static_cast<int>(read() * 100)) / 100.0;
    }

    size_t readActivations() const { return activations; }

    bool isRunning() const { return running; }
};

/*----------------------------------------------------*/

/**
 *
 * class for managing StopWatches
 *
 **/

class Timer {
protected:
    std::map<std::string, StopWatch> watches;
    std::map<std::string, size_t> counters;

public:
    void reset()
    {
        watches.clear();
        counters.clear();
    }
    void reset(const std::string& label)
    {
        auto it = watches.find(label);
        assert(it != watches.end());
        it->second.reset();

        auto counter = counters.find(label);
        assert(counter != counters.end());
        counters.erase(counter);
    }

    void start(const std::string& label) { watches[label].start(); }

    void stop(const std::string& label)
    {
        assert(watches.find(label) != watches.end());
        watches[label].stop();
    }

    void count(const std::string& label) { counters[label]++; }

    double read(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());
        return it->second.read();
    }

    double read100(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());
        return it->second.read100();
    }

    double cpuread(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());
        return it->second.readCPU();
    }
    double cpuread100(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());
        return it->second.readCPU100();
    }

    size_t counterread(const std::string& label) const
    {
        auto it = counters.find(label);
        assert(it != counters.end());
        return it->second;
    }

    void print(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());

        std::cout << std::setw(24) << std::left << label << std::setw(12) << std::right
                  << it->second.read() << std::setw(12) << std::right << it->second.readCPU()
                  << std::setw(12) << std::right << it->second.readPerActivations() << std::setw(12)
                  << std::right << it->second.readActivations() << std::endl;
    }

    void print100(const std::string& label) const
    {
        auto it = watches.find(label);
        assert(it != watches.end());

        std::cout << label << "\t" << it->second.read100() << "\t" << it->second.readCPU100()
                  << std::endl;
    }

    void print100tofile(const std::string& filename) const
    {
        std::ofstream watch_logfile(filename);
        watch_logfile.precision(12);
        for (const auto& it : watches) {
            watch_logfile << it.first << "\t" << it.second.read100() << "\t"
                          << it.second.readCPU100() << std::endl;
        }
        watch_logfile.close();
    }

    void print() const
    {
        std::cout << std::showpoint << std::fixed << std::setprecision(6);
        for (auto it : watches) {
            print(it.first);
        }
        for (auto it : counters) {
            std::cerr << std::setw(24) << std::left << it.first << std::setw(12) << std::right
                      << it.second << std::endl;
        }
        std::cout << std::endl
                  << std::noshowpoint;
    }

    void print100() const
    {
        for (auto it : watches) {
            print100(it.first);
        }
        std::cout << std::endl;
    }
};

template <bool Active>
class PerfTimerAlt {
public:
    PerfTimerAlt(const std::string& _name)
        : m_name(_name)
        , m_count(0)
        , m_total(0.0)
    {
    }

    ~PerfTimerAlt()
    {
        if constexpr (Active) {
            std::cout << m_name << " " << m_total << " " << m_total / m_count << " " << m_count << std::endl;
        }
    }

    void start()
    {
        if constexpr (Active) {
            m_start = std::chrono::high_resolution_clock::now();
        }
    }
    void stop()
    {
        if constexpr (Active) {
            // ensure that gpu tasks are finished
            auto end = std::chrono::high_resolution_clock::now();
            ++m_count;
            m_total += std::chrono::duration<double>(end - m_start).count();
        }
    }
    void print()
    {
        if constexpr (Active) {
            std::cout << m_name << " total: " << m_total << ", avg: " << m_total / m_count << "\n";
        }
    }

private:
    std::string m_name;
    int m_count;
    double m_total;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

} // namespace Nextsim

#endif /* __stopwatch_H */
