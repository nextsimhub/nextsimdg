#ifndef KOKKOSTIMER_HPP
#define KOKKOSTIMER_HPP

namespace Nextsim {

constexpr bool DETAILED_MEASUREMENTS = true;

template <bool Active> class KokkosTimer {
public:
    KokkosTimer(const std::string& _name)
        : m_name(_name)
        , m_count(0)
        , m_total(0.0)
    {
    }

    ~KokkosTimer()
    {
        if constexpr (Active) {
            std::cout << m_name << " " << m_total << " " << m_total / m_count << " " << m_count
                      << std::endl;
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
            Kokkos::fence();
            auto end = std::chrono::high_resolution_clock::now();
            ++m_count;
            m_total += std::chrono::duration<double>(end - m_start).count();
        }
    }
    void print()
    {
        if constexpr (Active) {
            std::cout << m_name << " total: " << m_total << ", avg: " << m_total / m_count
                      << ", #: " << m_count << "\n";
        }
    }

private:
    std::string m_name;
    int m_count;
    double m_total;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
};

}

#endif