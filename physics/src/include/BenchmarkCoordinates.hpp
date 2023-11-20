/*!
 * @file BenchmarkCoordinates.hpp
 *
 * @date 26 Sept 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef BENCHMARKCOORDINATES_HPP
#define BENCHMARKCOORDINATES_HPP

#include "include/ModelArray.hpp"

namespace Nextsim {

//! A class to calculate and store the coordinates used in the BenchmarkOcean and
//! BenchmarkAtmosphere classes
class BenchmarkCoordinates {
public:
    static double dx;
    static double dy;

    //! Sets the data for the benchmark after the initial restart file has been read
    static void setData();

    //! Returns the x dimension of the benchmark coordinate arrays
    static size_t nx() { return m_nx; }
    //! Returns the y dimension of the benchmark coordinate arrays
    static size_t ny() { return m_ny; }

    /*!
     * Returns the benchmark x coordinate array
     */
    static const ModelArray& x() { return m_x; }
    /*!
     * Returns the benchmark y coordinate array
     */
    static const ModelArray& y() { return m_y; }

    /*!
     * Returns the benchmark fractional x coordinate array x/(nx dx)
     */
    static const ModelArray fx() { return m_x / (m_nx * dx); }
    /*!
     * Returns the benchmark fractional y coordinate array y/(ny dy)
     */
    static const ModelArray fy() { return m_y / (m_ny * dy); }

private:
    BenchmarkCoordinates() = default;
    ~BenchmarkCoordinates() = default;

    static bool isInitialized;

    // Size of the array dimensions
    static size_t m_nx;
    static size_t m_ny;

    // x and y coordinate arrays
    static HField m_x;
    static HField m_y;
};

} /* namespace Nextsim */

#endif /* BENCHMARKCOORDINATES_HPP */
