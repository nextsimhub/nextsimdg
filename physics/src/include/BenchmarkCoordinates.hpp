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
    BenchmarkCoordinates() = 0;
    ~BenchmarkCoordinates() = 0;

    static const double dx = 25000.;
    static const double dy = 25000.;

    //! Sets the data for the benchmark after the initial restart file has been read
    static void setData();

    //! Returns the x dimension of the benchmark coordinate arrays
    static size_t nx() const { return m_nx; }
    //! Returns the y dimension of the benchmark coordinate arrays
    static size_t ny() const { return m_ny; }

    /*!
     * Returns the benchmark x coordinate array
     */
    static const ModelArray& x() const { return m_x; }
    /*!
     * Returns the benchmark y coordinate array
     */
    static const ModelArray& y() const { return m_y; }

private:
    static bool isInitialized = false;

    // Size of the array dimensions
    static size_t m_nx;
    static size_t m_ny;

    // x and y coordinate arrays
    static ModelArray m_x;
    static ModelArray m_y;
};

} /* namespace Nextsim */

#endif /* BENCHMARKCOORDINATES_HPP */
