/*!
 * @file BenchmarkCoordinates.cpp
 *
 * @date 26 Sept 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BenchmarkCoordinates.hpp"

namespace Nextsim {

bool BenchmarkCoordinates::isInitialized = false;

void BenchmarkCoordinates::setData()
{
    if (!isInitialized) {
        m_nx = ModelArray::dimensions(ModelArray::Type::H)[0];
        m_ny = ModelArray::dimensions(ModelArray::Type::H)[1];

        for (size_t j = 0; j < m_ny; ++j) {
            double yVal = j * dy;
            for (size_t i = 0; i < m_nx; ++i) {
                double xVal = i * dx;
                m_x(i, j) = xVal;
                m_y(i, j) = yVal;
            }
        }
        isInitialized = true;
    }
}
} /* namespace Nextsim */
