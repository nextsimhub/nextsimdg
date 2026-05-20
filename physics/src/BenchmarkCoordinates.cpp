/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 */

#include "include/BenchmarkCoordinates.hpp"

namespace Nextsim {

bool BenchmarkCoordinates::isInitialized = false;
FloatType BenchmarkCoordinates::dx = 25000.;
FloatType BenchmarkCoordinates::dy = 25000.;
size_t BenchmarkCoordinates::m_nx;
size_t BenchmarkCoordinates::m_ny;
HField BenchmarkCoordinates::m_x(ModelArray::Type::H);
HField BenchmarkCoordinates::m_y(ModelArray::Type::H);

void BenchmarkCoordinates::setData()
{
    if (!isInitialized) {
        m_nx = ModelArray::dimensions(ModelArray::Type::H)[0];
        m_ny = ModelArray::dimensions(ModelArray::Type::H)[1];

        m_x.reinitialize();
        m_y.reinitialize();

        dx = 512e3 / m_nx;
        dy = 512e3 / m_ny;

        for (size_t j = 0; j < m_ny; ++j) {
            FloatType yVal = j * dy;
            for (size_t i = 0; i < m_nx; ++i) {
                FloatType xVal = i * dx;
                m_x(i, j) = xVal;
                m_y(i, j) = yVal;
            }
        }
        isInitialized = true;
    }
}
} /* namespace Nextsim */
