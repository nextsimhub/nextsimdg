/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#ifndef VECTORROTATOR_HPP
#define VECTORROTATOR_HPP

#include "FloatType.hpp"
#include "ModelArray.hpp"
#include "cgVector.hpp"

#include <Eigen/src/Core/Matrix.h>
#include <vector>

namespace Nextsim {

/*!
 * A class to perform vector rotations for inputs and outputs between the nextsim, displaced-pole
 * coordinate system and that of an arbitrary input or output grid. The constructor calculates the
 * unit vectors of the input/output grid in the nextsim coordinate system. These are then used to
 * transform input vectors in toParametricMesh and fromParametricMesh.
 */
class VectorRotator {
public:
    VectorRotator() = delete;

    /*!
     * @brief Constructor for the VectorRotator class
     *
     * @param dimsIn The dimensions of the input grid
     * @param lon The longitudes of the grid points
     * @param lat The latitudes of the grid points
     * @param isSpherical Whether the grid is spherical (true, by default)
     */
    VectorRotator(const std::vector<size_t>& dimsIn, const std::vector<FloatType>& lon,
        const std::vector<FloatType>& lat, bool isSpherical = true);

    /*!
     * @brief Constructor for the VectorRotator class
     *
     * @param coords A ModelArray containing the coordinates of the grid points
     * @param isSpherical Whether the grid is spherical (true, by default)
     */
    explicit VectorRotator(const ModelArray& coords, bool isSpherical = true);

    /*!
     * @brief Transforms velocities to the parametric mesh. All vectors are at the grid cell centre.
     *
     * @param uIn The input u-velocities
     * @param vIn The input v-velocities
     * @param uOut The output u-velocities
     * @param vOut The output v-velocities
     */
    void toParametricMesh(const std::vector<FloatType>& uIn, const std::vector<FloatType>& vIn,
        std::vector<FloatType>& uOut, std::vector<FloatType>& vOut) const;

    /*!
     * @brief Transforms velocities from the parametric mesh. All vectors are at the grid cell
     * centre.
     *
     * @param uIn The input u-velocities
     * @param vIn The input v-velocities
     * @param uOut The output u-velocities
     * @param vOut The output v-velocities
     */
    void fromParametricMesh(const std::vector<FloatType>& uIn, const std::vector<FloatType>& vIn,
        std::vector<FloatType>& uOut, std::vector<FloatType>& vOut) const;

    /*!
     * @brief Transforms velocities to the parametric mesh. Input vectors are at the grid cell
     * centre, the outputs are CGVectors.
     *
     * @param uIn The input u-velocities
     * @param vIn The input v-velocities
     * @param uOut The output u-velocities
     * @param vOut The output v-velocities
     */
    template <int CG>
    void toParametricMesh(const std::vector<FloatType>& uIn, const std::vector<FloatType>& vIn,
        CGVector<CG>& uOut, CGVector<CG>& vOut) const;

private:
    std::vector<Eigen::Matrix<FloatType, 2, 1>> ex, ey;
    std::vector<FloatType> det;
    std::vector<size_t> dims;
};

} // namespace Nextsim

#endif // VECTORROTATOR_HPP
