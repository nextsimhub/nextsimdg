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
    /*!
     * @brief An eunumerator class for vector orientation
     *
     * @enum GRID: Vectors are aligned with the grid.
     * @enum EAST_NORTH: Vectors are aligned with the east and north directions.
     */
    enum class orientation { GRID, EAST_NORTH };

    VectorRotator() = delete;

    /*!
     * @brief Constructor for the VectorRotator class giving no rotation
     *
     * @param dimsIn The dimensions of the input grid
     */
    explicit VectorRotator(const std::vector<size_t>& dimsIn);

    /*!
     * @brief Constructor for the VectorRotator class based on vectors of lon/lat
     *
     * @param dimsIn The dimensions of the input grid
     * @param lon The longitudes of the grid points
     * @param lat The latitudes of the grid points
     * @param orient Vector orientation (GRID or EAST_NORTH)
     */
    VectorRotator(const std::vector<size_t>& dimsIn, const std::vector<FloatType>& lon,
        const std::vector<FloatType>& lat, const orientation orient);

    /*!
     * @brief Constructor for the VectorRotator class based on a the coordinates in a ModelArray
     *
     * @param coords A ModelArray containing the coordinates of the grid points
     * @param orient Vector orientation (GRID or EAST_NORTH)
     */
    explicit VectorRotator(const ModelArray& coords, const orientation orient);

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
    /*!
     * @brief Initialises the unit vectors using coordinate rotation from the geographic pole to the
     * displaced one over Greenland.
     *
     * @param lon The longitudes of the grid points
     * @param lat The latitudes of the grid points
     */
    void initENOrientation(const std::vector<FloatType>& lon, const std::vector<FloatType>& lat);

    std::vector<Eigen::Matrix<FloatType, 2, 1>> ex, ey;
    std::vector<FloatType> det;
    std::vector<size_t> dims;
};

} // namespace Nextsim

#endif // VECTORROTATOR_HPP
