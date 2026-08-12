/*!
 * @author  Thomas Richter <thomas.richter@ovgu.de>
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#include "include/VectorRotator.hpp"
#include "include/ParametricMesh.hpp"
#include "include/cgVector.hpp"
#include "include/constants.hpp"

namespace Nextsim {

/* An almost empty constructor with no rotation */
VectorRotator::VectorRotator(const std::vector<size_t>& dimsIn)
    : dims(dimsIn)
{
    det.assign(dims[0] * dims[1], 1.);
    ex.assign(det.size(), { 1, 0 });
    ey.assign(det.size(), { 0, 1 });
}

/* Constructor that uses two std::vectors of longitude and latitude to construct the unit vectors.
 * We construct a ParametricMesh object with the input coordinates at cell vertices. For most of the
 * domain, the unit vectors are calculated from the lower left corner, but the top row and last
 * column need to be handled separately.
 */
VectorRotator::VectorRotator(const std::vector<size_t>& dimsIn, const std::vector<FloatType>& lon,
    const std::vector<FloatType>& lat, const orientation orient)
    : dims(dimsIn)
{
    det.resize(dims[0] * dims[1]);
    ex.resize(det.size());
    ey.resize(det.size());

    switch (orient) {
    case orientation::EAST_NORTH:
        // Call the ENOrientation routine if we're in East-North orientation
        initENOrientation(lon, lat);
        break;
    case orientation::GRID: {
        // Build a smesh object for spherical coordinates
        ParametricMesh smesh(SPHERICAL);

        // Use the lon/lat coords as input for the ParametricMesh, convert, and rotate to Greenland
        smesh.coordinatesFromVectors(dims, lon, lat);
        smesh.TransformToRadians();
        smesh.RotatePoleToGreenland();

        /* Assemble the ex, ey, and det vectors needed by toParametricMesh and
         * fromParametricMesh. We start by constructing the element orientation everywhere except in
         * the last row and column by connecting the lower left corner of the grid cell with the
         * upper left and lower right.
         */

        // weights to connect lower left corner with its neighbours
        Eigen::Matrix<FloatType, 4, 1> ix({ -1, 1, 0, 0 });
        Eigen::Matrix<FloatType, 4, 1> iy({ -1, 0, 1, 0 });

        // Loop through the full smesh grid. This leaves the upper and right outer boundary.
#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid) {
            const Eigen::Matrix<FloatType, 4, 2> coe = smesh.coordinatesOfElement(eid);

            // The two vectors spanning the element give the direction of the ocean velocity.
            // NB! dimsIn != { smesh.nx, smesh.ny }
            const std::vector<size_t> ij = deIndexer({ smesh.nx, smesh.ny }, eid);
            const size_t k = indexer(dimsIn, ij);
            ex[k] = (coe.transpose() * ix).normalized();
            ey[k] = (coe.transpose() * iy).normalized();
            det[k] = ex[k](0) * ey[k](1) - ex[k](1) * ey[k](0);
        }

        /* Handle the edge cases by assuming a different connectivity within the smesh element */

        // Top row
        // weights to connect upper left corner with its neighbours.
        ix = { 0, 0, -1, 1 };
        iy = { -1, 0, 1, 0 };
#pragma omp parallel for
        for (size_t i = 0; i < smesh.nx; ++i) {
            const size_t j = smesh.ny - 1;

            const size_t eid = indexer({ smesh.nx, smesh.ny }, { i, j });
            const Eigen::Matrix<FloatType, 4, 2> coe = smesh.coordinatesOfElement(eid);

            // Place the results into i and j+1, because the reference is upper left corner
            const size_t k = indexer(dimsIn, { i, j + 1 });
            ex[k] = (coe.transpose() * ix).normalized();
            ey[k] = (coe.transpose() * iy).normalized();
            det[k] = ex[k](0) * ey[k](1) - ex[k](1) * ey[k](0);
        }

        //  Last column
        // weights to connect lower right corner with its neighbours.
        ix = { -1, 1, 0, 0 };
        iy = { 0, -1, 0, 1 };
#pragma omp parallel for
        for (size_t j = 0; j < smesh.ny; ++j) {
            const size_t i = smesh.nx - 1;

            const size_t eid = indexer({ smesh.nx, smesh.ny }, { i, j });
            const Eigen::Matrix<FloatType, 4, 2> coe = smesh.coordinatesOfElement(eid);

            // Place the results into i+1 and j, because the reference is lower right corner
            const size_t k = indexer(dimsIn, { i + 1, j });
            ex[k] = (coe.transpose() * ix).normalized();
            ey[k] = (coe.transpose() * iy).normalized();
            det[k] = ex[k](0) * ey[k](1) - ex[k](1) * ey[k](0);
        }

        // The remaining upper right corner
        // weights to connect upper right with its neighbours.
        ix = { 0, 0, -1, 1 };
        iy = { 0, -1, 0, 1 };

        // Upper right corner
        const size_t i = smesh.nx - 1;
        const size_t j = smesh.ny - 1;

        const size_t eid = indexer({ smesh.nx, smesh.ny }, { i, j });
        const Eigen::Matrix<FloatType, 4, 2> coe = smesh.coordinatesOfElement(eid);

        // Place the results into i+1 and j+1, because the reference is upper right corner
        const size_t k = indexer(dimsIn, { i + 1, j + 1 });
        ex[k] = (coe.transpose() * ix).normalized();
        ey[k] = (coe.transpose() * iy).normalized();
        det[k] = ex[k](0) * ey[k](1) - ex[k](1) * ey[k](0);

        break;
    }
    default:
        throw std::runtime_error(
            "VectorRotator::VectorRotator: Unknown orientation in constructor.\n");
    }
}

/* A constructor that uses a ModelArray with the model coordinates, and coordinates of cell vertices
 * to construct the unit vectors. Much simpler than the other one.
 */
VectorRotator::VectorRotator(const ModelArray& coords, const orientation orient)
{
    switch (orient) {
    case orientation::EAST_NORTH: {
        // Call the ENOrientation routine if we're in East-North orientation
        const auto lon = std::vector(
            coords.components(0).data(), coords.components(0).data() + coords.components(0).size());
        const auto lat = std::vector(
            coords.components(1).data(), coords.components(1).data() + coords.components(1).size());
        initENOrientation(lon, lat);
        break;
    }
    case orientation::GRID: {
        // Build a smesh object for spherical coordinates
        ParametricMesh smesh(SPHERICAL);

        dims = { ModelArray::size(ModelArray::Dimension::X),
            ModelArray::size(ModelArray::Dimension::Y) };

        // Build a ParametricMesh object and rotate to Greenland
        smesh.coordinatesFromModelArray(coords);
        smesh.RotatePoleToGreenland();

        /* Assemble the ex, ey, and det vectors needed by toParametricMesh and
         * fromParametricMesh. In this case, coords contains all the grid cell corners, making
         * things easy.
         */
        det.resize(dims[0] * dims[1]);
        ex.resize(det.size());
        ey.resize(det.size());

        // Connect the edge-midpoints to get the unit-vectors
        const Eigen::Matrix<FloatType, 4, 1> iix({ -0.5, 0.5, -0.5, 0.5 });
        const Eigen::Matrix<FloatType, 4, 1> iiy({ -0.5, -0.5, 0.5, 0.5 });

#pragma omp parallel for
        for (size_t eid = 0; eid < smesh.nelements; ++eid) {
            // construct the "direction" of the element, i.e. the ocean ex,ey-vectors
            const Eigen::Matrix<FloatType, 4, 2> coe = smesh.coordinatesOfElement(eid);

            ex[eid] = (coe.transpose() * iix).normalized();
            ey[eid] = (coe.transpose() * iiy).normalized();
            det[eid] = ex[eid](0) * ey[eid](1) - ex[eid](1) * ey[eid](0);
        }
        break;
    }
    default:
        throw std::runtime_error(
            "VectorRotator::VectorRotator: Unknown orientation in constructor.\n");
    }
}

/* We have a closed-form solution if the input/output vectors should be oriented along
 * east/north directions. It's just a rotation between the geographic and displaced poles with a
 * rather complicated angle:
 * \alpha = \atan2(\cos\phi_p \sin\Delta\lambda,
 *                                      \sin\phi_p \cos\phi - \cos\phi_p \sin\phi
 * \cos\Delta\lambda)
 */
void VectorRotator::initENOrientation(
    const std::vector<FloatType>& lon, const std::vector<FloatType>& lat)
{
    // TODO: The Greenland pole shouldn't be hardcoded!
    const FloatType polLon = radians(15.);
    const FloatType polLat = radians(40.);

#pragma omp parallel for
    for (size_t eid = 0; eid < det.size(); ++eid) {
        const std::vector<size_t> ij = deIndexer(dims, eid);
        const size_t i = ij[0];
        const size_t j = ij[1];

        const FloatType rLon = radians(lon[i]);
        const FloatType rLat = radians(lat[j]);

        // alpha = atan2(a, b)
        const FloatType deltaLon = polLon - rLon;
        const FloatType a = std::cos(polLat) * std::sin(deltaLon);
        const FloatType b = std::cos(rLat) * std::sin(polLat)
            - std::sin(rLat) * std::cos(polLat) * std::cos(deltaLon);

        // Instead of the atan2, cos, and sin functions
        const FloatType sinAlpha = a / std::hypot(a, b);
        const FloatType cosAlpha = b / std::hypot(a, b);

        ex[eid] = { cosAlpha, -sinAlpha };
        ey[eid] = { sinAlpha, cosAlpha };

        // det[eid] = ex[eid](0) * ey[eid](1) - ex[eid](1) * ey[eid](0);
        // The determinant is just one
        det[eid] = 1.;
    }
}

// A: From ocean to ParamMesh:
// ocean velocity is ox * ex + ey * ey. This can directly be evaluated:
void VectorRotator::toParametricMesh(const std::vector<FloatType>& uIn,
    const std::vector<FloatType>& vIn, std::vector<FloatType>& uOut,
    std::vector<FloatType>& vOut) const
{
#pragma omp parallel for
    for (size_t i = 0; i < uIn.size(); ++i) {
        const Eigen::Matrix<FloatType, 2, 1> pVelOut = ex[i] * uIn[i] + ey[i] * vIn[i];
        uOut[i] = pVelOut(0);
        vOut[i] = pVelOut(1);
    }
}

// B: from ParamMesh to ocean
// solve linear system such that ex * ox + ey * uy = v
void VectorRotator::fromParametricMesh(const std::vector<FloatType>& uIn,
    const std::vector<FloatType>& vIn, std::vector<FloatType>& uOut,
    std::vector<FloatType>& vOut) const
{
#pragma omp parallel for
    for (size_t i = 0; i < uIn.size(); ++i) {
        uOut[i] = -ey[i](0) / det[i] * vIn[i] + ey[i](1) / det[i] * uIn[i];
        vOut[i] = +ex[i](0) / det[i] * vIn[i] - ex[i](1) / det[i] * uIn[i];
    }
}

// A version of toParametricMesh which interpolates the output to CGVectors
template <int CG>
void VectorRotator::toParametricMesh(const std::vector<FloatType>& uIn,
    const std::vector<FloatType>& vIn, CGVector<CG>& uOut, CGVector<CG>& vOut) const
{
    uOut.setZero();
    vOut.setZero();

#pragma omp parallel for
    for (size_t i = 0; i < uIn.size(); ++i) {
        // compute center velocity in nextsim coordinate system
        const Eigen::Matrix<FloatType, 2, 1> Vcenter = ex[i] * uIn[i] + ey[i] * vIn[i];

        // average to CG-Vector
        // index of lower-left node in cell
        const std::vector<size_t> loc = deIndexer(dims, i);
        const size_t ix = loc[0];
        const size_t iy = loc[1];

        const size_t n0 = CG * (CG * dims[0] + 1) * iy + CG * ix;

        for (size_t cy = 0; cy <= CG; ++cy) {
            for (size_t cx = 0; cx <= CG; ++cx) {
                // weights for averaging from cell center to vertices
                constexpr double wgt[2][3] = { { 0.5, 0.5, 0 }, // CG1
                    { 0.5, 1.0, 0.5 } }; // CG2

                uOut(n0 + (CG * dims[0] + 1) * cy + cx)
                    += wgt[CG - 1][cx] * wgt[CG - 1][cy] * Vcenter(0);
                vOut(n0 + (CG * dims[0] + 1) * cy + cx)
                    += wgt[CG - 1][cx] * wgt[CG - 1][cy] * Vcenter(1);
            }
        }
    }
}

} // namespace Nextsim