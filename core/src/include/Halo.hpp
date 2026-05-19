/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Halo exchange class
 * @details
 *
 * Halo exchange class
 *
 * All functionality for halo exchange between MPI ranks is contained in this class.
 *
 * Halo supports the main data structures of NextSim e.g., ModelArray, DGVector, DGVectorHolder and
 * CGVector.
 *
 * The halos are exchange via one-sided MPI communication using Remote Memory Access (RMA).
 */

#ifndef HALO_HPP
#define HALO_HPP

#include <cstddef>
#include <numeric>
#include <vector>

#include "Slice.hpp"
#include "cgVector.hpp"
#include "dgVector.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelMPI.hpp"
#include "include/ModelMetadata.hpp"
#include "include/dgVectorHolder.hpp"
#include "mpi.h"

#ifndef HALOWIDTH
#define HALOWIDTH 1
#endif

namespace Nextsim {

/*!
 * @brief A class to facilitate halo exchange between MPI ranks
 */
class Halo {
public:
    /*!
     * @brief Constructs a halo object from ModelArray
     * @param ma ModelArray object to create halo from
     */
    Halo(ModelArray& ma)
    {
        m_numComps = ma.nComponents();
        isVertex = ma.getType() == ModelArray::Type::VERTEX;
        setSpatialDims();
        initializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from DGVectorHolder
     * @param dgvh DGVectorHolder object to create halo from
     */
    template <int N>
    Halo(DGVectorHolder<N>& dgvh)
        : Halo(static_cast<DGVector<N>&>(dgvh))
    {
    }

    /*!
     * @brief Constructs a halo object from DGVector
     * @param dgv DGVector object to create halo from
     */
    template <int N> Halo(DGVector<N>& dgv)
    {
        m_numComps = N;
        setSpatialDims();
        initializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from CGVector
     * @param cgv CGVector object to create halo from
     */
    template <int N> Halo(CGVector<N>& cgv)
    {
        m_numComps = 1;
        isCG = true;
        CGdegree = N;
        nCells = CGdegree;
        setSpatialDims();
        initializeHaloMetadata();
    }

    static const size_t haloWidth = HALOWIDTH; // how many cells wide is the halo region

private:
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;
    using Edge = ModelMetadata::Edge;
    using Corner = ModelMetadata::Corner;
    using VBounds = ArraySlicer::Slice::VBounds;

    const typedef ArraySlicer::SliceIter::MultiDim MultiDim;
    const MultiDim m_haloDims;

    size_t m_Nx; // local extent including halo cells in x-direction
    size_t m_Ny; // local extent including halo cells in y-direction
    size_t m_innerNx; // local extent in x-direction
    size_t m_innerNy; // local extent in y-direction
    size_t m_numHaloCells; // number of halo cells
    size_t m_numComps; // number of DG components
    int nCells = 1; // this is the number of points to grab per cell, per direction (1 for
                    // everything except CG Fields)

    std::array<size_t, Edge::N_EDGE> m_edgeLengths; // array containing length of each edge
    std::array<Edge, Edge::N_EDGE> edges = ModelMetadata::edges; // array of edge enums
    std::array<Corner, Corner::N_CORNER> corners = ModelMetadata::corners; // array of edge enums

    /**
     * @brief Sets the spatial dimensions of the domain
     *
     * This function sets the inner dimensions of the computational domain based on metadata,
     * adjusts dimensions for vertex fields if necessary, and calculates total dimensions
     * including halo cells.
     */
    void setSpatialDims();

    /**
     * @brief Initializes halo metadata including buffers and slice information
     *
     * Sets up the send/receive buffers for each component and defines the boundary
     * slices for halo exchanges. Handles both vertex and non-vertex cases.
     *
     * - Calculates total number of halo cells based on halo width
     * - Allocates send/receive buffers for each component
     * - Defines edge lengths
     * - Sets up inner and outer slices for Bottom, Right, Top, Left edges
     */
    void initializeHaloMetadata();

    /**
     * @brief Return true if the provided edge is vertical (LEFT or RIGHT).
     *
     * @param edge Edge enum to test
     * @return true if edge is LEFT or RIGHT, false otherwise
     */
    static inline bool isVertical(const Edge edge)
    {
        return (edge == Edge::LEFT) || (edge == Edge::RIGHT);
    }

    /**
     * @brief Return true if the provided edge is horizontal (TOP or BOTTOM).
     *
     * @param edge Edge enum to test
     * @return true if edge is TOP or BOTTOM, false otherwise
     */
    static inline bool isHorizontal(const Edge edge)
    {
        return (edge == Edge::TOP) || (edge == Edge::BOTTOM);
    }

    MPI_Win m_win; // RMA memory window object (used for sharing send buffers between ranks)

    bool isVertex = false; // some ModelArrays can be of type VERTEX
    bool isCG = false; // is layout CGVector
    size_t CGdegree = 0;
    size_t sendBufferSize = 0;
    size_t recvBufferSize = 0;

    std::vector<std::vector<double>>
        send; // buffer to store halo region that will be read by other ranks
    std::vector<std::vector<double>>
        recv; // buffer to store halo region which is read from other ranks

    /*!
     * @brief Open memory window to exchange send buffer between MPI ranks.
     *
     * @details this is not intended to be used manually. It should only be called as part of the
     * update method.
     *
     * @param idx Component index for which to open the memory window
     */
    void openMemoryWindow(size_t idx);

    /*!
     * @brief Initialise memory window to exchange send buffer between MPI ranks.
     *
     * @details this is not intended to be used manually. It should only be called as part of the
     * update method.
     */
    void closeMemoryWindow();

    /*!
     * @brief Populate send buffer with halo data of the specified ModelArray
     */
    template <typename S> void populateSendBuffers(S& source)
    {
        size_t buffer_len = sendBufferSize;
        if (isCG) {
            buffer_len = sendBufferSize / CGdegree;
        }
        for (size_t comp = 0; comp < m_numComps; ++comp) {

            // sourceToSendBuffer requires the `buffer_map` and `source_map` to have a similar
            // memory layout to take advantage of Eigen notation for copying the data, and they both
            // need to be eigen array-like objects. `Eigen::Map` allow us to represent a pointer to
            // an array of data as an Eigen::ArrayXXd with specified dimensions.
            //
            // `source_map` maps the raw data pointer from `source` to a 2D Eigen array with
            // dimensions m_Nx x m_Ny Uses dynamic inner stride of m_numComps to handle column-major
            // storage format
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<Eigen::Dynamic>> source_map(
                source.col(comp).data(), m_Nx, m_Ny, Eigen::InnerStride<>(m_numComps));

            // The `send` buffers are stored as a vector of vectors. There is a send buffer for each
            // component (e.g., DG, CG etc.). Each component will be of length `sendBufferSize`.
            //
            // In most cases `send` buffer is 1D, because we have a halo width of 1. But for
            // CGVectors we have `nCells = CGdegree` which effectively means our effective halo
            // width is `CGdegree`. This means the send buffer will be of dimensions:
            // buffer_len x nCells
            //
            // buffer_map maps the raw data pointer to a 2D Eigen array with the correct dimensions.
            Eigen::Map<Eigen::ArrayXXd> buffer_map(send[comp].data(), buffer_len, nCells);
            for (auto edge : edges) {
                // populate edge data into send buffer
                sourceToSendBuffer(edge, buffer_map, source_map);
            }
        }
    }

    /**
     * @brief Calculate which edge the send position sendPos would fall under in the send buffer
     *
     * @param sendPos Position in the send buffer
     * @param fromRank Rank from which the data is being sent
     * @return Edge enum indicating which edge the position corresponds to
     */
    Edge edgeFromSendPos(int sendPos, int fromRank);

    /**
     * @brief Calculate recv buffer positions and offsets (count, disp and recvOffset) for halo
     * transfer
     *
     * @param fromRank Rank from which to receive data
     * @param count Number of elements to receive
     * @param disp Displacement in the send buffer
     * @param recvOffset Offset in the receive buffer
     * @param edge Edge for which to calculate positions
     * @param neighbourIndex Index of the current neighbor
     * @param cell Cell index for CGVector fields
     * @param isPeriodic Whether the boundary is periodic
     */
    void recvPositions(int& fromRank, size_t& count, size_t& disp, size_t& recvOffset, Edge edge,
        const size_t neighbourIndex, const size_t cell, bool isPeriodic);

    /**
     * @brief Calculate recv buffer positions and offsets (count, disp and recvOffset) for halo
     * transfer of corner cells
     *
     * @param fromRank Rank from which to receive data
     * @param count Number of elements to receive
     * @param disp Displacement in the send buffer
     * @param recvOffset Offset in the receive buffer
     * @param corner Corner for which to calculate positions
     * @param neighbourIndex Index of the current neighbor
     * @param cell Cell index for CGVector fields
     * @param isPeriodic Whether the boundary is periodic
     */
    void recvPositions(int& fromRank, size_t& count, size_t& disp, size_t& recvOffset,
        Corner corner, const size_t cell, bool isPeriodic);

    /*!
     * @brief Populate recv buffer with halo data from other ranks send buffers via the memory
     * window
     */
    void populateRecvBuffers();

    /**
     * @brief Populate target array with data from the receive buffer
     *
     * @param target Target array to populate with received halo data
     */
    template <typename T> void populateTarget(T& target)
    {
        size_t buffer_len = recvBufferSize;
        if (isCG) {
            buffer_len = recvBufferSize / CGdegree;
        }
        for (size_t comp = 0; comp < m_numComps; ++comp) {
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<Eigen::Dynamic>> target_map(
                target.col(comp).data(), m_Nx, m_Ny, Eigen::InnerStride<>(m_numComps));
            Eigen::Map<Eigen::ArrayXXd> buffer_map(recv[comp].data(), buffer_len, nCells);
            for (auto edge : edges) {
                recvBufferToTarget(edge, buffer_map, target_map);
            }
            for (auto corner : corners) {
                recvBufferToTarget(corner, buffer_map, target_map);
            }
        }
    }

    /*!
     * @brief Transpose corners received from vertical edges
     *
     * If the sending edge was from a vertical side (e.g., LEFT or RIGHT) then it has been
     * transposed when flattened into the 1D send buffer. We need to account for that.
     */
    void transposeCorners();

    /*!
     * @brief Calculate send buffer positions for a given edge
     *
     * @param idx_a Starting index along given edge
     * @param idx_b Ending index along given edge
     * @param offset Offset position in the send buffer
     * @param edge Edge for which to calculate positions
     */
    void sendBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge);

    /**
     * @brief Populate send buffer with data from source array for a given edge
     *
     * Buffer maps are used to get both arrays in the same format.
     *
     * @param edge Edge for which to populate the send buffer
     * @param buffer_map maps send arrays to Eigen::Arrays
     * @param source_map maps source to similar Eigen::Array as buffer_map
     */
    template <typename S, typename T>
    void sourceToSendBuffer(const Edge edge, Eigen::Map<T>& buffer_map,
        Eigen::Map<S, 0, Eigen::InnerStride<Eigen::Dynamic>>& source_map)
    {
        int idx_a, idx_b, offset;
        sendBufferPositions(idx_a, idx_b, offset, edge);
        if (isVertical(edge)) {
            buffer_map.transpose()(Eigen::all, Eigen::seq(offset, offset + m_innerNy - 1))
                = source_map(Eigen::seq(idx_a, idx_b), Eigen::seq(nCells, Eigen::last - nCells));
        } else {
            buffer_map(Eigen::seq(offset, offset + m_innerNx - 1), Eigen::all)
                = source_map(Eigen::seq(nCells, Eigen::last - nCells), Eigen::seq(idx_a, idx_b));
        }
    }

    /*!
     * @brief Calculate receive buffer positions for a given edge
     *
     * @param idx_a Starting index along given edge
     * @param idx_b Ending index along given edge
     * @param offset Offset position in the receive buffer
     * @param edge Edge for which to calculate positions
     */
    void recvBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge);

    /**
     * @brief Copy data from receive buffer to target array for a given edge
     *
     * Buffer maps are used to get both arrays in the same format.
     *
     * @param edge Edge for which to copy data from receive buffer to target
     * @param buffer_map Mapped receive buffer containing the data
     * @param target_map Mapped target array to store the received data
     */
    template <typename S, typename T>
    void recvBufferToTarget(const Edge edge, Eigen::Map<T>& buffer_map,
        Eigen::Map<S, 0, Eigen::InnerStride<Eigen::Dynamic>>& target_map)
    {
        int idx_a, idx_b, offset;
        recvBufferPositions(idx_a, idx_b, offset, edge);
        if (isVertical(edge)) {
            target_map(Eigen::seq(idx_a, idx_b), Eigen::seq(nCells, Eigen::last - nCells))
                = buffer_map.transpose()(Eigen::all, Eigen::seq(offset, offset + m_innerNy - 1));
        } else {
            target_map(Eigen::seq(nCells, Eigen::last - nCells), Eigen::seq(idx_a, idx_b))
                = buffer_map(Eigen::seq(offset, offset + m_innerNx - 1), Eigen::all);
        }
    }

    /**
     * @brief Copy data from receive buffer to target array for a given corner
     *
     * Buffer maps are used to get both arrays in the same format.
     *
     * @param corner Corner for which to copy data from receive buffer to target
     * @param buffer_map Mapped receive buffer containing the data
     * @param target_map Mapped target array to store the received data
     */
    template <typename S, typename T>
    void recvBufferToTarget(const Corner corner, Eigen::Map<T>& buffer_map,
        Eigen::Map<S, 0, Eigen::InnerStride<Eigen::Dynamic>>& target_map)
    {
        int perimeter = 2 * (m_innerNx + m_innerNy);
        int offset = perimeter + corner * nCells;

        switch (corner) {
        case Corner::TOP_LEFT:
            target_map(Eigen::seq(0, nCells - 1), Eigen::seq(Eigen::last - nCells + 1, Eigen::last))
                = buffer_map(Eigen::seq(offset, offset + nCells - 1), Eigen::all);
            break;
        case Corner::TOP_RIGHT:
            target_map(Eigen::seq(Eigen::last - nCells + 1, Eigen::last),
                Eigen::seq(Eigen::last - nCells + 1, Eigen::last))
                = buffer_map(Eigen::seq(offset, offset + nCells - 1), Eigen::all);
            break;
        case Corner::BOTTOM_RIGHT:
            target_map(Eigen::seq(Eigen::last - nCells + 1, Eigen::last), Eigen::seq(0, nCells - 1))
                = buffer_map(Eigen::seq(offset, offset + nCells - 1), Eigen::all);
            break;
        case Corner::BOTTOM_LEFT:
            target_map(Eigen::seq(0, nCells - 1), Eigen::seq(0, nCells - 1))
                = buffer_map(Eigen::seq(offset, offset + nCells - 1), Eigen::all);
            break;
        default:
            throw std::runtime_error("Unrecognised corner type");
        }
    }

public:
    /**
     * @brief Check if dimension is in the lateral direction i.e., X, Y, XVERTEX, or YVERTEX
     *
     * Determines whether a dimension represents a lateral (horizontal) direction,
     * which includes both cell centers and vertices in X and Y dimensions.
     *
     * @param dim The dimension to check
     * @return true if dimension is X, Y, XVERTEX, or YVERTEX, false otherwise
     */
    static bool isDimLateral(const ModelArray::Dimension dim)
    {
        return dim == ModelArray::Dimension::X || dim == ModelArray::Dimension::Y
            || dim == ModelArray::Dimension::XVERTEX || dim == ModelArray::Dimension::YVERTEX;
    }

    /*!
     * @brief Returns size of the inner flattened array
     */
    size_t getInnerSize() { return m_innerNx * m_innerNy; }

    /*!
     * @brief Get inner block from source array and copy into target
     *
     * @param source Source ModelArray-like object to extract inner block from
     * @param target Target ModelArray-like object to store the extracted inner block
     */
    template <typename S, typename T = S> void getInnerBlock(S& source, T& target)
    {
        ArraySlicer::Slice::VBounds sourceSlice, targetSlice;
        sourceSlice = { { haloWidth, -haloWidth }, { haloWidth, -haloWidth } };
        targetSlice = { {}, {} };
        ArraySlicer::SliceIter sourceIter(sourceSlice, { m_Nx, m_Ny });
        ArraySlicer::SliceIter targetIter(targetSlice, { m_innerNx, m_innerNy });

        ModelArraySlice::copySliceWithIters(source, sourceIter, target, targetIter);
    }

    /*!
     * @brief Set inner block of the target array from a given source
     *
     * @param source Source ModelArray-like object containing data to set in the inner block
     * @param target Target ModelArray-like object to have its inner block set
     */
    template <typename S, typename T = S> void setInnerBlock(S& source, T& target)
    {
        ArraySlicer::Slice::VBounds sourceSlice, targetSlice;
        sourceSlice = { {}, {} };
        targetSlice = { { haloWidth, -haloWidth }, { haloWidth, -haloWidth } };
        ArraySlicer::SliceIter sourceIter(sourceSlice, { m_innerNx, m_innerNy });
        ArraySlicer::SliceIter targetIter(targetSlice, { m_Nx, m_Ny });

        target.setZero(); // everything outside the inner block should be initialized to zero.

        ModelArraySlice::copySliceWithIters(source, sourceIter, target, targetIter);
    }

    /**
     * @brief Exchanges halo regions between neighboring MPI ranks
     *
     * Performs a complete halo exchange operation by:
     * 1. Populating send buffers with local halo data
     * 2. Receiving halo data from neighboring ranks into receive buffers
     * 3. Transposing corner data for CGVector fields (if needed)
     * 4. Populating the target array with the received halo data
     *
     * @param target The target array to update with exchanged halo data
     */
    template <typename T> void exchangeHalos(T& target)
    {
        populateSendBuffers(target);
        populateRecvBuffers();
        if (isCG) {
            transposeCorners();
        }
        populateTarget(target);
    }
};
} // end of nextsim namespace

#endif /* HALO_HPP */
