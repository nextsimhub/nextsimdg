/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Halo exchange class
 * @details
 *
 * Halo exchange class
 *
 * All functionality for halo exchange between MPI ranks is contained in this class.
 *
 * Halo supports the main data structures of NextSim e.g., ModelArray, DGVector and CGVector.
 *
 * The halos are exchange via one-sided MPI communication using RMA.
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
        intializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from DGVector
     * @param dgvh DGVectorHolder object to create halo from
     */
    template <int N> Halo(DGVectorHolder<N>& dgvh)
    {
        m_numComps = N;
        setSpatialDims();
        intializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from DGVector
     * @param dgv DGVector object to create halo from
     */
    template <int N> Halo(DGVector<N>& dgv)
    {
        m_numComps = N;
        setSpatialDims();
        intializeHaloMetadata();
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
        intializeHaloMetadata();
    }

    /**
     * @brief Sets the spatial dimensions of the domain
     *
     * This function sets the inner dimensions of the computational domain based on metadata,
     * adjusts dimensions for vertex fields if necessary, and calculates total dimensions
     * including halo cells.
     */
    void setSpatialDims()
    {
        auto& metadata = ModelMetadata::getInstance();

        // spatial dimension of domain
        m_innerNx = metadata.getLocalExtentX();
        m_innerNy = metadata.getLocalExtentY();

        // extend dimensions by 1 for vertex fields
        if (isVertex) {
            m_innerNx += 1;
            m_innerNy += 1;
        }
        // extend dimensions for CGVectors
        if (isCG) {
            m_innerNy = m_innerNy * CGdegree + 1;
            m_innerNx = m_innerNx * CGdegree + 1;
        }

        // inner dimension of domain excluding the halo cells
        m_Nx = m_innerNx + 2 * haloWidth;
        m_Ny = m_innerNy + 2 * haloWidth;

        // each additional halo cell add CGdegree more points
        if (isCG) {
            m_Nx = m_innerNx + 2 * haloWidth * CGdegree;
            m_Ny = m_innerNy + 2 * haloWidth * CGdegree;
        }
    }

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
    void intializeHaloMetadata()
    {
        // number of halo cells (should be general for any halo width)
        if (not isCG) {
            m_numHaloCells = 2 * haloWidth * (m_innerNx + m_innerNy + 2 * haloWidth);
        } else {
            m_numHaloCells
                = 2 * haloWidth * CGdegree * (m_innerNx + m_innerNy + 2 * haloWidth * CGdegree);
        }

        // need send / recv buffers for each component (e.g., each DGCOMP)
        sendBufferSize = m_numHaloCells - 4 * nCells * nCells;
        recvBufferSize = m_numHaloCells;
        send.resize(m_numComps);
        recv.resize(m_numComps);
        for (size_t i = 0; i < m_numComps; i++) {
            // allocate size and initialize to zero
            send[i].resize(sendBufferSize, 0.0);
            recv[i].resize(recvBufferSize, 0.0);
        }

        // order is Bottom, Right, Top, Left
        m_edgeLengths = { m_innerNx, m_innerNy, m_innerNx, m_innerNy };
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
    void openMemoryWindow(size_t idx)
    {
        // create a RMA memory window which all ranks will be able to access
        auto& modelMPI = ModelMPI::getInstance();
        MPI_Win_create(&send[idx][0], sendBufferSize * sizeof(double), sizeof(double),
            MPI_INFO_NULL, modelMPI.getComm(), &m_win);
        // remove fence and check that no proceding RMA calls have been made
        MPI_Win_fence(MPI_MODE_NOPRECEDE, m_win);
    }

    /*!
     * @brief Initialise memory window to exchange send buffer between MPI ranks.
     *
     * @details this is not intended to be used manually. It should only be called as part of the
     * update method.
     */
    void closeMemoryWindow()
    {
        // enable fence i.e., disable future RMA calls until we re-open memory window
        MPI_Win_fence(MPI_MODE_NOSUCCEED, m_win);
        // free window object
        MPI_Win_free(&m_win);
    }

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
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<Eigen::Dynamic>> source_map(
                source.col(comp).data(), m_Nx, m_Ny, Eigen::InnerStride<>(m_numComps));
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
     */
    Edge edgeFromSendPos(int sendPos, int fromRank)
    {
        auto& metadata = ModelMetadata::getInstance();

        // extents of sending domain
        auto extentX = metadata.getRankExtentsX()[fromRank];
        auto extentY = metadata.getRankExtentsY()[fromRank];

        if (sendPos - (2 * extentX + extentY) >= 0) {
            return Edge::LEFT;
        } else if (sendPos - (extentX + extentY) >= 0) {
            return Edge::TOP;
        } else if (sendPos - extentX >= 0) {
            return Edge::RIGHT;
        } else {
            return Edge::BOTTOM;
        }
    }

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
        const size_t neighbourIndex, const size_t cell, bool isPeriodic)
    {
        auto& metadata = ModelMetadata::getInstance();
        if (isPeriodic) {
            fromRank = metadata.neighbourRanksPeriodic[edge][neighbourIndex];
            count = metadata.neighbourExtentsPeriodic[edge][neighbourIndex];
            disp = metadata.neighbourHaloSendPeriodic[edge][neighbourIndex];
            recvOffset = metadata.neighbourHaloRecvPeriodic[edge][neighbourIndex];
        } else {
            fromRank = metadata.neighbourRanks[edge][neighbourIndex];
            count = metadata.neighbourExtents[edge][neighbourIndex];
            disp = metadata.neighbourHaloSend[edge][neighbourIndex];
            recvOffset = metadata.neighbourHaloRecv[edge][neighbourIndex];
        }
        auto sendEdge = edgeFromSendPos(disp, fromRank);
        if (isVertex) {
            count = count + 1;
            disp = disp + sendEdge;
            recvOffset = recvOffset + edge;
        }
        if (isCG) {
            count = CGdegree * count + 1;
            disp = (disp > 0) ? CGdegree * disp + sendEdge : 0;
            recvOffset = (recvOffset > 0) ? CGdegree * recvOffset + edge : 0;

            // recvOffset is the offset in the recv buffer and this belongs to the current rank
            recvOffset = recvOffset + recvBufferSize / nCells * cell;

            // disp is the offset in the "sending" buffer which belongs to rank "fromRank"
            // Therefore we need to compute how many halo cells that rank has to work out the offset
            // for each cell
            auto extentX = CGdegree * metadata.getRankExtentsX()[fromRank] + 1;
            auto extentY = CGdegree * metadata.getRankExtentsY()[fromRank] + 1;
            // TODO this is too big. Should only be 2 x CGdegree * haloWidth (extentX + extentY)
            auto fromRankSendBufferSize = 2 * haloWidth * CGdegree * (extentX + extentY);
            disp = disp + fromRankSendBufferSize / nCells * cell;
        }
    }

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
        Corner corner, const size_t cell, bool isPeriodic)
    {
        count = 1; // we only have maximum of 1 corner neighbour for each corner
        auto& metadata = ModelMetadata::getInstance();
        if (isPeriodic) {
            fromRank = metadata.cornerRanksPeriodic[corner][0];
            disp = metadata.cornerHaloSendPeriodic[corner][0];
            recvOffset = metadata.cornerHaloRecvPeriodic[corner][0];
        } else {
            fromRank = metadata.cornerRanks[corner][0];
            disp = metadata.cornerHaloSend[corner][0];
            recvOffset = metadata.cornerHaloRecv[corner][0];
        }
        auto sendEdge = edgeFromSendPos(disp, fromRank);
        if (isVertex) {
            count = 1;
            disp = disp + sendEdge;
            recvOffset = recvOffset + 4;
            // this is messy
            // Essentially account for the fact that the vertex field is split differently to the
            // face centered fields
            if ((sendEdge == Edge::TOP or sendEdge == Edge::BOTTOM)
                and (corner == Corner::TOP_RIGHT or corner == Corner::BOTTOM_RIGHT)) {
                disp = disp + 1;
            }
            if ((sendEdge == Edge::LEFT or sendEdge == Edge::RIGHT)
                and (corner == Corner::TOP_RIGHT or corner == Corner::TOP_LEFT)) {
                disp = disp + 1;
            }
        }
        if (isCG) {
            count = CGdegree * count;
            disp = (disp > 0) ? CGdegree * disp + sendEdge : 0;
            recvOffset = CGdegree * recvOffset + 4;

            // recvOffset is the offset in the recv buffer and this belongs to the current rank
            recvOffset = recvOffset + recvBufferSize / nCells * cell;

            // disp is the offset in the "sending" buffer which belongs to rank "fromRank"
            // Therefore we need to compute how many halo cells that rank has to work out the offset
            // for each cell
            auto extentX = CGdegree * metadata.getRankExtentsX()[fromRank] + 1;
            auto extentY = CGdegree * metadata.getRankExtentsY()[fromRank] + 1;
            auto fromRankSendBufferSize = 2 * haloWidth * CGdegree * (extentX + extentY);
            disp = disp + fromRankSendBufferSize / nCells * cell;
            if ((sendEdge == Edge::TOP or sendEdge == Edge::BOTTOM)
                and (corner == Corner::TOP_RIGHT or corner == Corner::BOTTOM_RIGHT)) {
                disp = disp + 1;
            }
            if ((sendEdge == Edge::LEFT or sendEdge == Edge::RIGHT)
                and (corner == Corner::TOP_RIGHT or corner == Corner::TOP_LEFT)) {
                disp = disp + 1;
            }
        }
    }

    /*!
     * @brief Populate recv buffer with halo data from other ranks send buffers via the memory
     * window
     */
    void populateRecvBuffers()
    {

        int fromRank;
        size_t count, disp, recvOffset;
        // do halo exchange for each component
        for (size_t comp = 0; comp < m_numComps; comp++) {
            // open memory window to send buffer on other ranks
            openMemoryWindow(comp);
            auto& metadata = ModelMetadata::getInstance();
            // get non-periodic neighbours and populate recv buffer (if the exist)
            for (auto edge : edges) {

                // get neighbours (if they exist)
                auto numNeighbours = metadata.neighbourRanks[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        for (size_t cell = 0; cell < nCells; ++cell) {
                            recvPositions(fromRank, count, disp, recvOffset, edge, i, cell, false);
                            MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp,
                                count, MPI_DOUBLE, m_win);
                        }
                    }
                }

                // get periodic neighbours (if they exist)
                numNeighbours = metadata.neighbourRanksPeriodic[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        for (size_t cell = 0; cell < nCells; ++cell) {
                            recvPositions(fromRank, count, disp, recvOffset, edge, i, cell, true);
                            MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp,
                                count, MPI_DOUBLE, m_win);
                        }
                    }
                }
            }

            // get non-periodic corner neighbours and populate recv buffer (if the exist)
            for (auto corner : corners) {

                // get neighbours (if they exist)
                auto hasCorner = metadata.cornerRanks[corner].size();
                // hasCorner will either be 0 or 1
                if (hasCorner) {
                    for (size_t cell = 0; cell < nCells; ++cell) {
                        recvPositions(fromRank, count, disp, recvOffset, corner, cell, false);
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }

                // get periodic neighbours (if they exist)
                hasCorner = metadata.cornerRanksPeriodic[corner].size();
                if (hasCorner) {
                    for (size_t cell = 0; cell < nCells; ++cell) {
                        recvPositions(fromRank, count, disp, recvOffset, corner, cell, true);
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }
            }

            // close memory window (essentially make sure all communications are done before
            // moving on)
            closeMemoryWindow();
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
    static const bool isDimLateral(const ModelArray::Dimension dim)
    {
        return dim == ModelArray::Dimension::X || dim == ModelArray::Dimension::Y
            || dim == ModelArray::Dimension::XVERTEX || dim == ModelArray::Dimension::YVERTEX;
    }

    /*!
     * @brief Returns size of the inner flattened array
     */
    size_t getInnerSize() { return m_innerNx * m_innerNy; }

    /*!
     * @brief Calculate send buffer positions for a given edge
     *
     * @param idx_a Starting index along given edge
     * @param idx_b Ending index along given edge
     * @param offset Offset position in the send buffer
     * @param edge Edge for which to calculate positions
     */
    void sendBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge)
    {
        int vertexOffset = 0;
        if (isVertex || isCG) {
            vertexOffset = 1;
        }
        offset = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
        switch (edge) {
        case Edge::LEFT:
            idx_a = nCells + vertexOffset;
            idx_b = 2 * nCells - 1 + vertexOffset;
            break;
        case Edge::RIGHT:
            idx_a = m_Nx - 2 * nCells - vertexOffset;
            idx_b = m_Nx - nCells - 1 - vertexOffset;
            break;
        case Edge::BOTTOM:
            idx_a = nCells + vertexOffset;
            idx_b = 2 * nCells - 1 + vertexOffset;
            break;
        case Edge::TOP:
            idx_a = m_Ny - 2 * nCells - vertexOffset;
            idx_b = m_Ny - nCells - 1 - vertexOffset;
            break;
        default:
            throw std::runtime_error("Unrecognised edge type");
        }
    }

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
    void recvBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge)
    {
        offset = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
        switch (edge) {
        case Edge::LEFT:
            idx_a = 0;
            idx_b = nCells - 1;
            break;
        case Edge::RIGHT:
            idx_a = m_Nx - nCells;
            idx_b = m_Nx - 1;
            break;
        case Edge::BOTTOM:
            idx_a = 0;
            idx_b = nCells - 1;
            break;
        case Edge::TOP:
            idx_a = m_Ny - nCells;
            idx_b = m_Ny - 1;
            break;
        default:
            throw std::runtime_error("Unrecognised edge type");
        }
    }

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

    /*!
     * @brief Transpose corners received from vertical edges
     *
     * If the sending edge was from a vertical side (e.g., LEFT or RIGHT) then it has been
     * transposed when flattened into the 1D send buffer. We need to account for that.
     */
    void transposeCorners()
    {
        auto& metadata = ModelMetadata::getInstance();
        int fromRank;
        size_t disp;
        for (auto corner : corners) {
            // non-periodic corners
            bool hasCorner = false;
            bool isPeriodic = false;

            if (metadata.cornerRanks[corner].size() > 0) {
                // non-periodic case
                hasCorner = true;
                isPeriodic = false;
            } else if (metadata.cornerRanksPeriodic[corner].size() > 0) {
                // periodic case
                hasCorner = true;
                isPeriodic = true;
            }

            if (hasCorner) {
                if (isPeriodic) {
                    fromRank = metadata.cornerRanksPeriodic[corner][0];
                    disp = metadata.cornerHaloSendPeriodic[corner][0];
                } else {
                    fromRank = metadata.cornerRanks[corner][0];
                    disp = metadata.cornerHaloSend[corner][0];
                }

                auto sendEdge = edgeFromSendPos(disp, fromRank);

                if (isVertical(sendEdge)) {
                    auto buffer_len = recvBufferSize / CGdegree;
                    for (size_t comp = 0; comp < m_numComps; ++comp) {
                        Eigen::Map<Eigen::ArrayXXd> rmap(recv[comp].data(), buffer_len, nCells);
                        auto offset = buffer_len - (4 - corner) * nCells;
                        rmap.block(offset, 0, nCells, nCells).transposeInPlace();
                    }
                }
            }
        }
    }

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
