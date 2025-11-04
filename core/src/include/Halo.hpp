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
#include <iostream>
#include <numeric>
#include <ostream>
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
     */
    template <int N> Halo(DGVectorHolder<N>& dgvh)
    {
        m_numComps = N;
        setSpatialDims();
        intializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from DGVector
     */
    template <int N> Halo(DGVector<N>& dgv)
    {
        m_numComps = N;
        setSpatialDims();
        intializeHaloMetadata();
    }

    /*!
     * @brief Constructs a halo object from CGVector
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
        send.resize(m_numComps);
        recv.resize(m_numComps);
        for (size_t i = 0; i < m_numComps; i++) {
            // allocate size and initialize to zero
            send[i].resize(m_numHaloCells, 0.0);
            recv[i].resize(m_numHaloCells, 0.0);
        }

        // order is Bottom, Right, Top, Left
        m_edgeLengths = { m_innerNx, m_innerNy, m_innerNx, m_innerNy };
    }

    static const size_t haloWidth = HALOWIDTH; // how many cells wide is the halo region

private:
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;
    using Edge = ModelMetadata::Edge;
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
    std::map<Edge, Edge> oppositeEdge = {
        { Edge::LEFT, Edge::RIGHT },
        { Edge::RIGHT, Edge::LEFT },
        { Edge::TOP, Edge::BOTTOM },
        { Edge::BOTTOM, Edge::TOP },
    }; // map to opposite edge

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

    std::vector<std::vector<double>>
        send; // buffer to store halo region that will be read by other ranks
    std::vector<std::vector<double>>
        recv; // buffer to store halo region which is read from other ranks
    ModelArray::DataType tempBuffer;

    /*!
     * @brief Open memory window to exchange send buffer between MPI ranks.
     *
     * @details this is not intended to be used manually. It should only be called as part of the
     * update method.
     */
    void openMemoryWindow(size_t idx)
    {
        // create a RMA memory window which all ranks will be able to access
        auto& modelMPI = ModelMPI::getInstance();
        MPI_Win_create(&send[idx][0], m_numHaloCells * sizeof(double), sizeof(double),
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
        size_t buffer_len = m_numHaloCells;
        if (isCG) {
            buffer_len = m_numHaloCells / CGdegree;
        }
        for (size_t comp = 0; comp < m_numComps; ++comp) {
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<Eigen::Dynamic>> source_map(
                source.col(comp).data(), m_Nx, m_Ny, Eigen::InnerStride<>(m_numComps));
            Eigen::Map<Eigen::ArrayXXd> buffer_map(send[comp].data(), buffer_len, nCells);
            for (auto edge : edges) {
                sourceToSendBuffer(edge, buffer_map, source_map);
            }
        }
    }

    /**
     * @brief Calculate recv buffer positions and offsets for halo transfer
     *
     * @param[in,out] count The count of elements to be communicated, increased by haloWidth
     * @param[in,out] disp The displacement value, adjusted based on opposite edge
     * @param[in,out] recvOffset The receive offset, adjusted based on edge
     * @param[in] edge The edge along which communication occurs
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
        if (isVertex) {
            count = count + 1;
            disp = disp + oppositeEdge.at(edge);
            recvOffset = recvOffset + edge;
        }
        if (isCG) {
            count = CGdegree * count + 1;
            disp = (disp > 0) ? CGdegree * disp + oppositeEdge.at(edge) : 0;
            recvOffset = (recvOffset > 0) ? CGdegree * recvOffset + edge : 0;

            // recvOffset is the offset in the recv buffer and this belongs to the current rank
            recvOffset = recvOffset + m_numHaloCells / nCells * cell;

            // disp is the offset in the "sending" buffer which belongs to rank "fromRank"
            // Therefore we need to compute how many halo cells that rank has to work out the offset
            // for each cell
            auto extentX = CGdegree * metadata.getRankExtentsX()[fromRank] + 1;
            auto extentY = CGdegree * metadata.getRankExtentsY()[fromRank] + 1;
            auto fromRank_numHaloCells
                = 2 * haloWidth * CGdegree * (extentX + extentY + 2 * haloWidth * CGdegree);
            disp = disp + fromRank_numHaloCells / nCells * cell;
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

    template <typename T> void populateTarget(T& target)
    {
        size_t buffer_len = m_numHaloCells;
        if (isCG) {
            buffer_len = m_numHaloCells / CGdegree;
        }
        for (size_t comp = 0; comp < m_numComps; ++comp) {
            Eigen::Map<Eigen::ArrayXXd, 0, Eigen::InnerStride<Eigen::Dynamic>> target_map(
                target.col(comp).data(), m_Nx, m_Ny, Eigen::InnerStride<>(m_numComps));
            Eigen::Map<Eigen::ArrayXXd> buffer_map(recv[comp].data(), buffer_len, nCells);
            for (auto edge : edges) {
                recvBufferToTarget(edge, buffer_map, target_map);
            }
        }
    }

    /*!
     * @brief Get inner block of ModelArray from tempData
     *
     * @params ma ModelArray which we intend to update across MPI ranks
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
        populateTarget(target);
    }
};
} // end of nextsim namespace

#endif /* HALO_HPP */
