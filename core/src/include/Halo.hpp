/*!
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef HALO_HPP
#define HALO_HPP

#include <cstddef>
#include <iostream>
#include <numeric>
#include <vector>

#include "Slice.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelMPI.hpp"
#include "include/ModelMetadata.hpp"
#include "mpi.h"

#ifndef HALOWIDTH
#define HALOWIDTH 1
#endif

#ifndef DGCOMP
#define DGCOMP 3
#endif

#ifndef DGSTRESSCOMP
#define DGSTRESSCOMP 8
#endif

#ifndef CGDEGREE
#define CGDEGREE 2
#endif

namespace Nextsim {

/*!
 * @brief A class to facilitate halo exchange between MPI ranks
 *
 * @details
 */
class Halo {
public:
    /*!
     * @brief Constructs a halo object
     */
    Halo(size_t numComps, bool isVertex = false)
        : m_numComps(numComps)
        , isVertex(isVertex)
    {
        setSpatialDims();
        intializeHaloMetadata();
    }

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

        // inner dimension of domain excluding the halo cells
        m_Nx = m_innerNx + 2 * haloWidth;
        m_Ny = m_innerNy + 2 * haloWidth;
    }

    void intializeHaloMetadata()
    {
        // number of halo cells (should be general for any halo width)
        m_numHaloCells = 2 * haloWidth * (m_innerNx + m_innerNy + 2 * haloWidth);

        // need send / recv buffers for each componenet (e.g., each DGCOMP)
        send.resize(m_numComps);
        recv.resize(m_numComps);
        for (size_t i = 0; i < m_numComps; i++) {
            // allocate size and initialize to zero
            send[i].resize(m_numHaloCells, 0.0);
            recv[i].resize(m_numHaloCells, 0.0);
        }

        // order is Bottom, Right, Top, Left
        m_edgeLengths = { m_innerNx, m_innerNy, m_innerNx, m_innerNy };

        m_outerSlices = {
            { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 0 } }) },
            { Edge::RIGHT, VBounds({ { -1 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -1 } }) },
            { Edge::LEFT, VBounds({ { 0 }, { 1, m_innerNy + haloWidth } }) },
        };
        if (isVertex) {
            m_innerSlices = {
                { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 2 } }) },
                { Edge::RIGHT, VBounds({ { -3 }, { 1, m_innerNy + haloWidth } }) },
                { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -3 } }) },
                { Edge::LEFT, VBounds({ { 2 }, { 1, m_innerNy + haloWidth } }) },
            };
        } else {
            m_innerSlices = {
                { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 1 } }) },
                { Edge::RIGHT, VBounds({ { -2 }, { 1, m_innerNy + haloWidth } }) },
                { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -2 } }) },
                { Edge::LEFT, VBounds({ { 1 }, { 1, m_innerNy + haloWidth } }) },
            };
        }
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

    std::array<size_t, Edge::N_EDGE> m_edgeLengths; // array containing length of each edge
    std::array<Edge, Edge::N_EDGE> edges = ModelMetadata::edges; // array of edge enums
    std::map<Edge, Slice> m_outerSlices;
    std::map<Edge, Slice> m_innerSlices;
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

    bool isVertex; // some ModelArrays can be of type VERTEX

    std::vector<std::vector<double>>
        send; // buffer to store halo region that will be read by other ranks
    std::vector<std::vector<double>>
        recv; // buffer to store halo region which is read from other ranks
    ModelArray::DataType tempBuffer;

    /*!
     * @brief Open memory window to exchange send buffer between MPI ranks.
     *
     * @ details this is not intended to be used manually. It should only be called as part of the
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
     * @ details this is not intended to be used manually. It should only be called as part of the
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
        tempBuffer.resize(m_numHaloCells, m_numComps);
        tempBuffer = 0.;

        for (auto edge : edges) {
            size_t beg = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
            size_t num = m_edgeLengths.at(edge);

            Slice sourceSlice = m_innerSlices.at(edge);
            SliceIter sourceIter(sourceSlice, { m_Nx, m_Ny });

            if (isVertical(edge)) {
                Slice tempBufferSlice = { { {}, { beg, beg + num } } };
                // spatial dims are flattened into 1-D.
                SliceIter tempBufferIter(tempBufferSlice, { 1, m_numHaloCells });
                ModelArraySlice::copySliceWithIters(source, sourceIter, tempBuffer, tempBufferIter);
            } else {
                Slice tempBufferSlice = { { { beg, beg + num }, {} } };
                // spatial dims are flattened into 1-D.
                SliceIter tempBufferIter(tempBufferSlice, { m_numHaloCells, 1 });
                ModelArraySlice::copySliceWithIters(source, sourceIter, tempBuffer, tempBufferIter);
            }
        }
        // we need to copy into std vector send buffer because MPI doesn't work with Eigen
        // Arrays
        for (size_t i = 0; i < m_numComps; i++) {
            typedef Eigen::Map<ModelArray::DataType> MapType;
            MapType map(send[i].data(), m_numHaloCells, 1);
            // map is connected with the send buffer so the following line sets the data in send
            map = tempBuffer.col(i);
        }
    }

    /*!
     * @brief Populate recv buffer with halo data from other ranks send buffers via the memory
     * window
     */
    void populateRecvBuffers()
    {

        // do halo exchange for each component
        for (size_t comp = 0; comp < m_numComps; comp++) {
            // open memory window to send buffer on other ranks
            openMemoryWindow(comp);
            auto& metadata = ModelMetadata::getInstance();

            // get non-periodic neighbours and populate recv buffer (if the exist)
            for (auto edge : edges) {
                auto numNeighbours = metadata.neighbourRanks[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        int fromRank = metadata.neighbourRanks[edge][i];
                        size_t count = metadata.neighbourExtents[edge][i];
                        size_t disp = metadata.neighbourHaloSend[edge][i];
                        size_t recvOffset = metadata.neighbourHaloRecv[edge][i];
                        if (isVertex) {
                            vertexAdjustedPositions(
                                count = count, disp = disp, recvOffset = recvOffset, edge = edge);
                        }
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }
            }

            // get periodic neighbours and populate recv buffer (if they exist)
            for (auto edge : edges) {
                auto numNeighbours = metadata.neighbourRanksPeriodic[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        int fromRank = metadata.neighbourRanksPeriodic[edge][i];
                        size_t count = metadata.neighbourExtentsPeriodic[edge][i];
                        size_t disp = metadata.neighbourHaloSendPeriodic[edge][i];
                        size_t recvOffset = metadata.neighbourHaloRecvPeriodic[edge][i];
                        if (isVertex) {
                            vertexAdjustedPositions(
                                count = count, disp = disp, recvOffset = recvOffset, edge = edge);
                        }
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }
            }

            // close memory window (essentially make sure all communications are done before
            // moving on)
            closeMemoryWindow();
        }

        // copy from the recv std vector into an eigen array
        for (size_t i = 0; i < m_numComps; i++) {
            typedef Eigen::Map<ModelArray::DataType> MapType;
            MapType map(recv[i].data(), m_numHaloCells, 1);
            // map is connected with the recv buffer so the following line copies the data into
            // tempBuffer
            tempBuffer.col(i) = map;
        }
    }

    void vertexAdjustedPositions(size_t& count, size_t& disp, size_t& recvOffset, Edge& edge)
    {
        count = count + haloWidth;
        disp = disp + oppositeEdge.at(edge) * haloWidth;
        recvOffset = recvOffset + edge * haloWidth;
    }

public:
    /*!
     * @brief Returns size of the inner flattened array
     */
    size_t getInnerSize() { return m_innerNx * m_innerNy; }

    /*!
     * @brief Get inner block of ModelArray from tempData
     *
     * @params ma ModelArray which we intend to update across MPI ranks
     */
    template <typename S, typename T = S> void getInnerBlock(S& source, T& target)
    {
        ArraySlicer::Slice::VBounds sourceSlice, targetSlice;
        sourceSlice = { { HALOWIDTH, -HALOWIDTH }, { HALOWIDTH, -HALOWIDTH } };
        targetSlice = { {}, {} };
        ArraySlicer::SliceIter sourceIter(sourceSlice, { m_Nx, m_Ny });
        ArraySlicer::SliceIter targetIter(targetSlice, { m_innerNx, m_innerNy });

        ModelArraySlice::copySliceWithIters(source, sourceIter, target, targetIter);
    }

    template <typename S, typename T = S> void setInnerBlock(S& source, T& target)
    {
        ArraySlicer::Slice::VBounds sourceSlice, targetSlice;
        sourceSlice = { {}, {} };
        targetSlice = { { HALOWIDTH, -HALOWIDTH }, { HALOWIDTH, -HALOWIDTH } };
        ArraySlicer::SliceIter sourceIter(sourceSlice, { m_innerNx, m_innerNy });
        ArraySlicer::SliceIter targetIter(targetSlice, { m_Nx, m_Ny });

        // target = 0.; // everything outside the inner block should be initialized to zero.

        ModelArraySlice::copySliceWithIters(source, sourceIter, target, targetIter);
    }

    /*!
     * @brief Update a ModelArray with data from
     *
     * @params dgvec DGVector which we intend to update across MPI ranks based on halo cells
     * note that the start index is offset by 1 and the loop limit is size() - 2 because
     * the edge of each domain is 2 less than the length of the expanded halo region
     * (see diagram below - the empty cells are skipped by going from i+1 to size()-2)
     * ┌─┬─┬─┬─┐
     * │ │x│x│ │
     * ├─┼─┼─┼─┤
     * │x│o│o│x│
     * ├─┼─┼─┼─┤
     * │x│o│o│x│       o = original data
     * ├─┼─┼─┼─┤       x = mpi halo data (from recv)
     * │ │x│x│ │ (empty) = unused data in DGVector
     * └─┴─┴─┴─┘
     */
    template <typename T> void exchangeHalos(T& target)
    {
        populateSendBuffers(target);
        populateRecvBuffers();

        for (auto edge : edges) {
            size_t beg = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
            size_t num = m_edgeLengths.at(edge);

            Slice targetSlice = m_outerSlices.at(edge);
            SliceIter targetIter(targetSlice, { m_Nx, m_Ny });

            if (isVertical(edge)) {
                Slice tempBufferSlice = { { {}, { beg, beg + num } } };
                // spatial dims are flattened into 1-D.
                SliceIter tempBufferIter(tempBufferSlice, { 1, m_numHaloCells });
                ModelArraySlice::copySliceWithIters(tempBuffer, tempBufferIter, target, targetIter);
            } else {
                Slice tempBufferSlice = { { { beg, beg + num }, {} } };
                // spatial dims are flattened into 1-D.
                SliceIter tempBufferIter(tempBufferSlice, { m_numHaloCells, 1 });
                ModelArraySlice::copySliceWithIters(tempBuffer, tempBufferIter, target, targetIter);
            }
        }
    }
};
} // end of nextsim namespace

#endif /* HALO_HPP */
