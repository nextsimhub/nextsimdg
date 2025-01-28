/*!
 * @file Halo.hpp
 *
 * @date 29 Apr 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef HALO_HPP
#define HALO_HPP

#include <cstddef>
#include <iostream>
#include <memory>
#include <numeric>
#include <vector>

#include "Slice.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelMetadata.hpp"
#include "mpi.h"

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
    Halo(ModelMetadata& metadata, ModelArray& ma)
        : m_metadata(std::make_unique<ModelMetadata>(metadata))
        , m_ma(ma)
    {
        m_innerNx = ma.innerDimensions()[0];
        m_innerNy = ma.innerDimensions()[1];
        m_perimeterLength = 2 * m_innerNx + 2 * m_innerNy;
        m_numComps = m_ma.nComponents();
        send.resize(m_numComps);
        recv.resize(m_numComps);
        for (size_t i = 0; i < m_numComps; i++) {
            send[i].resize(m_perimeterLength, 0.0);
            recv[i].resize(m_perimeterLength, 0.0);
        }
        m_edgeLengths = { m_innerNx, m_innerNy, m_innerNx, m_innerNy }; // order is Bottom
        m_comm = metadata.mpiComm;

        m_outerSlices = {
            { Edge::LEFT, VBounds({ { 0 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::RIGHT, VBounds({ { -1 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -1 } }) },
            { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 0 } }) },
        };
        m_innerSlices = {
            { Edge::LEFT, VBounds({ { 1 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::RIGHT, VBounds({ { -2 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -2 } }) },
            { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 1 } }) },
        };
        m_innerSlicesVertexAdjusted = {
            { Edge::LEFT, VBounds({ { 2 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::RIGHT, VBounds({ { -3 }, { 1, m_innerNy + haloWidth } }) },
            { Edge::TOP, VBounds({ { 1, m_innerNx + haloWidth }, { -3 } }) },
            { Edge::BOTTOM, VBounds({ { 1, m_innerNx + haloWidth }, { 2 } }) },
        };
    }

    static const size_t haloWidth = 1; // how many cells wide is the halo region

private:
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;
    using Edge = ModelMetadata::Edge;
    using VBounds = ArraySlicer::Slice::VBounds;

    const typedef ArraySlicer::SliceIter::MultiDim MultiDim;
    const MultiDim m_haloDims;
    ModelArray& m_ma; // reference to modelarray
    size_t m_innerNx; // local extent in x-direction
    size_t m_innerNy; // local extent in y-direction
    size_t m_perimeterLength; // length of perimeter of domain
    size_t m_numComps; // number of DG components
    std::unique_ptr<ModelMetadata> m_metadata; // pointer to metadata
    std::array<size_t, Edge::N_EDGE> m_edgeLengths; // array containing length of each edge
    std::array<Edge, Edge::N_EDGE> edges = ModelMetadata::edges; // array of edge enums
    std::map<Edge, Slice> m_outerSlices;
    std::map<Edge, Slice> m_innerSlices;
    std::map<Edge, Slice> m_innerSlicesVertexAdjusted;
    MPI_Win m_win; // RMA memory window object (used for sharing send buffers between ranks)
    MPI_Comm m_comm; // RMA memory window object (used for sharing send buffers between ranks)

    std::map<Edge, Edge> oppositeEdge = {
        { Edge::LEFT, Edge::RIGHT },
        { Edge::RIGHT, Edge::LEFT },
        { Edge::TOP, Edge::BOTTOM },
        { Edge::BOTTOM, Edge::TOP },
    }; // map to opposite edge

    /*!
     * @brief Open memory window to exchange send buffer between MPI ranks.
     *
     * @ details this is not intended to be used manually. It should only be called as part of the
     * update method.
     */
    void openMemoryWindow(size_t idx)
    {
        // create a RMA memory window which all ranks will be able to access
        MPI_Win_create(&send[idx][0], m_perimeterLength * sizeof(double), sizeof(double),
            MPI_INFO_NULL, m_comm, &m_win);
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
    void populateSendBuffers()
    {
        tempBuffer.resize(m_perimeterLength, m_numComps);
        for (auto edge : edges) {
            size_t beg = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
            size_t num = m_edgeLengths.at(edge);
            if (m_ma.getType() == ModelArray::Type::VERTEX) {
                // because vertex points lie along the domain boundaries we need offset the
                // slices by and additional row/column
                // m_ma[m_innerSlicesVertexAdjusted.at(edge)].copyToBuffer(send[i], offset);
                tempBuffer(Eigen::seqN(beg, num), Eigen::all)
                    = static_cast<ModelArray::DataType>(m_ma[m_innerSlicesVertexAdjusted.at(edge)]);
            } else {
                // m_ma[m_innerSlices.at(edge)].copyToBuffer(send[i], offset);
                tempBuffer(Eigen::seqN(beg, num), Eigen::all)
                    = static_cast<ModelArray::DataType>(m_ma[m_innerSlices.at(edge)]);
            }
        }
        // we need to copy into std vector send buffer because MPI doesn't work with Eigen Arrays
        for (size_t i = 0; i < m_numComps; i++) {
            typedef Eigen::Map<ModelArray::DataType> MapType;
            MapType map(send[i].data(), m_perimeterLength, 1);
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

            // get non-periodic neighbours and populate recv buffer (if the exist)
            for (auto edge : edges) {
                auto numNeighbours = m_metadata->neighbourRanks[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        int fromRank = m_metadata->neighbourRanks[edge][i];
                        size_t count = m_metadata->neighbourExtents[edge][i];
                        size_t disp = m_metadata->neighbourHaloSend[edge][i];
                        size_t recvOffset = m_metadata->neighbourHaloRecv[edge][i];
                        if (m_ma.getType() == ModelArray::Type::VERTEX) {
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
                auto numNeighbours = m_metadata->neighbourRanksPeriodic[edge].size();
                if (numNeighbours) {
                    // get data for each neighbour that exists along a given edge
                    for (size_t i = 0; i < numNeighbours; ++i) {
                        int fromRank = m_metadata->neighbourRanksPeriodic[edge][i];
                        size_t count = m_metadata->neighbourExtentsPeriodic[edge][i];
                        size_t disp = m_metadata->neighbourHaloSendPeriodic[edge][i];
                        size_t recvOffset = m_metadata->neighbourHaloRecvPeriodic[edge][i];
                        if (m_ma.getType() == ModelArray::Type::VERTEX) {
                            vertexAdjustedPositions(
                                count = count, disp = disp, recvOffset = recvOffset, edge = edge);
                        }
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }
            }

            // close memory window (essentially make sure all communications are done before moving
            // on)
            closeMemoryWindow();
        }

        // copy from the recv std vector into an eigen array
        for (size_t i = 0; i < m_numComps; i++) {
            typedef Eigen::Map<ModelArray::DataType> MapType;
            MapType map(recv[i].data(), m_perimeterLength, 1);
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
    std::vector<std::vector<double>>
        send; // buffer to store halo region that will be read by other ranks
    std::vector<std::vector<double>>
        recv; // buffer to store halo region which is read from other ranks
    ModelArray::DataType tempBuffer;

    /*!
     * @brief Returns size of the inner flattened array
     */
    size_t getInnerSize() { return m_innerNx * m_innerNy; }

    /*!
     * @brief Populate inner block of ModelArray from tempData
     *
     * @params ma ModelArray which we intend to update across MPI ranks
     */
    void populateInnerBlock(ModelArray::DataType& tempData)
    {
        ArraySlicer::Slice::VBounds innerBlock, allBlock;
        if (m_ma.getType() == ModelArray::Type::Z) {
            innerBlock = { { 1, -1 }, { 1, -1 }, {} };
            allBlock = { {}, {}, {} };
        } else {
            innerBlock = { { 1, -1 }, { 1, -1 } };
            allBlock = { {}, {} };
        }
        ArraySlicer::SliceIter wholeBlock(allBlock, m_ma.innerDimensions());

        m_ma = 0.; // TODO -- should this be removed? It does mean that mask is zero by default

        // copy temp data to the central block of the main modelarray
        m_ma[innerBlock].copyFromSlicedBuffer(tempData, wholeBlock);
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
    void exchangeHalos()
    {
        populateSendBuffers();
        populateRecvBuffers();

        for (auto edge : edges) {
            size_t beg = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
            size_t num = m_edgeLengths.at(edge);
            // m_ma[m_outerSlices.at(edge)] = tempBuffer(Eigen::seqN(beg, num), Eigen::all);
            m_ma[m_outerSlices.at(edge)]
                = ModelArray::DataType(tempBuffer(Eigen::seqN(beg, num), Eigen::all));
        }
    }
};
} // end of nextsim namespace

#endif /* HALO_HPP */
