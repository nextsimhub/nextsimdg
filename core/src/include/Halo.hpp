/*!
 * @file Halo.hpp
 *
 * @date 28 Jan 2025
 * @author Tom Meltzer <tdm39@cam.ac.uk>
 */

#ifndef HALO_HPP
#define HALO_HPP

#include <memory>
#include <numeric>
#include <vector>

#include "Slice.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelArraySlice.hpp"
#include "include/ModelMetadata.hpp"
#include "include/dgVector.hpp"
#include "include/indexer.hpp"
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
    Halo(size_t localExtentX, size_t localExtentY, ModelMetadata& metadata)
        : m_localExtentX(localExtentX)
        , m_localExtentY(localExtentY)
        , m_metadata(std::make_unique<ModelMetadata>(metadata))
        , m_comm(metadata.mpiComm)
        , m_haloDims({ localExtentX + 2 * m_halo_width, localExtentY + 2 * m_halo_width })
    {
        m_perimeterLength = 2 * localExtentX + 2 * localExtentY;
        send.resize(m_perimeterLength, 0.0);
        recv.resize(m_perimeterLength, 0.0);
        m_edgeLengths
            = { localExtentX, localExtentY, localExtentX, localExtentY }; // order is Bottom
        // Right Top Left
    }

private:
    using Slice = ArraySlicer::Slice;
    using SliceIter = ArraySlicer::SliceIter;
    using Edge = ModelMetadata::Edge;
    using VBounds = ArraySlicer::Slice::VBounds;

    const size_t m_halo_width = 1; // how many cells wide is the halo region
    const typedef ArraySlicer::SliceIter::MultiDim MultiDim;
    const MultiDim m_haloDims;
    size_t m_localExtentX; // local extent in x-direction
    size_t m_localExtentY; // local extent in y-direction
    size_t m_perimeterLength; // length of perimeter of domain
    std::unique_ptr<ModelMetadata> m_metadata; // pointer to metadata
    std::array<size_t, Edge::N_EDGE> m_edgeLengths; // array containing length of each edge
    std::array<Edge, Edge::N_EDGE> edges = ModelMetadata::edges; // array of edge enums
    std::map<Edge, Slice> m_slices = {
        { Edge::LEFT, VBounds({ { 0 }, {} }) },
        { Edge::RIGHT, VBounds({ { -1 }, {} }) },
        { Edge::TOP, VBounds({ {}, { -1 } }) },
        { Edge::BOTTOM, VBounds({ {}, { 0 } }) },
    };

    MPI_Win m_win; // RMA memory window object (used for sharing send buffers between ranks)
    MPI_Comm m_comm; // RMA memory window object (used for sharing send buffers between ranks)

    /*!
     * @brief Open memory window to exchange send buffer between MPI ranks.
     *
     * @ details this is not intended to be used manually. It should only be called as part of the
     * update method.
     */
    void openMemoryWindow()
    {
        // create a RMA memory window which all ranks will be able to access
        MPI_Win_create(&send[0], m_perimeterLength * sizeof(double), sizeof(double), MPI_INFO_NULL,
            m_comm, &m_win);
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

public:
    std::vector<double> send; // buffer to store halo region that will be read by other ranks
    std::vector<double> recv; // buffer to store halo region which is read from other ranks

    /*!
     * @brief Populate send buffer with halo data of the specified ModelArray
     *
     * @params ma ModelArray which we intend to update across MPI ranks
     */
    void populateSendBuffer(ModelArray& ma)
    {
        for (auto edge : ModelMetadata::edges) {
            size_t offset = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);
            ma[m_slices.at(edge)].copyToBuffer(send, offset);
        }
    }

    /*!
     * @brief Populate recv buffer with halo data from other ranks send buffers via the memory
     * window
     */
    void populateRecvBuffer()
    {

        // open memory window to send buffer on other ranks
        openMemoryWindow();

        // get non-periodic neighbours and populate recv buffer (if the exist)
        for (auto edge : ModelMetadata::edges) {
            auto numNeighbours = m_metadata->neighbourRanks[edge].size();
            if (numNeighbours) {
                // get data for each neighbour that exists along a given edge
                for (size_t i = 0; i < numNeighbours; ++i) {
                    int fromRank = m_metadata->neighbourRanks[edge][i];
                    size_t count = m_metadata->neighbourExtents[edge][i];
                    size_t disp = m_metadata->neighbourHaloSend[edge][i];
                    size_t recvOffset = m_metadata->neighbourHaloRecv[edge][i];
                    MPI_Get(&recv[recvOffset], count, MPI_DOUBLE, fromRank, disp, count, MPI_DOUBLE,
                        m_win);
                }
            }
        }

        // get periodic neighbours and populate recv buffer (if they exist)
        for (auto edge : ModelMetadata::edges) {
            auto numNeighbours = m_metadata->neighbourRanksPeriodic[edge].size();
            if (numNeighbours) {
                // get data for each neighbour that exists along a given edge
                for (size_t i = 0; i < numNeighbours; ++i) {
                    int fromRank = m_metadata->neighbourRanksPeriodic[edge][i];
                    size_t count = m_metadata->neighbourExtentsPeriodic[edge][i];
                    size_t disp = m_metadata->neighbourHaloSendPeriodic[edge][i];
                    size_t recvOffset = m_metadata->neighbourHaloRecvPeriodic[edge][i];
                    MPI_Get(&recv[recvOffset], count, MPI_DOUBLE, fromRank, disp, count, MPI_DOUBLE,
                        m_win);
                }
            }
        }

        // close memory window (essentially make sure all communications are done before moving on)
        closeMemoryWindow();
    }

    /*!
     * @brief Update a DGVector with data from the recv buffer
     *
     * @params dgvec DGVector which we intend to update across MPI ranks based on halo cells
     */
    void updateDGVec(DGVector<DGCOMP>& dgvec)
    {
        for (auto edge : edges) {

            SliceIter sIter = SliceIter(m_slices.at(edge), m_haloDims);
            std::vector<size_t> edgeIndices;

            // populate edgeIndices with the indices along a given edge of the domain
            while (!sIter.isEnd()) {
                const size_t start = sIter.index();
                const size_t step = sIter.step(0);
                const size_t n = sIter.nElements(0);
                for (int i = 0; i < n; ++i) {
                    auto idx = start + i * step;
                    edgeIndices.push_back(idx);
                }
                sIter.incrementDim(1);
            }

            // calculate offset index for the recv buffer based on current edge
            const size_t offset
                = std::accumulate(m_edgeLengths.begin(), m_edgeLengths.begin() + edge, 0);

            // copy the halo region from recv buffer into the DGVector
            for (size_t i = 0; i < edgeIndices.size() - 2; ++i) {
                // note that the start index is offset by 1 and the loop limit is size() - 2 because
                // the edge of each domain is 2 less than the length of the expanded halo region
                // (see diagram below - the empty cells are skipped by going from i+1 to size()-2)
                // ┌─┬─┬─┬─┐
                // │ │x│x│ │
                // ├─┼─┼─┼─┤
                // │x│o│o│x│
                // ├─┼─┼─┼─┤
                // │x│o│o│x│       o = original data
                // ├─┼─┼─┼─┤       x = mpi halo data (from recv)
                // │ │x│x│ │ (empty) = unused data in DGVector
                // └─┴─┴─┴─┘
                dgvec(edgeIndices[i + 1], 0) = recv[offset + i];
            }
        }

        recv.clear();
        send.clear();
    }
};
} // end of nextsim namespace

#endif /* HALO_HPP */
