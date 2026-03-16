/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @brief Halo exchange class implementation
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

#ifdef USE_MPI
#include "include/Halo.hpp"

namespace Nextsim {

using Edge = ModelMetadata::Edge;
using Corner = ModelMetadata::Corner;

void Halo::setSpatialDims()
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

void Halo::intializeHaloMetadata()
{
    // number of halo cells (should be general for any halo width)
    if (not isCG) {
        m_numHaloCells = 2 * haloWidth * (m_innerNx + m_innerNy + 2 * haloWidth);
    } else {
        m_numHaloCells
            = 2 * haloWidth * CGdegree * (m_innerNx + m_innerNy + 2 * haloWidth * CGdegree);
    }

    // need send / recv buffers for each component (e.g., each DGCOMP)
    sendBufferSize = m_numHaloCells - Corner::N_CORNER * nCells * nCells;
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

void Halo::openMemoryWindow(size_t idx)
{
    // create a RMA memory window which all ranks will be able to access
    auto& modelMPI = ModelMPI::getInstance();
    MPI_Win_create(&send[idx][0], sendBufferSize * sizeof(double), sizeof(double), MPI_INFO_NULL,
        modelMPI.getComm(), &m_win);
    // remove fence and check that no proceding RMA calls have been made
    MPI_Win_fence(MPI_MODE_NOPRECEDE, m_win);
}

void Halo::closeMemoryWindow()
{
    // enable fence i.e., disable future RMA calls until we re-open memory window
    MPI_Win_fence(MPI_MODE_NOSUCCEED, m_win);
    // free window object
    MPI_Win_free(&m_win);
}

Edge Halo::edgeFromSendPos(int sendPos, int fromRank)
{
    auto& metadata = ModelMetadata::getInstance();

    // extents of sending domain
    int extentX = metadata.getRankExtentsX()[fromRank];
    int extentY = metadata.getRankExtentsY()[fromRank];

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

void Halo::recvPositions(int& fromRank, size_t& count, size_t& disp, size_t& recvOffset, Edge edge,
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
        auto fromRankSendBufferSize = 2 * haloWidth * CGdegree * (extentX + extentY);
        disp = disp + fromRankSendBufferSize / nCells * cell;
    }
}

void Halo::recvPositions(int& fromRank, size_t& count, size_t& disp, size_t& recvOffset,
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
        recvOffset = recvOffset + Edge::N_EDGE;
        // Account for the fact that the vertex field is split differently to the face centered
        // fields. We dont take the data directly adjacent to the halo, but the one after that.
        // e.g., if you have two adjacent domains, the vertex on the far right of the left-hand
        // domain is the same as the vertex on the far left of the right-hand domain. We need
        // the vertex which is the next one along for halo exchange.
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
        recvOffset = CGdegree * recvOffset + Edge::N_EDGE;

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

void Halo::populateRecvBuffers()
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
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
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
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
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

void Halo::sendBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge)
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

void Halo::recvBufferPositions(int& idx_a, int& idx_b, int& offset, const Edge edge)
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

void Halo::transposeCorners()
{
    auto& metadata = ModelMetadata::getInstance();
    for (auto corner : corners) {
        // non-periodic corners
        bool hasCorner = false;
        bool isPeriodic = false;

        // Check if both periodic and non-periodic corner neighbours exist
        if (metadata.cornerRanksPeriodic[corner].size() > 0
            && metadata.cornerRanks[corner].size() > 0) {
            throw std::runtime_error(
                "It is not possible to have a non-periodic corner neighbour that is also "
                "periodic. Please check your partition metadata file.");
        }

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
            int fromRank;
            size_t disp;
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
                    auto offset = buffer_len - (Corner::N_CORNER - corner) * nCells;
                    rmap.block(offset, 0, nCells, nCells).transposeInPlace();
                }
            }
        }
    }
}

} // end of nextsim namespace
#endif // USE_MPI
