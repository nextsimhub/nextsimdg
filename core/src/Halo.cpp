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

void Halo::setTripolarFlags()
{
    auto& metadata = ModelMetadata::getInstance();
    m_tripolarFold = metadata.needsTripolarFold();
}

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

void Halo::initializeHaloMetadata()
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

int Halo::recvPosFromEdge(Edge edge) const
{
    // extents of the local domain (buffer-map row unit, matching recvPositions)
    int extentX = m_innerNx;
    int extentY = m_innerNy;

    switch (edge) {
    case Edge::BOTTOM:
        return 0;
    case Edge::RIGHT:
        return extentX;
    case Edge::TOP:
        return extentX + extentY;
    case Edge::LEFT:
        return 2 * extentX + extentY;
    default:
        throw std::runtime_error("Halo :: Invalid edge enum");
    }
}

void Halo::recvPositions(int& fromRank, size_t& count, size_t& disp, size_t& recvOffset, Edge edge,
    const size_t neighbourIndex, const size_t cell)
{
    auto& metadata = ModelMetadata::getInstance();
    fromRank = metadata.neighbourRanks[edge][neighbourIndex];
    count = metadata.neighbourExtents[edge][neighbourIndex];
    disp = metadata.neighbourHaloSend[edge][neighbourIndex];
    recvOffset = metadata.neighbourHaloRecv[edge][neighbourIndex];
    auto sendEdge = edgeFromSendPos(disp, fromRank);
    if (isVertex) {
        recvOffset = recvOffset + edge;
        const bool isFirstTransaction = (recvOffset == recvPosFromEdge(edge));

        count = count + 1;
        disp = disp + sendEdge;

        if (!isFirstTransaction) {
            count = count - 1;
            recvOffset = recvOffset + 1;
            disp = disp + 1;

            // Explain the logic behind this.
            // I managed to convince myself there is one...
            // It is that we need to cut data from a different side
            // FIXME
            if (m_tripolarFold && edge == Edge::TOP) {
                disp = disp - 1;
            }
        }
    }
    if (isCG) {
        // Note that the CG memory transactions are overlapping
        // We need to make them non-overlapping to support th tripolar grid

        // We need to detect the first transaction along the edge (identified by the recv offset
        // matching the start of the edge). For the following transactions we need to shift
        // displacement by +1 and count by -1
        recvOffset = (recvOffset > 0) ? CGdegree * recvOffset + edge : 0;

        const bool isFirstTransaction = (recvOffset == recvPosFromEdge(edge));

        count = CGdegree * count + 1;
        disp = (disp > 0) ? CGdegree * disp + sendEdge : 0;

        if (!isFirstTransaction) {
            count = count - 1;
            recvOffset = recvOffset + 1;
            disp = disp + 1;

            // Explain the logic behind this.
            // I managed to convince myself there is one...
            // It is that we need to cut data from a different side
            // FIXME
            if (m_tripolarFold && edge == Edge::TOP) {
                disp = disp - 1;
            }
        }

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
    Corner corner, const size_t cell)
{
    count = 1; // we only have maximum of 1 corner neighbour for each corner
    auto& metadata = ModelMetadata::getInstance();
    fromRank = metadata.cornerRanks[corner][0];
    disp = metadata.cornerHaloSend[corner][0];
    recvOffset = metadata.cornerHaloRecv[corner][0];
    auto sendEdge = edgeFromSendPos(disp, fromRank);
    if (isVertex) {
        count = 1;
        disp = disp + sendEdge;
        recvOffset = recvOffset + Edge::N_EDGE;

        // FIXE: Dirty patch fix
        // This For the purpose of the patch below a Tripolar fold top corners
        // are flipped (since correct we communicate with the flipped image)
        // This monster below flips:
        //  TOP_LEFT -> BOTTOM_RIGHT
        //  TOP_RIGHT -> BOTTOM_LEFT
        // To be refactored. It makes eyes bleed and ears ringing as it is now...
        if (m_tripolarFold && (corner == Corner::TOP_LEFT || corner == Corner::TOP_RIGHT)) {
            corner = corner == Corner::TOP_LEFT ? Corner::BOTTOM_RIGHT : Corner::BOTTOM_LEFT;
        }

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
        // get neighbours and populate recv buffer (if the exist)
        for (auto edge : edges) {

            // get neighbours (if they exist)
            auto numNeighbours = metadata.neighbourRanks[edge].size();
            if (numNeighbours) {
                // get data for each neighbour that exists along a given edge
                for (size_t i = 0; i < numNeighbours; ++i) {
                    for (size_t cell = 0; cell < nCells; ++cell) {
                        recvPositions(fromRank, count, disp, recvOffset, edge, i, cell);
                        MPI_Get(&recv[comp][recvOffset], count, MPI_DOUBLE, fromRank, disp, count,
                            MPI_DOUBLE, m_win);
                    }
                }
            }
        }

        // get corner neighbours and populate recv buffer (if the exist)
        for (auto corner : corners) {

            // get neighbours (if they exist)
            auto hasCorner = metadata.cornerRanks[corner].size();
            // hasCorner will either be 0 or 1
            if (hasCorner) {
                for (size_t cell = 0; cell < nCells; ++cell) {
                    recvPositions(fromRank, count, disp, recvOffset, corner, cell);
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

void Halo::rotateTopEdgeInBuffer()
{
    // Protect against control flow errors.
    if (!m_tripolarFold) {
        throw std::runtime_error(
            "rotateTopEdgeInBuffer called on rank that does not need tripolar fold");
    }
    const auto& metadata = ModelMetadata::getInstance();
    const Edge edge = Edge::TOP;

    // We fixup the order of data component by component
    for (size_t comp = 0; comp < m_numComps; ++comp) {
        // We loop over all edge-based memory transactions on the top edge
        auto numNeighbours = metadata.neighbourRanks[edge].size();
        for (std::size_t i = 0; i < numNeighbours; ++i) {
            // Left -> Right flip
            for (size_t cell = 0; cell < nCells; ++cell) {
                int fromRank;
                std::size_t count, disp_ignore, recvOffset;
                recvPositions(fromRank, count, disp_ignore, recvOffset, edge, i, cell);
                m_tripolarFoldOp.flipCommTransaction(&recv[comp][recvOffset], count);
            }

            if (nCells > 1) {
                // In the buffer we have `nCells` number of rows that we need to swap
                // Hence we iterate over each pair, for odd `nCells`, the middle row is left
                // untouched.
                std::size_t midPoint = nCells / 2;
                for (std::size_t cell = 0; cell < midPoint; cell++) {
                    std::size_t otherCell = nCells - 1 - cell;

                    // Calculate the range in the buffer for each row
                    int fromRank;
                    std::size_t disp_ignore;
                    std::size_t count_1, recvOffset_1;
                    recvPositions(fromRank, count_1, disp_ignore, recvOffset_1, edge, i, cell);
                    std::size_t count_2, recvOffset_2;
                    recvPositions(fromRank, count_2, disp_ignore, recvOffset_2, edge, i, otherCell);

                    // Swap the data
                    FloatType* start_1 = &recv[comp][recvOffset_1];
                    FloatType* end_1 = start_1 + count_1;
                    FloatType* start_2 = &recv[comp][recvOffset_2];
                    std::swap_ranges(start_1, end_1, start_2);
                }
            }
        }
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
        bool hasCorner = metadata.cornerRanks[corner].size() > 0;

        if (hasCorner) {
            int fromRank = metadata.cornerRanks[corner][0];
            size_t disp = metadata.cornerHaloSend[corner][0];

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

namespace HaloExchange {

    void TripolarFold::flipCommTransaction(FloatType* data, std::size_t len) const
    {
        FloatType* start = data;
        FloatType* end = data + len;

        switch (m_dataType) {
        case DataType::VERTEX:
        case DataType::DG:
        case DataType::CG:
            std::reverse(start, end);
            break;
        default:
            // TODO: Use nextsim proper error handling conventions
            throw std::runtime_error("Unrecognised data type for tripolar fold");
        }
    }
}

} // end of nextsim namespace
#endif // USE_MPI
