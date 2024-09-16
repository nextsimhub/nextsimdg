/*!
 * @file ModelMetadata.cpp
 *
 * @date 10 Sep 2024
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/ModelMetadata.hpp"

#include "include/StructureModule.hpp"
#include "include/gridNames.hpp"

#ifdef USE_MPI
#include <ncDim.h>
#include <ncFile.h>
#include <ncGroup.h>
#include <ncVar.h>
#endif
#ifdef USE_OASIS
#include <oasis_c.h>
#endif

namespace Nextsim {

const std::string& ModelMetadata::structureName() const
{
    return Module::getImplementation<IStructure>().structureType();
}

#ifdef USE_MPI
ModelMetadata::ModelMetadata(std::string partitionFile, MPI_Comm comm)
{
    setMpiMetadata(comm);
    getPartitionMetadata(partitionFile);
}

void ModelMetadata::setMpiMetadata(MPI_Comm comm)
{
    mpiComm = comm;
    MPI_Comm_size(mpiComm, &mpiSize);
    MPI_Comm_rank(mpiComm, &mpiMyRank);
}

void ModelMetadata::getPartitionMetadata(std::string partitionFile)
{
    // TODO: Move the reading of the partition file to its own class
    netCDF::NcFile ncFile(partitionFile, netCDF::NcFile::read);
    int sizes = ncFile.getDim("L").getSize();
    int nBoxes = ncFile.getDim("P").getSize();
    if (nBoxes != mpiSize) {
        std::string errorMsg = "Number of MPI ranks " + std::to_string(mpiSize) + " <> "
            + std::to_string(nBoxes) + "\n";
        throw std::runtime_error(errorMsg);
    }
    globalExtentX = ncFile.getDim("NX").getSize();
    globalExtentY = ncFile.getDim("NY").getSize();
    netCDF::NcGroup bboxGroup(ncFile.getGroup(bboxName));
    std::vector<size_t> index(1, mpiMyRank);
    bboxGroup.getVar("domain_x").getVar(index, &localCornerX);
    bboxGroup.getVar("domain_y").getVar(index, &localCornerY);
    bboxGroup.getVar("domain_extent_x").getVar(index, &localExtentX);
    bboxGroup.getVar("domain_extent_y").getVar(index, &localExtentY);
    ncFile.close();
}

#endif

const ModelState& ModelMetadata::extractCoordinates(const ModelState& state)
{
    // More sophisticated grids include both vertex coordinates and grid azimuth values.
    if (state.data.count(coordsName) > 0) {
        m_vertexCoords = state.data.at(coordsName);
        m_gridAzimuth = state.data.at(gridAzimuthName);
        hasParameters = true;
    } else {
        // else don't resize the arrays
        hasParameters = false;
    }

    isCartesian = state.data.count(xName);
    if (isCartesian) {
        m_coord1 = state.data.at(xName);
        m_coord2 = state.data.at(yName);
    } else {
        m_coord1 = state.data.at(longitudeName);
        m_coord2 = state.data.at(latitudeName);
    }

    return state;
}

ModelState& ModelMetadata::affixCoordinates(ModelState& state) const
{
    if (hasParameters) {
        state.data[coordsName] = m_vertexCoords;
        state.data[gridAzimuthName] = m_gridAzimuth;
    }

    if (isCartesian) {
        state.data[xName] = m_coord1;
        state.data[yName] = m_coord2;
    } else {
        state.data[longitudeName] = m_coord1;
        state.data[latitudeName] = m_coord2;
    }
    return state;
}

#ifdef USE_OASIS
void ModelMetadata::initOasis(const bool writeOasisGrid)
{
    // Set the partitioning
    /* From the manual: "[ig_paral is a] vector of integers describing the local grid partition
     * in the global index space; has a different expression depending on the type of the
     * partition; in OASIS3-MCT, 5 types of partition are supported: Serial (no partition),
     * Apple, Box, Orange, and Points" - it looks like we should use "Box", so partInfo[0] = 2
     * (aka. ig_paral).
     *
     * #Box partition#
     * Each partition is a rectangular region of the global domain, described by the global
     * offset of its upper left corner, and its local extents in the X and Y dimensions. The
     * global extent in the X dimension must also be given. In this case, we have ig_paral(1:5):
     *  - ig_paral(1) = 2 (indicates a Box partition)
     *  - ig_paral(2) = the upper left corner global offset
     *  - ig paral(3) = the local extent in x
     *  - ig_paral(4) = the local extent in y
     *  - ig_paral(5) = the global extent in x.
     *
     * metdatata contains: localCornerX, localCornerY, localExtentX, localExtentY,
     * globalExtentX, globalExtentY;
     */
    // TODO: The contents of metadata is not certain!
    const int offset = globalExtentX * localCornerY + localCornerX;
    const std::vector<int> partInfo
        = { OASIS_Box, offset, localExtentX, localExtentY, globalExtentX };

    const int globalSize = globalExtentX * globalExtentY;
    const std::string compName = "nextsim"; // Not useful for any setups we have in mind
    OASIS_CHECK_ERR(oasis_c_def_partition(
        &OASISPartitionId, OASIS_Box_Params, &partInfo[0], globalSize, compName.c_str()));

    // TODO: Writing out grid information should be possible, but optional
    if (writeOasisGrid) {
        /* This needs to be figured out, but it's not a priority. Grid writing is
         * not necessary for the type of coupling we'll start with.

        const std::string gridName = "nxts";

        int flag = 1;
        OASIS_CHECK_ERR(oasis_c_start_grids_writing(&flag));

        OASIS_CHECK_ERR(oasis_c_write_grid(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, lon, lat, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_corner(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, clo, cla, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_area(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, area, OASISPartitionId));
        OASIS_CHECK_ERR(oasis_c_write_mask(
            gridName.c_str(), nx, ny, nx_loc, ny_loc, angle, OASISPartitionId));

        std::string companion = "land area fraction";
        OASIS_CHECK_ERR(oasis_c_write_frac(
                            gridName.c_str(), nx, ny, nx_loc, ny_loc, mask, OASISPartitionId),
            companion.c_str());
        companion = "land sea mask";
        OASIS_CHECK_ERR(oasis_c_write_mask(
                            gridName.c_str(), nx, ny, nx_loc, ny_loc, mask, OASISPartitionId),
            companion.c_str());
        */
    }
}
#endif

} /* namespace Nextsim */
