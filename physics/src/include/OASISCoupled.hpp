/*!
 * @file OASISCoupled.hpp
 *
 * @date 29 Aug 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLED_HPP
#define OASISCOUPLED_HPP

#include "include/ModelMetadata.hpp"
#ifdef USE_OASIS
#include <oasis_c.h>
#endif

namespace Nextsim {

class OASISCoupled {
public:
#ifdef USE_OASIS
    ~OASISCoupled() { OASIS_CHECK_ERR(!oasis_c_terminate()); }
#else
    ~OASISCoupled() { }
#endif

    virtual std::string getName() const { return "OASISCoupled"; }

#ifdef USE_OASIS
    int partitionID, OASISTime;
    void setMetadata(const ModelMetadata& metadata)
    {
        // Set the "OASIS time" (seconds since start) to zero
        // TODO: We need to figure out if this is ok to do when restarting
        OASISTime = 0;

        // Set the communicators
        // Some of these three may need to be global to the class
        int compID;
        MPI_Comm localComm;
        MPI_Comm coupledComm;

        // compName could be configurable, but I don't see the need
        const std::string compName = "nextsim";
        OASIS_CHECK_ERR(oasis_c_init_comp_with_comm(
            &compID, compName.c_str(), OASIS_COUPLED, metadata.mpiComm));

        // Eric didn't include those, but it seems they should be necessary
        OASIS_CHECK_ERR(oasis_c_get_localcomm(&localComm));
        OASIS_CHECK_ERR(oasis_c_create_couplcomm(OASIS_COUPLED, localComm, &coupledComm));

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
        const int offset = metadata.globalExtentX * metadata.localCornerY + metadata.localCornerX;
        const std::vector<int> partInfo = { OASIS_Box, offset, metadata.localExtentX,
            metadata.localExtentY, metadata.globalExtentX };

        const int globalSize = metadata.globalExtentX * metadata.globalExtentY;
        OASIS_CHECK_ERR(oasis_c_def_partition(
            &partitionID, OASIS_Box_Params, &partInfo[0], globalSize, compName.c_str()));

        // TODO: Writing out grid information should be possible, but optional

        // def_var and end_def calls are called by the child class
    }

    // Increment the "OASIS" time by the number of seconds in the time step
    // Must be called at the end of the child class' update or updateAfter call.
    void updateOASISTime(const TimestepTime& tst) { OASISTime += tst.step.seconds(); }
#else
    virtual void setMetadata(const ModelMetadata& metadata)
    {
        std::string message = __func__ + std::string(": OASIS support not compiled in.\n");
        throw std::runtime_error(message);
    }
#endif
};

}

#endif // OASISCOUPLED_HPP
