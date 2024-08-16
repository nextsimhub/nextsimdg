/*!
 * @file OASISCoupled.hpp
 *
 * @date Aug 15, 2024
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
    int compID;

#ifdef USE_OASIS
    ~OASISCoupled()
    {
        const std::string functionName = __func__;
        const std::string message = "Couldn't terminate OASIS";
        const std::string file = __FILE__;
        if (!oasis_c_terminate())
            oasis_c_abort(compID, functionName.c_str(), message.c_str(), file.c_str(), 28);
    }
#else
    ~OASISCoupled() { }
#endif

    virtual std::string getName() const { return "OASISCoupled"; }

#ifdef USE_OASIS
    virtual void setMetadata(const ModelMetadata& metadata)
    {
        // TODO: Insert OASIS initialisation calls here
        // Set the communicators
        // These two may need to be global to the class
        MPI_Comm localComm;
        MPI_Comm coupledComm;

        const std::string compName = "nextsim";
        const std::string functionName = __func__;
        const std::string message = "Couldn't initialise component";
        const std::string file = __FILE__;
        const bool coupled = true;
        if (!oasis_c_init_comp_with_comm(&compID, compName.c_str(), coupled, metadata.mpiComm))
            oasis_c_abort(compID, functionName.c_str(), message.c_str(), file.c_str(), __LINE__);

        if (!oasis_c_get_localcomm(&localComm))
            oasis_c_abort(compID, functionName.c_str(), message.c_str(), file.c_str(), __LINE__);

        const int icpl = 1;
        if (!oasis_c_create_couplcomm(icpl, localComm, &coupledComm))
            oasis_c_abort(compID, functionName.c_str(), message.c_str(), file.c_str(), __LINE__);

        // Set the partitioning
        /* From the manual: "vector of integers describing the local grid partition in the global
         * index space; has a different expression depending on the type of the partition; in
         * OASIS3-MCT, 5 types of partition are supported: Serial (no partition), Apple, Box,
         * Orange, and Points" - it looks like we should use "Box", so partInfo[0] = 2.
         *
         * metdatata contains: localCornerX, localCornerY, localExtentX, localExtentY,
         * globalExtentX, globalExtentY;
         */
        int partitionID;
        const int offset = metadata.localExtentX * metadata.localCornerY + metadata.localCornerX;
        const std::vector<int> partInfo
            = { 2, offset, metadata.localExtentX, metadata.localExtentY, metadata.globalExtentX };

        /* igSize and name are optional and seem not relevant for our setup. Using -1 and "" in the
         * c-interface means they are not sent further when the Fortran code is called. */
        const int igSize = -1;
        const std::string name = "";
        if (!oasis_c_def_partition(
                &partitionID, partInfo.size(), &partInfo[0], igSize, name.c_str()))
            oasis_c_abort(compID, functionName.c_str(), message.c_str(), file.c_str(), __LINE__);

        // def_var and end_def calls are called by the child class
    }
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
