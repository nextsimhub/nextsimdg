/*!
 * @file OASISCoupled.hpp
 *
 * @date Aug 15, 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 */

#ifndef OASISCOUPLED_HPP
#define OASISCOUPLED_HPP

#include "include/ModelMetadata.hpp"

namespace Nextsim {

class OASISCoupled {
public:
    ~OASISCoupled()
    {
        // TODO: Insert the oasis terminate call here
        /*
         * const std::string functionName = "~OASISCoupled";
         * const std::string message = "couldn't terminate OASIS";
         * const std::string file = "OASISCoupled.hpp";
         * if ( ! oasis_c_terminate() )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 26);
         */
    }

    virtual std::string getName() const { return "OASISCoupled"; }

    virtual void setMetadata(const ModelMetadata& metadata)
    {
        // TODO: Insert OASIS initialisation calls here
        // Set the communicators
        /* This is commented for now, as I don't have OASIS ready on my system
         * These three may need to be global to the class
         * int compID;
         * MPI_Comm localComm;
         * MPI_Comm coupledComm;
         *
         * const std::string compName = "nextsim";
         * const std::string functionName = getName()+"::setMetadata";
         * const std::string message = "couldn't initialise component";
         * const std::string file = "OASISCoupled.hpp";
         * const bool coupled = true;
         * if ( ! oasis_c_init_comp_with_comm(&compID, &compName.c_str(), coupled, metadata.mpiComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 47);
         *
         * if ( ! oasis_c_get_localcomm(&localComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 50);
         *
         * const int icpl = 1;
         * if ( ! oasis_c_create_couplcomm(icpl, &localComm, &coupledComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 54);
         */

        // Set the partitioning
        /* This is commented for now, as I don't have OASIS ready on my system
         *
         * From the manual: "vector of integers describing the local grid partition in the global
         * index space; has a different expression depending on the type of the partition; in
         * OASIS3-MCT, 5 types of partition are supported: Serial (no partition), Apple, Box,
         * Orange, and Points" - it looks like we should use "Box", so partInfo[0] = 2.
         *
         * metdatata contains: localCornerX, localCornerY, localExtentX, localExtentY, globalExtentX,
         * globalExtentY;
         *
         * int partitionID;
         * const int offset = metadata.localExtentX*metadata.localCornerY + metadata.localCornerX;
         * const std::vector<int> partInfo {2, offset, metadata.localExtentX, metadata.localExtentY, metadata.globalExtentX};
         *
         * igSize and name are optional and seem not relevant for our setup. Using -1 and "" in the
         * c-interface means they are not sent further when the Fortran code is called.
         * const int igSize = -1;
         * if ( ! oasis_c_def_partition(&partitionID, &paral.data(), paral.size(), igSize, &compName.c_str()) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 76);
         *
         * def_var and end_def calls are called by the child class
         */
    }
};

}

#endif // OASISCOUPLED_HPP
