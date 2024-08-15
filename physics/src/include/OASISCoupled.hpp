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
         * const std::string file = "OASISCoupledOcean.cpp";
         * const bool coupled = true;
         * if ( ! oasis_c_init_comp_with_comm(&compID, &compName.c_str(), coupled, metadata.mpiComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 33);
         *
         * if ( ! oasis_c_get_localcomm(&localComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 36);
         *
         * if ( ! oasis_c_create_couplcomm(metadata.mpiMyRank, &localComm, &coupledComm) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 39);
         */

        // Set the partitioning
        /* This is commented for now, as I don't have OASIS ready on my system
         * I don't know which values paral and igSize should take(!)
         * const std::vector<double> paral;
         * const int igSize = -1;
         * if ( ! oasis_c_def_partition(&partitionID, &paral.data(), paral.size(), igSize, &compName.c_str()) )
         *     oasis_c_abort(compID, &functionName.c_str(), &message.c_str(), &file.c_str(), 48);
         *
         * def_var and end_def calls are called by the child class
         */
    }
};

}

#endif // OASISCOUPLED_HPP
