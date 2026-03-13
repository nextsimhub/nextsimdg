/*!
 * @file StructureFactory.hpp
 *
 * @date Mar 11, 2026
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef MOCK_INCLUDES_STRUCTUREFACTORY_HPP
#define MOCK_INCLUDES_STRUCTUREFACTORY_HPP

#include "include/ModelState.hpp"
#include "include/StructureFactory.hpp"
class StructureFactory {
public:
    static Nextsim::ModelState& getState()
    {
        static Nextsim::ModelState storedState;
        return storedState;
    }
    static void fileFromState(const Nextsim::ModelState& state, const std::string& filePath, bool isRestart)
    {
        getState() = state;
    }
};

#endif /* MOCK_INCLUDES_STRUCTUREFACTORY_HPP */
