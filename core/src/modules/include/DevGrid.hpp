/*!
 * @file DevGrid.hpp
 *
 * @date Dec 20, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#ifndef DEVGRID_HPP
#define DEVGRID_HPP

#include "include/IStructure.hpp"

#include "include/IDevGridIO.hpp"
#include "include/ModelState.hpp"

#include <map>

namespace Nextsim {

class DevGridIO;

//! A class to hold a grid of ElementData instances in a fixed sized square grid.
class DevGrid : public IStructure {
public:
    DevGrid()
        : pio(nullptr)
    {
    }

    //! Destructor. The lifetime of pio should be the lifetime of the instance.
    virtual ~DevGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    const static int nx;
    const static std::string structureName;

    // Read/write override functions
#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, const std::string& partitionFile) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#else
    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#endif

    void dumpModelState(
        const ModelState& state, const ModelMetadata& metadata, const std::string& filePath, bool isRestart = false) const override
    {
        if (pio)
            pio->dumpModelState(state, metadata, filePath, isRestart);
    }
    const std::string& structureType() const override { return structureName; };

    int nIceLayers() const override { return 1; };

    //! Sets the pointer to the class that will perform the IO. Should be an instance of DevGridIO
    void setIO(IDevGridIO* p) { pio = p; }

    const static std::string xDimName;
    const static std::string yDimName;
    const static std::string nIceLayersName;

private:
    IDevGridIO* pio;

    friend DevGridIO;
};

} /* namespace Nextsim */

#endif /* DEVGRID_HPP */
