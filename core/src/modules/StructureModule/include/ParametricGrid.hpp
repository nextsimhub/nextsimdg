/*!
 * @file ParametricGrid.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 * @author Kacper Kornet <kk562@cam.ac.uk>
 */

#ifndef PARAMETRICGRID_HPP
#define PARAMETRICGRID_HPP

#include "include/IStructure.hpp"

#include "include/ModelState.hpp"

#include <set>

namespace Nextsim {

class ParaGridIO;

//! A class to hold the grid data for parameterised rectangular grids.
class ParametricGrid : public IStructure {
public:
    typedef ModelArray::MultiDim GridDimensions;
    ParametricGrid()
        : pio(nullptr)
    {
    }
    virtual ~ParametricGrid()
    {
        if (pio) {
            delete pio;
        }
    }

    // Read/write override functions
#ifdef USE_MPI
    ModelState getModelState(const std::string& filePath, const std::string& partitionFile,
        ModelMetadata& metadata) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#else
    ModelState getModelState(const std::string& filePath) override
    {
        return pio ? pio->getModelState(filePath) : ModelState();
    }
#endif

    void dumpModelState(const ModelState& state, const ModelMetadata& metadata,
        const std::string& filePath, bool isRestart = false) const override
    {
        if (pio) {
            if (isRestart) {
                pio->dumpModelState(state, metadata, filePath);
            } else {
                pio->writeDiagnosticTime(state, metadata, filePath);
            }
        }
    }
    const std::string& structureType() const override { return structureName; };

    int nIceLayers() const override
    {
        return ModelArray::definedDimensions.at(ModelArray::Dimension::Z).length;
    };

    class IParaGridIO {
    public:
        IParaGridIO(ParametricGrid& grid)
            : grid(grid)
        {
        }
        virtual ~IParaGridIO() = default;

        virtual ModelState getModelState(const std::string& filePath) = 0;
        virtual void dumpModelState(
            const ModelState& state, const ModelMetadata& metadata, const std::string& filePath)
            = 0;
        virtual ModelState readForcingTime(const std::set<std::string>& forcings,
            const TimePoint& time, const std::string& filePath)
            = 0;
        virtual void writeDiagnosticTime(
            const ModelState& state, const ModelMetadata& meta, const std::string& filePath)
            = 0;

    protected:
        IParaGridIO() = delete;

    private:
        ParametricGrid& grid;
    };

    /*!
     * Sets the IParaGridIO instance to use to drive the input and output.
     * @param p Pointer to the instance to use.
     */
    void setIO(IParaGridIO* p) { pio = p; }

    const static std::string structureName;

private:
    IParaGridIO* pio;

    friend ParaGridIO;
};

} /* namespace Nextsim */

#endif /* PARAMETRICGRID_HPP */
