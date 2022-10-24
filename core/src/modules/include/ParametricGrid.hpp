/*!
 * @file ParametricGrid.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PARAMETRICGRID_HPP
#define PARAMETRICGRID_HPP

#include "IStructure.hpp"

#include "include/ModelState.hpp"

namespace Nextsim {

class ParaGridIO;

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

    class IParaGridIO {
    public:
        IParaGridIO(ParametricGrid& grid)
            : grid(grid)
        {
        }
        virtual ~IParaGridIO() = default;

        virtual ModelState getModelState(const std::string& filePath) = 0;
        virtual const ModelState& dumpModelState(const ModelState& state,
            const ModelMetadata& metadata, const std::string& filePath, bool isRestart) const = 0;

    protected:
        IParaGridIO() = default;

    private:
        ParametricGrid& grid;
    };

    void setIO(IParaGridIO* p) { pio = p; }

    const static std::string structureName;

private:
    IParaGridIO* pio;

    friend ParaGridIO;
};

} /* namespace Nextsim */

#endif /* PARAMETRICGRID_HPP */
