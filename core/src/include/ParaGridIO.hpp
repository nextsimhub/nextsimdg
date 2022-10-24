/*!
 * @file ParaGridIO.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PARAGRIDIO_HPP
#define PARAGRIDIO_HPP

#include "include/ParametricGrid.hpp"

namespace Nextsim {

class ParaGridIO : public ParametricGrid::IParaGridIO {
public:
    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
    {
    }
    virtual ~ParaGridIO() = default;

    ModelState getModelState(const std::string& filePath) override;

    void dumpModelState(const ModelState& state, const ModelMetadata& meta,
        const std::string& filePath, bool isRestart) const override;

private:
    ParaGridIO() = default;
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
