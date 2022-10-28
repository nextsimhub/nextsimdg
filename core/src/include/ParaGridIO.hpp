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
    ParaGridIO() = delete;

    static const std::map<std::string, ModelArray::Type> dimensionKeys;

    static const std::map<ModelArray::Dimension, bool> isDG;
    static const std::map<ModelArray::Dimension, ModelArray::Type> dimCompMap;
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
