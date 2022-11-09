/*!
 * @file ParaGridIO.hpp
 *
 * @date Oct 24, 2022
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef PARAGRIDIO_HPP
#define PARAGRIDIO_HPP

#include "include/ParametricGrid.hpp"

#include <map>
#include <ncFile.h>
#include <string>

namespace Nextsim {

class ParaGridIO : public ParametricGrid::IParaGridIO {
public:
    ParaGridIO(ParametricGrid& grid)
        : IParaGridIO(grid)
    {
    }
    virtual ~ParaGridIO();

    ModelState getModelState(const std::string& filePath) override;

    void dumpModelState(const ModelState& state, const ModelMetadata& meta,
        const std::string& filePath, bool isRestart) override;

    void close(const std::string& filePath);

private:
    ParaGridIO() = delete;
    ParaGridIO(const ParaGridIO& other) = delete;
    ParaGridIO& operator=(const ParaGridIO& other) = delete;

    static const std::map<std::string, ModelArray::Type> dimensionKeys;

    static const std::map<ModelArray::Dimension, bool> isDG;
    static const std::map<ModelArray::Dimension, ModelArray::Type> dimCompMap;

    std::map<std::string, netCDF::NcFile> openFiles;
    std::map<std::string, size_t> timeIndexByFile;
};

} /* namespace Nextsim */

#endif /* PARAGRIDIO_HPP */
