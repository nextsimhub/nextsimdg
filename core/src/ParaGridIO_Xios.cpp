/*!
 * @author  Tim Spain <timothy.spain@nersc.no>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 */

#include "include/ParaGridIO.hpp"

#include "include/CommonRestartMetadata.hpp"
#include "include/FileCallbackCloser.hpp"
#include "include/Finalizer.hpp"
#include "include/Logged.hpp"
#include "include/MissingData.hpp"
#include "include/ModelMPI.hpp"

#ifdef USE_XIOS
#include "include/Xios.hpp"
#endif
#include "include/gridNames.hpp"

#include <ncDim.h>
#include <ncException.h>
#include <ncFile.h>
#include <ncVar.h>

#include <algorithm>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <string>

namespace Nextsim {

ParaGridIO::ParaGridIO(ParametricGrid& grid)
    : IParaGridIO(grid)
{
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.setupDomains();
    xiosHandler.setupGrids();
}

ParaGridIO::~ParaGridIO() = default;

ModelState ParaGridIO::getModelState(const std::string& filePath)
{
    // Close the XIOS context definition, if it hasn't already been closed
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.close_context_definition();

    ModelMetadata& metadata = ModelMetadata::getInstance();
    if (metadata.initialFileName != filePath) {
        throw std::runtime_error("ParaGridIO::getModelState: file path '" + filePath
            + "' is inconsistent with model.init_file '" + metadata.initialFileName + "'");
    }

    // Get all variables in the file and load them into a new ModelState
    ModelState state;
    for (const std::string& fieldId : xiosHandler.inputRestartFieldNames) {
        const ModelArray::Type& type = xiosHandler.getFieldType(fieldId);
        if (type == ModelArray::Type::H) {
            HField field(ModelArray::Type::H);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else if (type == ModelArray::Type::VERTEX) {
            VertexField field(ModelArray::Type::VERTEX);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else if (type == ModelArray::Type::DG) {
            DGField field(ModelArray::Type::DG);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else if (type == ModelArray::Type::DGSTRESS) {
            DGSField field(ModelArray::Type::DGSTRESS);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else if (type == ModelArray::Type::CG) {
            CGField field(ModelArray::Type::CG);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else {
            throw std::runtime_error("ParaGridIO::getModelState: field type for field " + fieldId
                + " is not supported.");
        }
    }

    // Assume that all fields in the supplied ModelState are necessary, and so read them from file.
    for (auto& [fieldId, modelarray] : state.data) {
        const std::string inputFieldId = fieldId + "_input";
        if (xiosHandler.inputRestartFieldNames.count(fieldId) == 0) {
            Logged::warning(
                "ParaGridIO::getModelState: field " + fieldId + " is not configured as a restart.");
            continue;
        }
        xiosHandler.read(inputFieldId, modelarray);
    }
    return state;
}

ModelState ParaGridIO::readForcingTimeStatic(
    const std::set<std::string>& forcings, const TimePoint& time, const std::string& filePath)
{
    // Close the XIOS context definition, if it hasn't already been closed
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.close_context_definition();

    if (xiosHandler.forcingFilename != filePath) {
        throw std::runtime_error("ParaGridIO::readForcingTimeStatic: file path '" + filePath
            + "' is inconsistent with XiosForcing.filename '" + xiosHandler.forcingFilename + "'");
    }

    // Get all forcings and load them into a new ModelState
    ModelState state;
    for (const std::string& fieldId : forcings) {
        if (xiosHandler.forcingFieldNames.count(fieldId) == 0) {
            throw std::runtime_error("ParaGridIO::readForcingTimeStatic: field " + fieldId
                + " is not configured as a forcing.");
        }
        const std::string forcingFieldId = fieldId + "_forcing";
        const ModelArray::Type& type = xiosHandler.getFieldType(forcingFieldId);
        // ASSUME all forcings are HFields: finite volume fields on the same
        // grid as ice thickness
        if (type == ModelArray::Type::H) {
            HField field(ModelArray::Type::H);
            field.resize();
            state.merge(ModelState { { { fieldId, field } }, {} });
        } else {
            throw std::runtime_error("ParaGridIO::readForcingTimeStatic: field type for field "
                + fieldId + " is not supported.");
        }
    }

    // Read all forcings from file
    for (auto& [fieldId, modelarray] : state.data) {
        const std::string forcingFieldId = fieldId + "_forcing";
        if (forcings.count(fieldId)) {
            xiosHandler.read(forcingFieldId, modelarray);
        }
    }
    return state;
}

void ParaGridIO::dumpModelState(const ModelState& state, const std::string& filePath)
{
    // Close the XIOS context definition, if it hasn't already been closed
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.close_context_definition();

    ModelMetadata& metadata = ModelMetadata::getInstance();

    if (filePath.find(xiosHandler.outputFileId) == std::string::npos) {
        throw std::runtime_error("ParaGridIO::dumpModelState: file path '" + filePath
            + "' is inconsistent with model.restart_file '" + metadata.finalFileName + "'");
    }

    // Assume that all fields in the supplied ModelState are necessary, and so write them to file.
    for (const auto& [fieldId, modelarray] : state.data) {
        if (xiosHandler.outputRestartFieldNames.count(fieldId) == 0) {
            Logged::warning("ParaGridIO::dumpModelState: field " + fieldId
                + " is not configured as a restart.");
            continue;
        }
        xiosHandler.write(fieldId, modelarray);
    }
}

void ParaGridIO::writeDiagnosticTime(const ModelState& state, const std::string& filePath)
{
    // Close the XIOS context definition, if it hasn't already been closed
    Xios& xiosHandler = Xios::getInstance();
    xiosHandler.close_context_definition();

    if (filePath.find(xiosHandler.diagnosticFileId) == std::string::npos) {
        throw std::runtime_error("ParaGridIO::writeDiagnosticTime: file path '" + filePath
            + "' is inconsistent with XiosDiagnostic.filename '" + xiosHandler.diagnosticFilename
            + "'");
    }

    // Assume that all fields in the supplied ModelState are necessary, and so write them to file.
    for (const auto& [fieldId, modelarray] : state.data) {
        if (xiosHandler.diagnosticFieldNames.count(fieldId) == 0) {
            throw std::runtime_error("ParaGridIO::writeDiagnosticTime: field " + fieldId
                + " is not configured as a diagnostic.");
        }
        xiosHandler.write(fieldId, modelarray);
    }
}

} /* namespace Nextsim */
