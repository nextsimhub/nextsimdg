# Sources for the main neXtSIM model

if(ENABLE_XIOS)
    set(XIOS_INCLUDE_LIST
        "${xios_INCLUDES}"
        "${xios_EXTERNS}/blitz/"
        "${xios_EXTERNS}/rapidxml/include"
    )
    set(NextsimIncludeDirs "${NextsimIncludeDirs}" "${XIOS_INCLUDE_LIST}")
    set(ParaGridIOImplSource "${CMAKE_CURRENT_SOURCE_DIR}/ParaGridIO_Xios.cpp")
else()
    set(ParaGridIOImplSource "${CMAKE_CURRENT_SOURCE_DIR}/ParaGridIO.cpp")
endif()

# Create a stable generated include dir and copy the chosen ModelArrayDetails.hpp
# This ensures a single canonical header, which can be included as
# "include/ModelArrayDetails.hpp" from any target.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/${ModelArrayStructure}/include/ModelArrayDetails.hpp"
    "${CMAKE_CURRENT_BINARY_DIR}/include/ModelArrayDetails.hpp"
    COPYONLY
)
set(NextsimIncludeDirs "${NextsimIncludeDirs}" "${CMAKE_CURRENT_BINARY_DIR}/include")

set(BaseSources
    "ModelMPI.cpp"
    "Logged.cpp"
    "Model.cpp"
    "Iterator.cpp"
    "Configurator.cpp"
    "ConfigurationHelpPrinter.cpp"
    "ConfiguredModule.cpp"
    "CommandLineParser.cpp"
    "CommonRestartMetadata.cpp"
    "RectGridIO.cpp"
    "FieldAdvection.cpp"
    "FileCallbackCloser.cpp"
    "Finalizer.cpp"
    "DevStep.cpp"
    "StructureFactory.cpp"
    "ModelArray.cpp"
    "ModelArraySlice.cpp"
    "ModelComponent.cpp"
    "ModelMetadata.cpp"
    "NetcdfMetadataConfiguration.cpp"
    "PrognosticData.cpp"
    "PDWriter.cpp"
    "Time.cpp"
    "${ModelArrayStructure}/ModelArrayDetails.cpp"
    "Xios.cpp"
)

set(ParallelNetCDFSources "${CMAKE_CURRENT_SOURCE_DIR}/ParallelNetcdfFile.cpp")

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

add_subdirectory("${ModuleDir}")

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(NextsimSourcesTemp
    "${BaseSources}"
    "${NextsimSources}"
    "${ModuleSources}"
    "${ParaGridIOImplSource}"
)

if(ENABLE_MPI)
    list(APPEND NextsimSourcesTemp "${ParallelNetCDFSources}")
endif()

set(NextsimSources "${NextsimSourcesTemp}" PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${ModelArrayStructure}"
    "${CMAKE_CURRENT_BINARY_DIR}/include"
    PARENT_SCOPE
)
