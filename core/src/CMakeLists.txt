# Sources for the main neXtSIM model

set(BaseSources
    "Logged.cpp"
    "Timer.cpp"
    "ModelConfig.cpp"
    "Model.cpp"
    "Iterator.cpp"
    "Configurator.cpp"
    "ConfigurationHelpPrinter.cpp"
    "ConfiguredModule.cpp"
    "CommandLineParser.cpp"
    "CommonRestartMetadata.cpp"
    "RectGridIO.cpp"
    "ParaGridIO.cpp"
    "FileCallbackCloser.cpp"
    "DevStep.cpp"
    "StructureFactory.cpp"
    "MissingData.cpp"
    "ModelArray.cpp"
    "ModelComponent.cpp"
    "ModelMetadata.cpp"
    "NetcdfMetadataConfiguration.cpp"
    "NZLevels.cpp"
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

set(NextsimSourcesTemp "${BaseSources}" "${NextsimSources}" "${ModuleSources}")

if(ENABLE_MPI)
    list(APPEND NextsimSourcesTemp "${ParallelNetCDFSources}")
endif()

set(NextsimSources "${NextsimSourcesTemp}" PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${ModelArrayStructure}"
    PARENT_SCOPE
)
