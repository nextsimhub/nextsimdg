# Sources for the main neXtSIM model

set(BaseSources
    "Logged.cpp"
    "Timer.cpp"
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
    "Time.cpp"
    "${ModelArrayStructure}/ModelArrayDetails.cpp"
    )

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

add_subdirectory("${ModuleDir}")

set(StructureModuleSources
    "RectangularGrid.cpp"
    "ParametricGrid.cpp"
    "IStructureModule.cpp"
    )

#set(FreezingPointModuleSources
#    "IFreezingPointModule.cpp"
#    )
#
#set(DynamicsModuleSources
#    "DynamicsModule.cpp"
#    "Dynamics.cpp"
#    )
#

set(ModuleSources
    "${StructureModuleSources}"
    "${FreezingPointModuleSources}"
    "${DynamicsModuleSources}"
    )

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(NextsimSources
    "${BaseSources}"
    "${NextsimSources}"
    "${ModuleSources}"
#    "${DiagnosticOutputModuleSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${ModelArrayStructure}"
    PARENT_SCOPE)
