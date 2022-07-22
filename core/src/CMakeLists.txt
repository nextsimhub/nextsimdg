# Sources for the main neXtSIM model

set(BaseSources
    "main.cpp"
    "Logged.cpp"
    "Timer.cpp"
    "Model.cpp"
    "Iterator.cpp"
    "Configurator.cpp"
    "ConfiguredModule.cpp"
    "CommandLineParser.cpp"
    "CommonRestartMetadata.cpp"
    "DevGridIO.cpp"
    "RectGridIO.cpp"
    "DevStep.cpp"
    "StructureFactory.cpp"
    "MissingData.cpp"
    "ModelArray.cpp"
    "ModelComponent.cpp"
    "ModelData.cpp"
    "ModelMetadata.cpp"
    "PrognosticData.cpp"
    "Time.cpp"
    )

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

# Coupling (or not) source files
set("CouplingSources" "")
if(USE_OASIS)
set(CouplingSources
    "OasisCoupler.cpp"
    "OasisData.cpp"
    )
list(TRANSFORM CouplingSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

#elseif(USE_SOMETHING_ELSE) For future coupling systems
else()
set(CouplingSources
    "NoCoupler.cpp"
    )
list(TRANSFORM CouplingSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")
endif()

set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(StructureModuleSources
    "DevGrid.cpp"
    "RectangularGrid.cpp"
    "IStructureModule.cpp"
    )

set(FreezingPointModuleSources
    "IFreezingPointModule.cpp"
    )

set(DiagnosticOutputModuleSources
    "DiagnosticOutputModule.cpp"
    "SimpleOutput.cpp"
    )

set(ModuleSources
    "${StructureModuleSources}"
    "${FreezingPointModuleSources}"
    "${DiagnosticOutputModuleSources}"
    )

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(NextsimSources
    "${NextsimSources}"
    "${BaseSources}"
    "${ModuleSources}"
    "${CouplingSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
