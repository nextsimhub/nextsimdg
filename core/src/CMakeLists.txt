# Sources for the main neXtSIM model

set(BaseSources
    "main.cpp"
    "Logged.cpp"
    "Timer.cpp"
    "Model.cpp"
    "Iterator.cpp"
    "SimpleIterant.cpp"
    "Configurator.cpp"
    "ConfiguredModule.cpp"
    "CommandLineParser.cpp"
    "ElementData.cpp"
    "PrognosticElementData.cpp"
    "ExternalData.cpp"
    "DevGridIO.cpp"
    "RectGridIO.cpp"
    "DevStep.cpp"
    "StructureFactory.cpp"
    "ModelArray.cpp"
    )

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(StructureModuleSources
    "DevGrid.cpp"
    "RectangularGrid.cpp"
    "IStructureModule.cpp"
    )

set(FreezingPointModuleSources
    "IFreezingPointModule.cpp"
    )

set(ModuleSources
    "${StructureModuleSources}"
    "${FreezingPointModuleSources}"
    )

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(NextsimSources
    "${NextsimSources}"
    "${BaseSources}"
    "${ModuleSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
