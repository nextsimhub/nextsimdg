set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(IceAlbedoModuleSources
    "IIceAlbedoModule.cpp"
    "CCSMIceAlbedo.cpp"
    "SMUIceAlbedo.cpp"
    "SMU2IceAlbedo.cpp"
    )

set(IceOceanHeatFluxModuleSources
    "IIceOceanHeatFluxModule.cpp"
    "BasicIceOceanHeatFlux.cpp"
    
    )

set(ColumnPhysicsModuleSources
    "IPhysics1dModule.cpp"
    "NextsimPhysics.cpp"
    )

set(ThermodynamicsModuleSources
    "IThermodynamicsModule.cpp"
    "ThermoIce0.cpp"
    )

set(ConcentrationModelModuleSources
    "IConcentrationModelModule.cpp"
    "HiblerConcentration.cpp"
    )
    
set(ModuleSources
    "${IceAlbedoModuleSources}"
    "${IceOceanHeatFluxModuleSources}"
    "${ColumnPhysicsModuleSources}"
    "${ThermodynamicsModuleSources}"
    "${ConcentrationModelModuleSources}"
    )

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(Sources
    ""
    )

list(TRANSFORM Sources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(NextsimSources
    "${NextsimSources}"
    "${Sources}"
    "${ModuleSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
