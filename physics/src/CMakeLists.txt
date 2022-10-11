set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(IceAlbedoModuleSources
    "IIceAlbedoModule.cpp"
    "CCSMIceAlbedo.cpp"
    "SMUIceAlbedo.cpp"
    "SMU2IceAlbedo.cpp"
    "MonthlyIceAlbedo.cpp"
    )

set(IceOceanHeatFluxModuleSources
    "IceOceanHeatFluxModule.cpp"
    "BasicIceOceanHeatFlux.cpp"
    
    )

set(IceThermodynamicsModuleSources
    "IceThermodynamicsModule.cpp"
    "ThermoIce0.cpp"
    )
    
set(LateralIceSpreadModuleSources
    "LateralIceSpreadModule.cpp"
    "HiblerSpread.cpp"
    )

set(FluxCalculationModuleSources
    "FluxCalculationModule.cpp"
    "FiniteElementFluxes.cpp"
    "FiniteElementSpecHum.cpp"
    "MonthlyFluxes.cpp"
    )

set(AtmosphereBoundaryModuleSources
    "AtmosphereBoundaryModule.cpp"
    "ConstantAtmosphereBoundary.cpp"
    "ConfiguredAtmosphere.cpp"
    "FluxConfiguredAtmosphere.cpp"
    )

set(OceanBoundaryModuleSources
    "OceanBoundaryModule.cpp"
    "ConstantOceanBoundary.cpp"
    "ConfiguredOcean.cpp"
    "FluxConfiguredOcean.cpp"
    )

set(ModuleSources
    "${IceAlbedoModuleSources}"
    "${IceOceanHeatFluxEDModuleSources}"
    "${IceOceanHeatFluxModuleSources}"
    "${ColumnPhysicsModuleSources}"
    "${ThermodynamicsModuleSources}"
    "${IceThermodynamicsModuleSources}"
    "${LateralIceSpreadModuleSources}"
    "${FluxCalculationModuleSources}"
    "${StateModuleSources}"
    "${AtmosphereBoundaryModuleSources}"
    "${OceanBoundaryModuleSources}"
    )

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(Sources
    "IceGrowth.cpp"
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
