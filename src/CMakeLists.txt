# Build the main neXtSIM model

set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Set the location where the module ipp files will be stored
set(ModuleLoaderIppTargetDirectory
"${CMAKE_CURRENT_SOURCE_DIR}/modules/")

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
    "ModuleLoader.cpp"
    "PrognosticData.cpp"
    "ExternalData.cpp"
    "PhysicsData.cpp"
    )

set(IceAlbedoModuleSources
    "CCSMIceAlbedo.cpp"
    "SMUIceAlbedo.cpp"
    "SMU2IceAlbedo.cpp"
    )

set(IceOceanHeatFluxModuleSources
    "BasicIceOceanHeatFlux.cpp"
    )

set(ColumnPhysicsModuleSources
    "NextsimPhysics.cpp"
    )

set(ThermodynamicsModuleSources
    "ThermoIce0.cpp"
    )

set(ConcentrationModelModuleSources
    "HiblerConcentration.cpp"
    )
    
set(ModuleSources
    "${IceAlbedoModuleSources}"
    "${IceOceanHeatFluxModuleSources}"
    "${ColumnPhysicsModuleSources}"
    "${ThermodynamicsModuleSources}"
    "${ConcentrationModelModuleSources}"
    )

set(BaseHeaders
    "${INC}/Logged.hpp"
    "${INC}/Timer.hpp"
    "${INC}/Model.hpp"
    "${INC}/Iterator.hpp"
    "${INC}/SimpleIterant.hpp"
    )

set(NextsimSources
    "${BaseSources}"
    "${ModuleSources}"
    )

add_executable(nextsim ${NextsimSources})

target_include_directories(nextsim PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${Boost_INCLUDE_DIRS}
    ${ModuleLoaderIppTargetDirectory}
    )
target_link_libraries(nextsim LINK_PUBLIC ${Boost_LIBRARIES})

# Build the module loader files. Defines the 'parse_modules' target
add_subdirectory(modules)

add_dependencies(nextsim parse_modules)
