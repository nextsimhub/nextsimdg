# Sources for the main neXtSIM model

#set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")

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
    )


#set(IceAlbedoModuleSources
#    "CCSMIceAlbedo.cpp"
#    "SMUIceAlbedo.cpp"
#    "SMU2IceAlbedo.cpp"
#    )

#set(IceOceanHeatFluxModuleSources
#    "BasicIceOceanHeatFlux.cpp"
#    )

#set(ColumnPhysicsModuleSources
#    "NextsimPhysics.cpp"
#    )

#set(ThermodynamicsModuleSources
#    "ThermoIce0.cpp"
#    )

#set(ConcentrationModelModuleSources
#    "HiblerConcentration.cpp"
#    )
    
#set(ModuleSources
#    "${IceAlbedoModuleSources}"
#    "${IceOceanHeatFluxModuleSources}"
#    "${ColumnPhysicsModuleSources}"
#    "${ThermodynamicsModuleSources}"
#    "${ConcentrationModelModuleSources}"
#    )

#set(BaseHeaders
#    "${INC}/Logged.hpp"
#    "${INC}/Timer.hpp"
#    "${INC}/Model.hpp"
#    "${INC}/Iterator.hpp"
#    "${INC}/SimpleIterant.hpp"
#    )

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(NextsimSources
    "${NextsimSources}"
    "${BaseSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
#add_executable(nextsim ${NextsimSources})

#target_include_directories(nextsim PRIVATE
#    ${CMAKE_CURRENT_SOURCE_DIR}
#    ${Boost_INCLUDE_DIRS}
#    ${ModuleLoaderIppTargetDirectory}
#    )
#target_link_libraries(nextsim LINK_PUBLIC ${Boost_LIBRARIES})

# Build the module loader files. Defines the 'parse_modules' target
add_subdirectory(modules)

#add_dependencies(nextsim parse_modules)
