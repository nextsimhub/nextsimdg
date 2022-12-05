set(CoreSrc "../core/src")
set(CoreTest "../core/test")
set(CoreInclude "${CoreSrc}/include")
set(CoreModulesDir "${CoreSrc}/modules")
set(ModelArrayDetails "testmodelarraydetails")
set(PhysicsSrc "../physics/src")
set(PhysicsModulesDir "${PhysicsSrc}/modules")

set(NextsimSources
    "${NextsimSources}"
    PARENT_SCOPE)

#Remove main.cpp from nextsim sources to create a generic list of sources for integration testing
set(mainSrc "${CMAKE_SOURCE_DIR}/core/src/main.cpp;")
string(REPLACE "${mainSrc}" "" TestSources "${NextsimSources}")

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    PARENT_SCOPE)

add_executable(testSimpleConfigIntegration
    "simple_config_integration_test.cpp"
    "${CoreTest}/ArgV.cpp"
    "${TestSources}")

add_definitions(-DBOOST_ALL_DYN_LINK)

target_include_directories(testSimpleConfigIntegration PRIVATE "${CoreSrc}" "${CoreInclude}" "${CoreSrc}/${ModelArrayStructure}" "${PhysicsSrc}" "${PhysicsModulesDir}" "${CoreModulesDir}" "${CoreTest}" "${netCDF_INCLUDE_DIR}")
target_link_directories(testSimpleConfigIntegration PUBLIC "${netCDF_LIB_DIR}") 
target_link_libraries(testSimpleConfigIntegration LINK_PUBLIC Boost::program_options Boost::log Catch2::Catch2 Eigen3::Eigen "${NSDG_NetCDF_Library}")



