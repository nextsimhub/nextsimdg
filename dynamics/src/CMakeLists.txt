set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")


set(IncludeDir "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Build the main neXtSIM model

set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")


### Set the list of all cpp files
FILE(GLOB Nextsim_Sources *.cpp)



set(NextsimSources
    "${NextsimSources}"
    "${Sources}"
    "${ModuleSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${IncludeDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
