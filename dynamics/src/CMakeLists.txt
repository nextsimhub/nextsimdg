set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(IncludeDir "${CMAKE_CURRENT_SOURCE_DIR}/include")


### Set the list of all cpp files
FILE(GLOB DynamicsSources *.cpp)
FILE(GLOB GPUSources gpu/*.cpp)

set(NextsimSources
    "${NextsimSources}"
    "${Sources}"
    "${ModuleSources}"
    "${DynamicsSources}"
    "${GPUSources}"
    PARENT_SCOPE)



set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${IncludeDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)
