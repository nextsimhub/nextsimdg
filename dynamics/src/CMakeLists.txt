set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

set(IncludeDir "${CMAKE_CURRENT_SOURCE_DIR}/include")

### Set the list of all cpp files
file(GLOB DynamicsSources *.cpp)

set(NextsimSources
    "${NextsimSources}"
    "${Sources}"
    "${ModuleSources}"
    "${DynamicsSources}"
    PARENT_SCOPE
)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${IncludeDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE
)
