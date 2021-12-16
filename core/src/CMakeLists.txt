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
    "ModuleLoader.cpp"
    "PrognosticData.cpp"
    "ExternalData.cpp"
    )

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(NextsimSources
    "${NextsimSources}"
    "${BaseSources}"
    PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    PARENT_SCOPE)

# Build the module loader files. Defines the 'parse_modules' target
add_subdirectory(modules)
