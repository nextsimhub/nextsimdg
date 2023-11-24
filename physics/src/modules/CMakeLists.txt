set(AllModuleSources "")

# Get all module subdirectories. These must be suffix with "Module", and
# nothing else in this directory can be.
file(GLOB ModuleSubDirs LIST_DIRECTORIES TRUE "*Module")

foreach(ModuleSubDir IN ITEMS "${ModuleSubDirs}")
    IF (NOT EXISTS("include/${ModuleSubDir}.hpp"))
        # Run the header python script
        EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "${ModuleHeaderScript}" WORKING_DIRECTORY "${ModuleSubDir}")
    ENDIF()
    
    # Run the source python script
    EXECUTE_PROCESS(COMMAND "${PYTHON_EXECUTABLE}" "${ModuleBuilderScript}" WORKING_DIRECTORY "${ModuleSubDir}")
    
    # List all sources in the subdirectory
    file(GLOB ThisModuleSources LIST_DIRECTORIES FALSE "${ModuleSubDir}/*.cpp")
    
    set(AllModuleSources
        "${AllModuleSources}"
        "${ThisModuleSources}"
        )
endforeach()

set(NextsimSources
    "${NextsimSources}"
    "${AllModuleSources}"
    PARENT_SCOPE)
