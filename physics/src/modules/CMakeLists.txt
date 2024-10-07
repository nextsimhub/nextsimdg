set(AllModuleSources "")

# Get all module subdirectories. These must be suffix with "Module", and
# nothing else in this directory can be.
file(GLOB ModuleSubDirs LIST_DIRECTORIES TRUE "*Module")

foreach(ModuleSubDir IN LISTS ModuleSubDirs)
    # Ignore directories without a module.cfg file
    if(NOT EXISTS "${ModuleSubDir}/module.cfg")
        continue()
    endif()

    if(NOT EXISTS ("include/${ModuleSubDir}.hpp"))
        # Run the header python script
        execute_process(
            COMMAND "${Python_EXECUTABLE}" "${ModuleHeaderScript}"
            WORKING_DIRECTORY "${ModuleSubDir}"
        )
    endif()

    # Run the source python script
    execute_process(
        COMMAND "${Python_EXECUTABLE}" "${ModuleBuilderScript}"
        WORKING_DIRECTORY "${ModuleSubDir}"
    )

    # List all sources in the subdirectory
    file(GLOB ThisModuleSources LIST_DIRECTORIES FALSE "${ModuleSubDir}/*.cpp")

    set(AllModuleSources "${AllModuleSources}" "${ThisModuleSources}")
endforeach()

set(NextsimSources "${NextsimSources}" "${AllModuleSources}" PARENT_SCOPE)
