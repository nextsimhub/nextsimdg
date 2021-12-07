find_package( PythonInterp 3 REQUIRED )

set(ModuleLoaderIncludes
    "moduleLoaderHeaders.ipp"
    "moduleLoaderFunctions.ipp"
    "moduleLoaderNames.ipp"
    "moduleLoaderAssignments.ipp"
)

# Modules for the model infrastructure are defined in this directory
set(ModuleLoaderFileCommon ${CMAKE_CURRENT_SOURCE_DIR}/modules.json)

# Other modules are defined in component subdirectories
set(ModuleSubdirectories
"${ModularComponents}"
)

# Transform the list of component subdirectories to a list of module file paths
# Each component that provides modules should have the directories 
# ${PROJECT_SOURCE_DIR}/{component}/src/modules (for source files)
# ${PROJECT_SOURCE_DIR}/{component}/src/modules/include (for header files)
list(TRANSFORM ModuleSubdirectories PREPEND ${PROJECT_SOURCE_DIR})
list(TRANSFORM ModuleSubdirectories APPEND "/src/modules")


# Locations of the module loader files
set(ModuleLoaderFiles
${CMAKE_CURRENT_SOURCE_DIR}
${ModuleSubdirectories}
)
# And the files themselves
list(TRANSFORM ModuleLoaderFiles APPEND "/modules.json")

add_custom_target(
parse_modules ALL
COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/moduleloader_builder.py "--ipp" ${ModuleLoaderIppTargetDirectory} ${ModuleLoaderFiles}
BYPRODUCTS ${ModuleLoaderIncludes}
COMMENT "Generating inclusion files for ModuleLoader.cpp"
)
