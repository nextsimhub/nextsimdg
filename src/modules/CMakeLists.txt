find_package( PythonInterp 3 REQUIRED )

set(ModuleLoaderIncludes
    "moduleLoaderHeaders.ipp"
    "moduleLoaderFunctions.ipp"
    "moduleLoaderNames.ipp"
    "moduleLoaderAssignments.ipp"
)

set(ModuleLoaderFile
${CMAKE_CURRENT_SOURCE_DIR}/modules.json)

add_custom_target(
parse_modules ALL
COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/moduleloader_builder.py "--ipp" ${ModuleLoaderIppTargetDirectory} ${ModuleLoaderFile}
BYPRODUCTS ${ModuleLoaderIncludes}
COMMENT "Generating inclusion files for ModuleLoader.cpp"
)
