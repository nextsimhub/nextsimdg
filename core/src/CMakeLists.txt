# Sources for the main neXtSIM model

if(ENABLE_XIOS)
    set(XIOS_INCLUDE_LIST
        "${xios_INCLUDES}"
        "${xios_EXTERNS}/blitz/"
        "${xios_EXTERNS}/rapidxml/include"
    )
    set(NextsimIncludeDirs "${NextsimIncludeDirs}" "${XIOS_INCLUDE_LIST}")
    set(ParaGridIOImplSource "${CMAKE_CURRENT_SOURCE_DIR}/ParaGridIO_Xios.cpp")
else()
    set(ParaGridIOImplSource "${CMAKE_CURRENT_SOURCE_DIR}/ParaGridIO.cpp")
endif()

set(BaseSources
    "ModelMPI.cpp"
    "Logged.cpp"
    "Model.cpp"
    "Iterator.cpp"
    "Configurator.cpp"
    "ConfigurationHelpPrinter.cpp"
    "ConfiguredModule.cpp"
    "CommandLineParser.cpp"
    "CommonRestartMetadata.cpp"
    "RectGridIO.cpp"
    "FieldAdvection.cpp"
    "FileCallbackCloser.cpp"
    "Finalizer.cpp"
    "DevStep.cpp"
    "StructureFactory.cpp"
    "ModelArray.cpp"
    "ModelArraySlice.cpp"
    "ModelArrayStore.cpp"
    "ModelComponent.cpp"
    "ModelMetadata.cpp"
    "NetcdfMetadataConfiguration.cpp"
    "PrognosticData.cpp"
    "PDWriter.cpp"
    "Time.cpp"
    "${ModelArrayStructure}/ModelArrayDetails.cpp"
    "Xios.cpp"
)

set(ParallelNetCDFSources "${CMAKE_CURRENT_SOURCE_DIR}/ParallelNetcdfFile.cpp")

list(TRANSFORM BaseSources PREPEND "${CMAKE_CURRENT_SOURCE_DIR}/")

set(ModuleDir "${CMAKE_CURRENT_SOURCE_DIR}/modules")

add_subdirectory("${ModuleDir}")

list(TRANSFORM ModuleSources PREPEND "${ModuleDir}/")

set(NextsimSourcesTemp
    "${BaseSources}"
    "${NextsimSources}"
    "${ModuleSources}"
    "${ParaGridIOImplSource}"
)

if(ENABLE_MPI)
    list(APPEND NextsimSourcesTemp "${ParallelNetCDFSources}")
endif()

if(WITH_KOKKOS)
    list(APPEND NextsimSourcesTemp "${CMAKE_CURRENT_SOURCE_DIR}/kokkos/KokkosUtils.cpp")
endif()

set(NextsimSources "${NextsimSourcesTemp}" PARENT_SCOPE)

set(NextsimIncludeDirs
    "${NextsimIncludeDirs}"
    "${ModuleDir}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/${ModelArrayStructure}"
    PARENT_SCOPE
)
