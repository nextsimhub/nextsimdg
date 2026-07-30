set(COMMON_INCLUDE_DIRS "../src/include" "../../core/src" "../../core/src/modules"
                        "../../core/src/include" "../../dynamics/src")

set(PHYSICS_INCLUDE_DIRS "../../physics/src" "../../physics/src/modules")

set(ModulesRoot "../src/modules")

# definitions for tests that do not link to nextsimlib
add_compile_definitions("NEXTSIM_FLOAT_TYPE=${FloatType}")

# Generate netCDF files for tests
execute_process(COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/alt_names.py"
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/ParaGridIO_input_test.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

include_directories(${COMMON_INCLUDE_DIRS})
set(MODEL_INCLUDE_DIR "../../core/src/discontinuousgalerkin")

if(ENABLE_MPI)

    # This function takes a list of test executables and adds them to the test suite
    function(add_mpi_test test_exec_list)
        foreach(test_exec IN LISTS test_exec_list)
            # Extract the number of MPI ranks from the test name (e.g., TestName_MPI3 -> 3)
            string(REGEX MATCH "_MPI([0-9]+)$" _match "${test_exec}")
            set(num_ranks "${CMAKE_MATCH_1}")
            if(num_ranks)
                # In the MPI tests it's important to include the "-ns" flag which tells doctest to
                # fail if test has been skipped. Without this, the test will technically pass if we
                # run with fewer ranks than requested.
                add_test(NAME ${test_exec} COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG}
                                                   ${num_ranks} ./${test_exec} -ns)
            else()
                message(WARNING "Could not extract MPI rank from test name: ${test_exec}")
            endif()
        endforeach()
    endfunction()

    # Generate mask files needed for MPI test from respective cdl files
    add_custom_command(
        OUTPUT rectgrid_test_mask.nc
        COMMAND ncgen -b -o rectgrid_test_mask.nc ${CMAKE_CURRENT_SOURCE_DIR}/rectgrid_test_mask.cdl
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/rectgrid_test_mask.cdl)
    add_custom_command(
        OUTPUT paragrid_test_mask.nc
        COMMAND
            ncgen -b -o paragrid_test_mask.nc ${CMAKE_CURRENT_SOURCE_DIR}/paragrid_test_mask.cdl
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/paragrid_test_mask.cdl
    )
    add_custom_command(
        OUTPUT halo_test_mask.nc
        COMMAND
            ncgen -b -o halo_test_mask.nc ${CMAKE_CURRENT_SOURCE_DIR}/halo_test_mask.cdl
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/halo_test_mask.cdl
    )
    add_custom_target(
        generate_masks
        ALL
        DEPENDS rectgrid_test_mask.nc paragrid_test_mask.nc halo_test_mask.nc
    )

    # Generate partition files using domain_decomp
    add_custom_command(
        OUTPUT rectgrid_test_partition_metadata_3.nc
        COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 decomp -m mask -o xy -g
                rectgrid_test_mask.nc --output-prefix rectgrid_test
        DEPENDS rectgrid_test_mask.nc)
    add_custom_command(
        OUTPUT paragrid_test_partition_metadata_2.nc
        COMMAND
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 decomp -m mask -o xy
                -g paragrid_test_mask.nc --output-prefix paragrid_test
        DEPENDS paragrid_test_mask.nc
    )
    add_custom_command(
        OUTPUT halo_cb_test_partition_metadata_3.nc
        COMMAND
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 decomp -m mask -o xy
                -g halo_test_mask.nc --output-prefix halo_cb_test
        DEPENDS halo_test_mask.nc
    )
    add_custom_command(
        OUTPUT halo_pb_test_partition_metadata_3.nc
        COMMAND
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 decomp -m mask -o xy
                -g halo_test_mask.nc --output-prefix halo_pb_test --px --py
        DEPENDS halo_test_mask.nc
    )
    add_custom_target(
        generate_partition_files
        ALL
        DEPENDS rectgrid_test_partition_metadata_3.nc
        paragrid_test_partition_metadata_2.nc
        halo_cb_test_partition_metadata_3.nc
        halo_pb_test_partition_metadata_3.nc
    )

    if(ENABLE_XIOS)
        set(XIOS_INCLUDE_LIST
            "${xios_INCLUDES}" "${xios_EXTERNS}/blitz/" "${xios_EXTERNS}/boost_extraction/include"
            "${xios_EXTERNS}/rapidxml/include")
        set(MODEL_INCLUDE_DIR "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}")

        # Generate XML file for XIOS configuration
        add_custom_command(
            OUTPUT iodef.xml
            COMMAND
                python3 ${CMAKE_CURRENT_SOURCE_DIR}/../src/generate_iodef.py ${DGComp}
                ${DGStressComp} ${CMAKE_CURRENT_SOURCE_DIR}/../src/iodef.xml.jinja
                ${CMAKE_CURRENT_BINARY_DIR}/iodef.xml
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/../src/iodef.xml.jinja)
        add_custom_target(generate_xios_iodef_file ALL DEPENDS iodef.xml)

        # Generate input and partition files needed for XIOS tests from respective cdl files
        add_custom_command(
            OUTPUT xios_test_input.nc
            COMMAND ncgen -k netCDF-4 -o xios_test_input.nc
                    ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_input.cdl
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_input.cdl)
        add_custom_target(generate_xios_input_file ALL DEPENDS xios_test_input.nc)
        add_custom_command(
            OUTPUT xios_test_era5_forcing.nc
            COMMAND ncgen -k netCDF-4 -o xios_test_era5_forcing.nc
                    ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_era5_forcing.cdl
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_era5_forcing.cdl
        )
        add_custom_target(generate_xios_era5_forcing_file ALL DEPENDS xios_test_era5_forcing.nc)
        add_custom_command(
            OUTPUT xios_test_topaz_forcing.nc
            COMMAND ncgen -k netCDF-4 -o xios_test_topaz_forcing.nc
                    ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_topaz_forcing.cdl
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_topaz_forcing.cdl
        )
        add_custom_target(generate_xios_topaz_forcing_file ALL DEPENDS xios_test_topaz_forcing.nc)
        add_custom_command(
            OUTPUT xios_test_mask.nc
            COMMAND ncgen -b -o xios_test_mask.nc ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_mask.cdl
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_mask.cdl)
        add_custom_target(generate_xios_masks ALL DEPENDS xios_test_mask.nc)
        add_custom_command(
            OUTPUT xios_test_partition_metadata_1.nc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 1 decomp -m mask -o xy -g
                    xios_test_mask.nc --output-prefix xios_test
            DEPENDS xios_test_mask.nc)
        add_custom_command(
            OUTPUT xios_test_partition_metadata_2.nc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 decomp -m mask -o xy -g
                    xios_test_mask.nc --output-prefix xios_test
            DEPENDS xios_test_mask.nc)
        add_custom_command(
            OUTPUT xios_test_partition_metadata_3.nc
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 3 decomp -m mask -o xy -g
                    xios_test_mask.nc --output-prefix xios_test
            DEPENDS xios_test_mask.nc)
        add_custom_target(
            generate_xios_partition_metadata ALL
            DEPENDS xios_test_partition_metadata_1.nc xios_test_partition_metadata_2.nc
                    xios_test_partition_metadata_3.nc)

        add_executable(testXiosCalendar_MPI1 "XiosCalendar_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosCalendar_MPI1 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosCalendar_MPI1 PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}"
                                          "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosCalendar_MPI1 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosAxis_MPI3 "XiosAxis_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosAxis_MPI3 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosAxis_MPI3 PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}"
                                      "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosAxis_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosReadForcing_MPI2 "XiosReadForcing_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosReadForcing_MPI2
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosReadForcing_MPI2
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosReadForcing_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosReadDiagnostic_MPI3 "XiosReadDiagnostic_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosReadDiagnostic_MPI3
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosReadDiagnostic_MPI3
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosReadDiagnostic_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosReadRestart_MPI2 "XiosReadRestart_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosReadRestart_MPI2
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosReadRestart_MPI2
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosReadRestart_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosWriteDiagnostic_MPI2 "XiosWriteDiagnostic_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosWriteDiagnostic_MPI2
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosWriteDiagnostic_MPI2
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosWriteDiagnostic_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosWriteRestart_MPI2 "XiosWriteRestart_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosWriteRestart_MPI2
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosWriteRestart_MPI2
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosWriteRestart_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosReadWriteRestart_MPI2 "XiosReadWriteRestart_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosReadWriteRestart_MPI2
                                   PRIVATE USE_XIOS TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testXiosReadWriteRestart_MPI2
            PRIVATE ${PHYSICS_INCLUDE_DIRS} "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}"
                    "${ModulesRoot}/StructureModule")
        target_link_libraries(testXiosReadWriteRestart_MPI2 PRIVATE nextsimlib doctest::doctest)

        set(XIOS_TESTS
            testXiosCalendar_MPI1
            testXiosAxis_MPI3
            testXiosReadForcing_MPI2
            testXiosReadDiagnostic_MPI3
            testXiosReadRestart_MPI2
            testXiosWriteDiagnostic_MPI2
            testXiosWriteRestart_MPI2
            testXiosReadWriteRestart_MPI2)

        add_mpi_test("${XIOS_TESTS}")

        # Ensure XiosWrite<Diagnostic/Restart>_test are run before XiosRead<Diagnostic/Restart>_test
        set_tests_properties(testXiosReadRestart_MPI2 PROPERTIES DEPENDS testXiosWriteRestart_MPI2)
        set_tests_properties(testXiosReadDiagnostic_MPI3 PROPERTIES DEPENDS
                                                                    testXiosWriteDiagnostic_MPI2)
    else()
        add_executable(testHaloExchangeCB_MPI3 "HaloExchangeCB_test.cpp" "MainMPI.cpp"
                                               "../src/ModelArray.cpp" "../src/ModelArraySlice.cpp")
        target_compile_definitions(testHaloExchangeCB_MPI3
                                   PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(
            testHaloExchangeCB_MPI3
            PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule" "../src"
                    "../../dynamics/src" "../../dynamics/src/include")
        target_link_libraries(testHaloExchangeCB_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testParaGrid_MPI2 "ParaGrid_test.cpp" "MainMPI.cpp")
        target_compile_definitions(
            testParaGrid_MPI2 PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
                                      TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\")
        target_include_directories(
            testParaGrid_MPI2 PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule"
                                      "../../dynamics/src/include")
        target_link_libraries(testParaGrid_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testRectGrid_MPI3 "RectGrid_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testRectGrid_MPI3
                                   PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(testRectGrid_MPI3 PRIVATE ${MODEL_INCLUDE_DIR}
                                                             "${ModulesRoot}/StructureModule")
        target_link_libraries(testRectGrid_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testModelMetadataCB_MPI3 "ModelMetadataCB_test.cpp" "MainMPI.cpp")
        target_compile_definitions(
            testModelMetadataCB_MPI3 PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
                                             TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\")
        target_include_directories(testModelMetadataCB_MPI3
                                   PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule")
        target_link_libraries(testModelMetadataCB_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testModelMetadataPB_MPI3 "ModelMetadataPB_test.cpp" "MainMPI.cpp")
        target_compile_definitions(
            testModelMetadataPB_MPI3 PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
                                             TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\")
        target_include_directories(testModelMetadataPB_MPI3
                                   PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule")
        target_link_libraries(testModelMetadataPB_MPI3 PRIVATE nextsimlib doctest::doctest)

        add_executable(testConfigOutput_MPI2 "ConfigOutput_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testConfigOutput_MPI2
                                   PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
        target_include_directories(testConfigOutput_MPI2 PRIVATE ${MODEL_INCLUDE_DIR}
                                                                 "../../dynamics/src/include")
        target_link_libraries(testConfigOutput_MPI2 PRIVATE nextsimlib doctest::doctest)

        set(MPI_TESTS
            testHaloExchangeCB_MPI3
            testParaGrid_MPI2
            testRectGrid_MPI3
            testModelMetadataCB_MPI3
            testModelMetadataPB_MPI3
            testConfigOutput_MPI2)
        add_mpi_test("${MPI_TESTS}")
    endif()

else()
    add_executable(testRectGrid "RectGrid_test.cpp")
    target_compile_definitions(testRectGrid PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
    target_include_directories(testRectGrid PRIVATE ${MODEL_INCLUDE_DIR}
                                                    "${ModulesRoot}/StructureModule")
    target_link_libraries(testRectGrid PRIVATE nextsimlib doctest::doctest)

    add_executable(testParaGrid "ParaGrid_test.cpp")
    target_compile_definitions(
        testParaGrid PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
                             TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\")
    target_include_directories(testParaGrid PRIVATE ${MODEL_INCLUDE_DIR}
                                                    "${ModulesRoot}/StructureModule")
    target_link_libraries(testParaGrid PRIVATE nextsimlib doctest::doctest)

    add_executable(testConfigOutput "ConfigOutput_test.cpp")
    target_compile_definitions(testConfigOutput
                               PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
    target_include_directories(testConfigOutput PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testConfigOutput PRIVATE nextsimlib doctest::doctest)

    add_executable(testIterator "Iterator_test.cpp")
    target_link_libraries(testIterator PRIVATE nextsimlib doctest::doctest)

    add_executable(testCommandLineParser "CommandLineParser_test.cpp" "ArgV.cpp")
    target_link_libraries(testCommandLineParser PRIVATE nextsimlib doctest::doctest)

    add_executable(testConfigurator "Configurator_test.cpp" "ArgV.cpp")
    target_link_libraries(testConfigurator PRIVATE nextsimlib doctest::doctest)

    add_executable(testConfiguredModule "ConfiguredModule_test.cpp" "ArgV.cpp")
    target_link_libraries(testConfiguredModule PRIVATE nextsimlib doctest::doctest)

    add_executable(testTimeClasses "Time_test.cpp")
    target_link_libraries(testTimeClasses PRIVATE nextsimlib doctest::doctest)

    add_executable(testModelComponent "ModelComponent_test.cpp")
    target_include_directories(testModelComponent PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testModelComponent PRIVATE nextsimlib doctest::doctest)

    add_executable(testModelArrayAccessor "ModelArrayAccessor_test.cpp")
    target_include_directories(testModelArrayAccessor PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testModelArrayAccessor PRIVATE nextsimlib doctest::doctest)

    # PrognosticData (and hopefully that alone) requires code from the physics tree
    add_executable(testPrognosticData "PrognosticData_test.cpp" "DynamicsModuleForPDtest.cpp"
                                      "MainTest.cpp")
    target_include_directories(
        testPrognosticData
        PRIVATE ${PHYSICS_INCLUDE_DIRS} ${MODEL_INCLUDE_DIR} "${ModulesRoot}/FreezingPointModule"
                "${ModulesRoot}/DynamicsModule")
    target_link_libraries(testPrognosticData PRIVATE nextsimlib doctest::doctest)

    # TODO: itegrate with the above test.
    add_executable(testPrognosticDataIO "PrognosticDataIO_test.cpp" "DynamicsModuleForPDtest.cpp")
    target_compile_definitions(testPrognosticDataIO PRIVATE ISDG=${isDG})
    target_include_directories(
        testPrognosticDataIO
        PRIVATE ${PHYSICS_INCLUDE_DIRS} ${MODEL_INCLUDE_DIR} "${ModulesRoot}/FreezingPointModule"
                "${ModulesRoot}/DynamicsModule" "${ModulesRoot}/StructureModule")
    target_link_libraries(testPrognosticDataIO PRIVATE nextsimlib doctest::doctest)

    add_executable(testMonthlyCubicBSpline "MonthlyCubicBSpline_test.cpp")
    target_link_libraries(testMonthlyCubicBSpline LINK_PUBLIC Boost::boost doctest::doctest)

    add_executable(testFileCallbackCloser "FileCallbackCloser_test.cpp"
                                          "../../core/src/FileCallbackCloser.cpp")
    target_link_libraries(testFileCallbackCloser PRIVATE doctest::doctest)

    add_executable(testFinalizer "Finalizer_test.cpp" "../../core/src/Finalizer.cpp")
    target_link_libraries(testFinalizer PRIVATE doctest::doctest)

    set(MODEL_INCLUDE_DIR "./testmodelarraydetails")

    add_executable(
        testModelArray "ModelArray_test.cpp" "../src/ModelArray.cpp" "../src/ModelArraySlice.cpp"
                       "${MODEL_INCLUDE_DIR}/ModelArrayDetails.cpp")
    target_include_directories(testModelArray PRIVATE "../src" "${MODEL_INCLUDE_DIR}")
    target_link_libraries(testModelArray PRIVATE doctest::doctest Eigen3::Eigen)

    add_executable(testSlice "Slice_test.cpp")
    target_include_directories(testSlice PRIVATE "../src")
    target_link_libraries(testSlice PRIVATE doctest::doctest)

    add_executable(
        testModelArraySlice
        "ModelArraySlice_test.cpp" "../src/ModelArray.cpp" "../src/ModelArraySlice.cpp"
        "${MODEL_INCLUDE_DIR}/ModelArrayDetails.cpp")
    target_include_directories(testModelArraySlice PRIVATE "../src" "${MODEL_INCLUDE_DIR}")
    target_link_libraries(testModelArraySlice PRIVATE doctest::doctest Eigen3::Eigen)

    set(TEST_EXECUTABLES
        testRectGrid
        testParaGrid
        testConfigOutput
        testIterator
        testCommandLineParser
        testConfigurator
        testConfiguredModule
        testTimeClasses
        testModelComponent
        testModelArrayAccessor
        testPrognosticData
        testPrognosticDataIO
        testMonthlyCubicBSpline
        testFileCallbackCloser
        testFinalizer
        testModelArray
        testSlice
        testModelArraySlice)

    foreach(test_exec IN LISTS TEST_EXECUTABLES)
        add_test(NAME ${test_exec} COMMAND ${test_exec})
    endforeach()

    # The help config test needs access to the full binary
    file(CREATE_LINK ${NEXTSIM_BINARY_PATH} ${CMAKE_CURRENT_BINARY_DIR}/nextsim SYMBOLIC)
endif()
