set(COMMON_INCLUDE_DIRS
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "../../core/src/include"
)

set(PHYSICS_INCLUDE_DIRS "../../physics/src" "../../physics/src/modules")

set(ModulesRoot "../src/modules")

# Generate netCDF files for tests
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/old_names.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)

include_directories(${COMMON_INCLUDE_DIRS})
set(MODEL_INCLUDE_DIR "../../core/src/discontinuousgalerkin")

if(ENABLE_MPI)
    # Generate partition files needed for MPI test from respective cdl files
    add_custom_command(
        OUTPUT partition_metadata_3.nc
        COMMAND
            ncgen -b -o partition_metadata_3.nc ${CMAKE_CURRENT_SOURCE_DIR}/partition_metadata_3.cdl
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/partition_metadata_3.cdl
    )
    add_custom_command(
        OUTPUT partition_metadata_2.nc
        COMMAND
            ncgen -b -o partition_metadata_2.nc ${CMAKE_CURRENT_SOURCE_DIR}/partition_metadata_2.cdl
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/partition_metadata_2.cdl
    )
    add_custom_target(
        generate_partition_files
        ALL
        DEPENDS partition_metadata_3.nc partition_metadata_2.nc
    )

    add_executable(testRectGrid_MPI3 "RectGrid_test.cpp" "MainMPI.cpp")
    target_compile_definitions(
        testRectGrid_MPI3
        PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
    )
    target_include_directories(
        testRectGrid_MPI3
        PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule"
    )
    target_link_libraries(testRectGrid_MPI3 PRIVATE nextsimlib doctest::doctest)

    add_executable(testParaGrid_MPI2 "ParaGrid_test.cpp" "MainMPI.cpp")
    target_compile_definitions(
        testParaGrid_MPI2
        PRIVATE
            USE_MPI
            TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
            TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\"
    )
    target_include_directories(
        testParaGrid_MPI2
        PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule"
    )
    target_link_libraries(testParaGrid_MPI2 PRIVATE nextsimlib doctest::doctest)

    add_executable(testConfigOutput_MPI2 "ConfigOutput_test.cpp" "MainMPI.cpp")
    target_compile_definitions(
        testConfigOutput_MPI2
        PRIVATE USE_MPI TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
    )
    target_include_directories(testConfigOutput_MPI2 PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testConfigOutput_MPI2 PRIVATE nextsimlib doctest::doctest)

    if(ENABLE_XIOS)
        file(
            CREATE_LINK
                "${CMAKE_SOURCE_DIR}/core/test/iodef.xml"
                "${CMAKE_CURRENT_BINARY_DIR}/iodef.xml"
            SYMBOLIC
        )
        add_custom_command(
            OUTPUT xios_test_input.nc
            COMMAND ncgen -o xios_test_input.nc ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_input.cdl
            DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/xios_test_input.cdl
        )
        add_custom_target(generate_xios_input_file ALL DEPENDS xios_test_input.nc)
        set(XIOS_INCLUDE_LIST
            "${xios_INCLUDES}"
            "${xios_EXTERNS}/blitz/"
            "${xios_EXTERNS}/rapidxml/include"
        )

        add_executable(testXiosCalendar_MPI2 "XiosCalendar_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosCalendar_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosCalendar_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosCalendar_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosAxis_MPI2 "XiosAxis_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosAxis_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosAxis_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosAxis_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosDomain_MPI2 "XiosDomain_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosDomain_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosDomain_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosDomain_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosGrid_MPI2 "XiosGrid_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosGrid_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosGrid_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosGrid_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosField_MPI2 "XiosField_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosField_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosField_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosField_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosFile_MPI2 "XiosFile_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosFile_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosFile_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosFile_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosWrite_MPI2 "XiosWrite_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosWrite_MPI2 PRIVATE USE_XIOS)
        target_include_directories(
            testXiosWrite_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosWrite_MPI2 PRIVATE nextsimlib doctest::doctest)

        add_executable(testXiosRead_MPI2 "XiosRead_test.cpp" "MainMPI.cpp")
        target_compile_definitions(testXiosRead_MPI2 PRIVATE USE_XIOS USE_MPI)
        target_include_directories(
            testXiosRead_MPI2
            PRIVATE "${MODEL_INCLUDE_DIR}" "${XIOS_INCLUDE_LIST}" "${ModulesRoot}/StructureModule"
        )
        target_link_libraries(testXiosRead_MPI2 PRIVATE nextsimlib doctest::doctest)
    endif()
else()
    add_executable(testRectGrid "RectGrid_test.cpp")
    target_compile_definitions(testRectGrid PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\")
    target_include_directories(
        testRectGrid
        PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule"
    )
    target_link_libraries(testRectGrid PRIVATE nextsimlib doctest::doctest)

    add_executable(testParaGrid "ParaGrid_test.cpp")
    target_compile_definitions(
        testParaGrid
        PRIVATE
            TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
            TEST_FILE_SOURCE=\"${CMAKE_CURRENT_SOURCE_DIR}\"
    )
    target_include_directories(
        testParaGrid
        PRIVATE ${MODEL_INCLUDE_DIR} "${ModulesRoot}/StructureModule"
    )
    target_link_libraries(testParaGrid PRIVATE nextsimlib doctest::doctest)

    add_executable(testConfigOutput "ConfigOutput_test.cpp")
    target_compile_definitions(
        testConfigOutput
        PRIVATE TEST_FILES_DIR=\"${CMAKE_CURRENT_BINARY_DIR}\"
    )
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

    add_executable(testTimer "Timer_test.cpp")
    target_link_libraries(testTimer PRIVATE nextsimlib doctest::doctest)

    add_executable(testScopedTimer "ScopedTimer_test.cpp" "../src/ScopedTimer.cpp")
    target_link_libraries(testScopedTimer PRIVATE nextsimlib doctest::doctest)

    add_executable(testTimeClasses "Time_test.cpp")
    target_link_libraries(testTimeClasses PRIVATE nextsimlib doctest::doctest)

    add_executable(testModelComponent "ModelComponent_test.cpp")
    target_include_directories(testModelComponent PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testModelComponent PRIVATE nextsimlib doctest::doctest)

    add_executable(testModelArrayRef "ModelArrayRef_test.cpp")
    target_include_directories(testModelArrayRef PRIVATE ${MODEL_INCLUDE_DIR})
    target_link_libraries(testModelArrayRef PRIVATE nextsimlib doctest::doctest)

    # PrognosticData (and hopefully that alone) requires code from the physics tree
    add_executable(testPrognosticData "PrognosticData_test.cpp" "DynamicsModuleForPDtest.cpp")
    target_include_directories(
        testPrognosticData
        PRIVATE
            ${PHYSICS_INCLUDE_DIRS}
            ${MODEL_INCLUDE_DIR}
            "${ModulesRoot}/FreezingPointModule"
            "${ModulesRoot}/DynamicsModule"
    )
    target_link_libraries(testPrognosticData PRIVATE nextsimlib doctest::doctest)

    # TODO: itegrate with the above test.
    add_executable(testPrognosticDataIO "PrognosticDataIO_test.cpp" "DynamicsModuleForPDtest.cpp")
    target_compile_definitions(testPrognosticDataIO PRIVATE ISDG=${isDG})
    target_include_directories(
        testPrognosticDataIO
        PRIVATE
            ${PHYSICS_INCLUDE_DIRS}
            ${MODEL_INCLUDE_DIR}
            "${ModulesRoot}/FreezingPointModule"
            "${ModulesRoot}/DynamicsModule"
            "${ModulesRoot}/StructureModule"
    )
    target_link_libraries(testPrognosticDataIO PRIVATE nextsimlib doctest::doctest)

    add_executable(testMonthlyCubicBSpline "MonthlyCubicBSpline_test.cpp")
    target_link_libraries(testMonthlyCubicBSpline LINK_PUBLIC Boost::boost doctest::doctest)

    add_executable(
        testFileCallbackCloser
        "FileCallbackCloser_test.cpp"
        "../../core/src/FileCallbackCloser.cpp"
    )
    target_link_libraries(testFileCallbackCloser PRIVATE doctest::doctest)

    set(MODEL_INCLUDE_DIR "./testmodelarraydetails")

    add_executable(
        testModelArray
        "ModelArray_test.cpp"
        "../src/ModelArray.cpp"
        "${MODEL_INCLUDE_DIR}/ModelArrayDetails.cpp"
    )
    target_include_directories(testModelArray PRIVATE "../src" "${MODEL_INCLUDE_DIR}")
    target_link_libraries(testModelArray PRIVATE doctest::doctest Eigen3::Eigen)

    # The help config test needs access to the full binary
    file(CREATE_LINK ${NEXTSIM_BINARY_PATH} ${CMAKE_CURRENT_BINARY_DIR}/nextsim SYMBOLIC)
endif()
