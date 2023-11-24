set(COMMON_INCLUDE_DIRS
    "../src/include"
    "../../core/src"
    "../../core/src/modules"
    "../../core/src/include"
    )

set(PHYSICS_INCLUDE_DIRS
    "../../physics/src"
    "../../physics/src/modules"
    )

include_directories(${COMMON_INCLUDE_DIRS})

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

set(MODEL_INCLUDE_DIR "./testmodelarraydetails")

add_executable(testModelArray "ModelArray_test.cpp" "../src/ModelArray.cpp" "${MODEL_INCLUDE_DIR}/ModelArrayDetails.cpp")
target_include_directories(testModelArray PRIVATE "../src" "${MODEL_INCLUDE_DIR}")
target_link_libraries(testModelArray PRIVATE doctest::doctest Eigen3::Eigen)

set(MODEL_INCLUDE_DIR "../../core/src/discontinuousgalerkin")

add_executable(testRectGrid "RectGrid_test.cpp")
target_include_directories(testRectGrid PRIVATE ${MODEL_INCLUDE_DIR})
target_link_libraries(testRectGrid PRIVATE nextsimlib doctest::doctest)

add_executable(testParaGrid "ParaGrid_test.cpp")
target_include_directories(testParaGrid PRIVATE ${MODEL_INCLUDE_DIR})
target_link_libraries(testParaGrid PRIVATE nextsimlib doctest::doctest)

add_executable(testModelComponent "ModelComponent_test.cpp")
target_include_directories(testModelComponent PRIVATE ${MODEL_INCLUDE_DIR})
target_link_libraries(testModelComponent PRIVATE nextsimlib doctest::doctest)

add_executable(testModelArrayRef "ModelArrayRef_test.cpp")
target_include_directories(testModelArrayRef PRIVATE ${MODEL_INCLUDE_DIR})
target_link_libraries(testModelArrayRef PRIVATE nextsimlib doctest::doctest)

# PrognosticData (and hopefully that alone) requires code from the physics tree
add_executable(testPrognosticData "PrognosticData_test.cpp" "DynamicsModuleForPDtest.cpp")
target_include_directories(testPrognosticData PRIVATE ${PHYSICS_INCLUDE_DIRS} ${MODEL_INCLUDE_DIR})
target_link_libraries(testPrognosticData PRIVATE nextsimlib doctest::doctest)

add_executable(testConfigOutput "ConfigOutput_test.cpp")
target_include_directories(testConfigOutput PRIVATE ${MODEL_INCLUDE_DIRS})
target_link_libraries(testConfigOutput PRIVATE nextsimlib doctest::doctest)

add_executable(testMonthlyCubicBSpline
    "MonthlyCubicBSpline_test.cpp"
    )
target_include_directories(testMonthlyCubicBSpline PUBLIC ${MODEL_INCLUDE_DIRS})
target_link_libraries(testMonthlyCubicBSpline LINK_PUBLIC Boost::boost doctest::doctest)

add_executable(testFileCallbackCloser
    "FileCallbackCloser_test.cpp"
    "../../core/src/FileCallbackCloser.cpp"
    )
target_include_directories(testFileCallbackCloser PRIVATE ${MODEL_INCLUDE_DIRS})
target_link_libraries(testFileCallbackCloser PRIVATE doctest::doctest)
