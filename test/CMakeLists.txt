# Generate the example feature and target data files for testing
set(BASE_NAMES
    "init_25km_NH"
    "25km_NH.ERA5_2010-01-01T000000_test_data"
    "25km_NH.TOPAZ4_2010-01-01T000000_test_data"
)
foreach(BASE_NAME IN LISTS BASE_NAMES)
    set(OUTPUT_FILE ${BASE_NAME}.nc)
    set(INPUT_FILE ${BASE_NAME}.cdl)
    add_custom_command(
        OUTPUT ${OUTPUT_FILE}
        COMMAND ncgen -k netCDF-4 -o ${OUTPUT_FILE} ${INPUT_FILE}
        DEPENDS ${INPUT_FILE}
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    )
    list(APPEND GENERATED_FILES ${OUTPUT_FILE})
endforeach()
add_custom_target(generate_test_input_files ALL DEPENDS ${GENERATED_FILES})
