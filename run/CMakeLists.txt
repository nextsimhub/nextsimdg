# Generate the example restart files that don't depend on external data
# TODO: Convert these to the CDL approach
execute_process(
  COMMAND "${Python_EXECUTABLE}"
          "${CMAKE_CURRENT_SOURCE_DIR}/make_init_benchmark.py"
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")

# Generate NetCDF files from CDL files
set(BASE_NAMES "init_column" "init_rect20x30" "init_rect30x30"
               "init_para24x30" "init_para128x128")
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
add_custom_target(generate_run_input_files ALL DEPENDS ${GENERATED_FILES})

# Copy init_column.nc to test directory
add_custom_command(
    TARGET generate_run_input_files POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/init_column.nc"
            "${CMAKE_CURRENT_SOURCE_DIR}/../test/init_column.nc"
)
