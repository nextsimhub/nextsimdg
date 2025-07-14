# Generate the example restart files that don't depend on external data
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init_rect20x30.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_custom_command(
    OUTPUT init_column.nc
    COMMAND ncgen -k netCDF-4 -o init_column.nc
            ${CMAKE_CURRENT_SOURCE_DIR}/init_column.cdl
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/init_column.cdl
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init_para24x30.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init_benchmark.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_custom_target(generate_input_files ALL DEPENDS init_column.nc)
