# Generate the example restart files that don't depend on external data
add_custom_command(
    OUTPUT init_rect20x30.nc
    COMMAND ncgen -k netCDF-4 -o init_rect20x30.nc init_rect20x30.cdl
    DEPENDS init_rect20x30.cdl
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_custom_command(
    OUTPUT init_column.nc
    COMMAND ncgen -k netCDF-4 -o init_column.nc init_column.cdl
    DEPENDS init_column.cdl
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init_para24x30.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
execute_process(
    COMMAND "${Python_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/make_init_benchmark.py"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
)
add_custom_target(generate_run_input_files ALL DEPENDS init_rect20x30.nc init_column.nc)
