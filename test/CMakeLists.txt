# Generate input files
add_custom_command(
  OUTPUT init_25km_NH.nc
  COMMAND ncgen -k netCDF-4 -o init_25km_NH.nc init_25km_NH.cdl
  DEPENDS init_25km_NH.cdl
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
add_custom_command(
  OUTPUT 25km_NH.ERA5_2010-01-01T000000_test_data.nc
  COMMAND ncgen -k netCDF-4 -o 25km_NH.ERA5_2010-01-01T000000_test_data.nc
          25km_NH.ERA5_2010-01-01T000000_test_data.cdl
  DEPENDS 25km_NH.ERA5_2010-01-01T000000_test_data.cdl
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
add_custom_command(
  OUTPUT 25km_NH.TOPAZ4_2010-01-01T000000_test_data.nc
  COMMAND ncgen -k netCDF-4 -o 25km_NH.TOPAZ4_2010-01-01T000000_test_data.nc
          25km_NH.TOPAZ4_2010-01-01T000000_test_data.cdl
  DEPENDS 25km_NH.TOPAZ4_2010-01-01T000000_test_data.cdl
  WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
add_custom_target(generate_test_input_files ALL DEPENDS init_25km_NH.nc)
