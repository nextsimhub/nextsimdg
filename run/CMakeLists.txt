add_custom_target(netCDF_files ALL
    DEPENDS
    init_rect10x10.nc
)

add_custom_command(
    OUTPUT init_rect10x10.nc
    COMMAND ${Python3_EXECUTABLE} make_init_rect10x10.py
    MAIN_DEPENDENCY make_init_rect10x10.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

