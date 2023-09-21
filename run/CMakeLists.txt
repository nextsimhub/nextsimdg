add_custom_target(netCDF_files ALL
    DEPENDS
    init_rect10x10.nc
    init_column.nc
    init_rect20x30.nc
)

add_custom_command(
    OUTPUT init_rect10x10.nc
    COMMAND ${Python3_EXECUTABLE} make_init_rect10x10.py
    MAIN_DEPENDENCY make_init_rect10x10.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_command(
    OUTPUT init_column.nc
    COMMAND ${Python3_EXECUTABLE} make_init_column.py
    MAIN_DEPENDENCY make_init_column.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_custom_command(
    OUTPUT init_rect20x30.nc
    COMMAND ${Python3_EXECUTABLE} make_init_rect20x30.py
    MAIN_DEPENDENCY make_init_rect20x30.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
