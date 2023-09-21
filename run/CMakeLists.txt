add_custom_target(netCDF_files ALL
    DEPENDS
    init_rect10x10.nc
    init_para128x128.nc
    init_column.nc
)

add_custom_command(
    OUTPUT init_rect10x10.nc
    COMMAND ${Python3_EXECUTABLE} make_init_rect10x10.py
    MAIN_DEPENDENCY make_init_rect10x10.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

add_custom_command(
    OUTPUT init_para128x128.nc
    COMMAND ${Python3_EXECUTABLE} init_para128x128.py
    MAIN_DEPENDENCY init_para128x128.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )

    add_custom_command(
    OUTPUT init_column.nc
    COMMAND ${Python3_EXECUTABLE} make_init_column.py
    MAIN_DEPENDENCY make_init_column.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
