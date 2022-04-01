### Set the list of all cpp files
file(GLOB NextsimDGLib_Sources *.cpp)
add_library(${NextsimDGLib} ${NextsimDGLib_Sources})
# Create an alias for our library
add_library(${PROJECT_NAME}::${NextsimDGLib} ALIAS ${NextsimDGLib})
target_include_directories(${NextsimDGLib}
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
target_link_libraries(${NextsimDGLib} PUBLIC Eigen3::Eigen)

# Set the header files that need to be installed
file(GLOB NextsimDGLib_Headers *.hpp)
set_target_properties(${NextsimDGLib} PROPERTIES PUBLIC_HEADER "${NextsimDGLib_Headers}")

install(TARGETS ${NextsimDGLib}
  EXPORT ${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
