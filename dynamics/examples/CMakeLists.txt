# Build the main neXtSIM model

INCLUDE_DIRECTORIES(${Nextsim_INCLUDE_DIRS})
LINK_LIBRARIES(Eigen3::Eigen ${NextsimDGLib})

add_executable(multithreading_square multithreading_square.cpp)
add_executable(multithreading_rectangle1 multithreading_rectangle1.cpp)
add_executable(multithreading_rectangle2 multithreading_rectangle2.cpp)

add_executable(example1 example1.cpp)

add_executable(example2a example2a.cpp)
add_executable(example2b example2b.cpp)
add_executable(example2c example2c.cpp)



