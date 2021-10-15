# Build the main neXtSIM model



set(SRC_DIR "../dynamics")
set(INCLUDE_DIR "${SRC_DIR}/include")

include_directories(${INCLUDE_DIR})

add_executable(example1 example1.cpp)

target_link_libraries (example1 Eigen3::Eigen nextsimdyn)

