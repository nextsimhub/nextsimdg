# Build the main neXtSIM model

set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(BaseSources
	"dgtransport.cpp"
	)

set(BaseHeaders
	"${INC}/dgbasisfunctions_gausspoints.hpp"
	"${INC}/dgtransport.hpp"
	"${INC}/mesh.hpp"
	"${INC}/dginitial.hpp"
	"${INC}/dgvector.hpp"
	"${INC}/timemesh.hpp"
  "${INC}/dglimiters.hpp"
	"${INC}/dgtimestepping.hpp"
	"${INC}/dgvector_manipulations.hpp"
	"${INC}/dgvisu.hpp"
	)


set(NextsimDynSources
	${BaseHeaders}
	${BaseSources}
	)

add_library(nextsimdyn ${NextsimDynSources})

target_include_directories(nextsimdyn PRIVATE "include")

target_link_libraries(nextsimdyn Eigen3::Eigen)

