# Build the main neXtSIM model

set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(BaseSources
	"main.cpp"
	"Logged.cpp"
	"Timed.cpp"
	"CountTimer.cpp" # TODO: basic timer for testing only
	"Model.cpp"
	"Iterator.cpp"
	"SimpleIterant.cpp"
	)

set(BaseHeaders
	"${INC}/Logged.hpp"
	"${INC}/Timed.hpp"
	"${INC}/Timer.hpp"
	"${INC}/Model.hpp"
	"${INC}/Iterator.hpp"
	"${INC}/SimpleIterant.hpp"
	)
	
set(NextsimSources
	${BaseHeaders}
	${BaseSources}
	)

add_executable(nextsim ${NextsimSources})
target_include_directories(nextsim PRIVATE "include" ${Boost_INCLUDE_DIRS})
target_link_libraries(nextsim LINK_PUBLIC ${Boost_Libraries})