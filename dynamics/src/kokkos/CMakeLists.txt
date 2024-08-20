FILE(GLOB KokkosSources *.cpp)

set(NextsimSources
	"${NextsimSources}"
    "${KokkosSources}"
    PARENT_SCOPE)