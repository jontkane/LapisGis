include_guard(GLOBAL)

add_library(LapisGis STATIC ${LAPISGIS_SOURCES})
#add_library(lazperf STATIC ${LAZPERF_FILES})
#add_library(whereami STATIC ${WHEREAMI_SOURCES})

target_include_directories(LapisGis PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_precompile_headers(LapisGis PRIVATE ${LAPISGIS_DIR}/src/gis_pch.hpp)

if (MSVC)
	target_compile_options(LapisGis PRIVATE /W3 /WX)
	#target_compile_options(lazperf PRIVATE /W0)
else()
	target_compile_options(LapisGis PRIVATE -Wall -Wextra -Werror)
endif()
