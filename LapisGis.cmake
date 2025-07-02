include_guard(GLOBAL)

set(LAPISGIS_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB LAPISGIS_SOURCES
	${LAPISGIS_DIR}/src/*.hpp
	${LAPISGIS_DIR}/src/*.cpp)

#not using add_subdirectory here because lazperf generates a very annoying number of targets
file(GLOB_RECURSE LAZPERF_FILES 
	${LAPISGIS_DIR}/src/lazperf/cpp/lazperf/*.cpp
	${LAPISGIS_DIR}/src/lazperf/cpp/lazperf/*.hpp
)

file(GLOB WHEREAMI_SOURCES
	${LAPISGIS_DIR}/src/whereami/*.h
	${LAPISGIS_DIR}/src/whereami/*.c)

add_library(LapisGis STATIC ${LAPISGIS_SOURCES})
add_library(lazperf STATIC ${LAZPERF_FILES})
add_library(whereami STATIC ${WHEREAMI_SOURCES})

find_package(GDAL REQUIRED)
find_package(GeoTIFF REQUIRED)
find_package(PROJ REQUIRED)
find_package(xtl REQUIRED)

set(LAPISGIS_EXTERNAL_INCLUDES
	${GDAL_INCLUDE_DIRS}
	${GeoTIFF_INCLUDE_DIRS}
	${PROJ_INCLUDE_DIR}
	${xtl_INCLUDE_DIR}
	${LAPISGIS_DIR}/src/lazperf/cpp
	${LAPISGIS_DIR}/src/whereami
	)

set(LAPISGIS_EXTERNAL_LINKS
	${GDAL_LIBRARIES}
	${GeoTIFF_LIBRARIES}
	${PROJ_LIBRARIES}
	lazperf
	whereami
	)

target_include_directories(LapisGis PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_precompile_headers(LapisGis PRIVATE ${LAPISGIS_DIR}/src/gis_pch.hpp)

set(LAPISGIS_INCLUDES
	${LAPISGIS_EXTERNAL_INCLUDES}
	${LAPISGIS_DIR}/src
	)
set(LAPISGIS_LINKS
	${LAPISGIS_EXTERNAL_LINKS}
	LapisGis
	)

if (MSVC)
	target_compile_options(LapisGis PRIVATE /W3 /WX)
	target_compile_options(lazperf PRIVATE /W0)
else()
	target_compile_options(LapisGis PRIVATE -Wall -Wextra -Werror)
endif()
