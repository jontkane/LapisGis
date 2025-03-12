cmake_minimum_required(VERSION 3.24)
project(LapisGis VERSION 0.9)

if(" ${CMAKE_SOURCE_DIR}" STREQUAL " ${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "
FATAL: In-source builds are not allowed.
       You should create a separate directory for build files.
")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})

file(GLOB LAPIS_GIS_SOURCES
	${CMAKE_SOURCE_DIR}/src/*.hpp
	${CMAKE_SOURCE_DIR}/src/*.cpp)

#not using add_subdirectory here because lazperf generates a very annoying number of targets
file(GLOB_RECURSE LAZPERF_FILES 
	${CMAKE_SOURCE_DIR}/src/lazperf/cpp/lazperf/*.cpp
	${CMAKE_SOURCE_DIR}/src/lazperf/cpp/lazperf/*.hpp
)
add_library(lazperf STATIC ${LAZPERF_FILES})

add_library(Lapis_gis STATIC ${LAPIS_GIS_SOURCES})

find_package(GDAL REQUIRED)
find_package(Boost COMPONENTS program_options REQUIRED)
find_package(GeoTIFF REQUIRED)
find_package(PROJ REQUIRED)
find_package(xtl REQUIRED)

set(LAPISGIS_EXTERNAL_INCLUDES
	${GDAL_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${GeoTIFF_INCLUDE_DIRS}
	${PROJ_INCLUDE_DIR}
	${xtl_INCLUDE_DIR}
	${CMAKE_SOURCE_DIR}/src/lazperf/cpp
	)

target_include_directories(Lapis_gis PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_precompile_headers(Lapis_gis PRIVATE ${CMAKE_SOURCE_DIR}/src/gis_pch.hpp)

if (MSVC)
	target_compile_options(Lapis_gis PRIVATE /W3 /WX)
else()
	target_compile_options(Lapis_gis PRIVATE -Wall -Wextra -Werror)
endif()
