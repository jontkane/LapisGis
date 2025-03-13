set(LAPIS_GIS_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB LAPIS_GIS_SOURCES
	${LAPIS_GIS_DIR}/src/*.hpp
	${LAPIS_GIS_DIR}/src/*.cpp)

file(GLOB LAPIS_GIS_TEST_SOURCES
	${LAPIS_GIS_DIR}/src/test/*.hpp
	${LAPIS_GIS_DIR}/src/test/*.cpp)

#not using add_subdirectory here because lazperf generates a very annoying number of targets
file(GLOB_RECURSE LAZPERF_FILES 
	${LAPIS_GIS_DIR}/src/lazperf/cpp/lazperf/*.cpp
	${LAPIS_GIS_DIR}/src/lazperf/cpp/lazperf/*.hpp
)


file(GLOB WHEREAMI_SOURCES
	${LAPIS_GIS_DIR}/src/whereami/*.h
	${LAPIS_GIS_DIR}/src/whereami/*.c)

add_library(Lapis_gis STATIC ${LAPIS_GIS_SOURCES})
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
	${LAPIS_GIS_DIR}/src/lazperf/cpp
	${LAPIS_GIS_DIR}/src/whereami
	)

set(LAPISGIS_EXTERNAL_LINKS
	${GDAL_LIBRARIES}
	${GeoTIFF_LIBRARIES}
	${PROJ_LIBRARIES}
	lazperf
	whereami
	)

target_include_directories(Lapis_gis PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_precompile_headers(Lapis_gis PRIVATE ${LAPIS_GIS_DIR}/src/gis_pch.hpp)

add_compile_definitions(LAPISTESTFILES="${LAPIS_GIS_DIR}/src/test/testfiles/")
add_executable(Lapis_gis_test ${LAPIS_GIS_TEST_SOURCES})
find_package(GTest REQUIRED)
target_include_directories(Lapis_gis_test PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_link_libraries(Lapis_gis_test PRIVATE ${LAPISGIS_EXTERNAL_LINKS})
target_link_libraries(Lapis_gis_test PRIVATE Lapis_gis)
target_include_directories(Lapis_gis_test PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(Lapis_gis_test PRIVATE ${GTEST_BOTH_LIBRARIES})

if (MSVC)
	target_compile_options(Lapis_gis PRIVATE /W3 /WX)
	target_compile_options(Lapis_gis_test PRIVATE /W3 /WX)
	target_compile_options(lazperf PRIVATE /W0)
else()
	target_compile_options(Lapis_gis PRIVATE -Wall -Wextra -Werror)
	target_compile_options(Lapis_gis_test PRIVATE -Wall -Wextra -Werror)
endif()
