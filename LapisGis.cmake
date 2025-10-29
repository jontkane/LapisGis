include_guard(GLOBAL)

set(PROJ_DB_PATH "" CACHE FILEPATH "Path to proj.db")
if(NOT EXISTS "${PROJ_DB_PATH}")
    message(WARNING "proj.db not found at: ${PROJ_DB_PATH}\n"
	"If you need LapisGis to handle proj.db, please specify it using -DPROJ_DB_PATH")
endif()

function(copy_proj_db_after_build target)
    add_custom_command(TARGET ${target} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${PROJ_DB_PATH}"
            "$<TARGET_FILE_DIR:${target}>/proj.db"
        COMMENT "Copying proj.db to output directory for ${target}"
    )
	target_compile_definitions(LapisGis PRIVATE LAPISGIS_PROJDB_IN_EXE_DIR)
endfunction()

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

find_package(GDAL 3.11 REQUIRED)
find_package(GeoTIFF 1.7.4 REQUIRED)
find_package(PROJ 9.6.1 REQUIRED)
find_package(xtl 0.8.0 REQUIRED)

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
