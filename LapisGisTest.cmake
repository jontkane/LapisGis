file(GLOB LAPISGIS_TEST_SOURCES
	${LAPISGIS_DIR}/src/test/*.hpp
	${LAPISGIS_DIR}/src/test/*.cpp)

add_compile_definitions(LAPISGISTESTFILES="${LAPISGIS_DIR}/src/test/testfiles/")
add_executable(LapisGis_test ${LAPISGIS_TEST_SOURCES})
find_package(GTest REQUIRED)
target_include_directories(LapisGis_test PRIVATE ${LAPISGIS_EXTERNAL_INCLUDES})
target_link_libraries(LapisGis_test PRIVATE ${LAPISGIS_EXTERNAL_LINKS})
target_link_libraries(LapisGis_test PRIVATE LapisGis)
target_include_directories(LapisGis_test PRIVATE ${GTEST_INCLUDE_DIRS})
target_link_libraries(LapisGis_test PRIVATE ${GTEST_BOTH_LIBRARIES})


if (MSVC)
	target_compile_options(LapisGis_test PRIVATE /W3 /WX)
else()
	target_compile_options(LapisGis_test PRIVATE -Wall -Wextra -Werror)
endif()