cmake_minimum_required(VERSION 2.6)
project(SUCCINCT)

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  enable_testing()
endif()

find_package(Boost 1.42.0 COMPONENTS unit_test_framework iostreams system thread filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
link_directories (${Boost_LIBRARY_DIRS})

include_directories(${PROJECT_SOURCE_DIR})

set(SUCCINCT_SOURCES
  rs_bit_vector.cpp
  bp_vector.cpp
  )

add_library(succinct STATIC ${SUCCINCT_SOURCES})

file(GLOB SUCCINCT_TEST_SOURCES test_*.cpp)
foreach(TEST_SRC ${SUCCINCT_TEST_SOURCES})
  get_filename_component (TEST_SRC_NAME ${TEST_SRC} NAME_WE)
  add_executable(${TEST_SRC_NAME} ${TEST_SRC})
  target_link_libraries(${TEST_SRC_NAME}
    succinct
    ${Boost_LIBRARIES}
    )
  add_test(${TEST_SRC_NAME} ${TEST_SRC_NAME})
endforeach(TEST_SRC)

add_subdirectory(perftest)

