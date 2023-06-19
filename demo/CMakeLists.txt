# SPDX-License-Identifier: Unlicense
#
# LibTomMath, a free open source portable number theoretic multiple-precision
# integer (MPI) library written entirely in C.
#

cmake_minimum_required(VERSION 3.10)

set(TFM_TEST test-tfm)

# This file can be included from the top level or used stand-alone
if(PROJECT_NAME)
    set(LIBRARY_NAME ${PROJECT_NAME})
else()
    # Define an independent project and all the necessary stuff around
    project(${TFM_TEST}
        LANGUAGES C)
    set(LIBRARY_NAME tomsfastmath)
    find_package(${LIBRARY_NAME})
    include(CTest)
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        set(CMAKE_BUILD_TYPE "Release")
    endif()
endif()

#-----------------------------------------------------------------------------
# Compose CFLAGS etc.
#-----------------------------------------------------------------------------

set(CMAKE_C_FLAGS "-DTFM_DEMO_TEST_VS_MTEST=0")
set(CMAKE_C_FLAGS_DEBUG "-g3 -O1")

#-----------------------------------------------------------------------------
# demo target
#-----------------------------------------------------------------------------

add_executable(${TFM_TEST}
    ${CMAKE_CURRENT_SOURCE_DIR}/test.c
)

target_include_directories(${TFM_TEST} PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
    $<$<STREQUAL:${PROJECT_NAME},${TFM_TEST}>:${CMAKE_CURRENT_SOURCE_DIR}/..>
)

target_link_libraries(${TFM_TEST} PRIVATE
    ${LIBRARY_NAME}
)

target_compile_options(${TFM_TEST} PRIVATE
    ${TFM_C_FLAGS}
)
target_link_options(${TFM_TEST} BEFORE PUBLIC
    ${TFM_LD_FLAGS}
)

#-----------------------------------------------------------------------------
# CTest
#-----------------------------------------------------------------------------
add_test(NAME ${TFM_TEST} COMMAND ${TFM_TEST})

find_program(MEMORYCHECK_COMMAND valgrind)
set(MEMORYCHECK_COMMAND_OPTIONS "--trace-children=yes --leak-check=full")
