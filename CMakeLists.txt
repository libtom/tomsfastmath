# SPDX-License-Identifier: Unlicense
#
# TomsFastMath, a fast ISO C bignum library.

cmake_minimum_required(VERSION 3.10)

set(LIBRARY_NAME "libtfm")
project(tomsfastmath
    VERSION 0.13.1
    DESCRIPTION "A fast public domain, open source, large integer arithmetic library written in portable ISO C."
    HOMEPAGE_URL "https://www.libtom.net/TomsFastMath"
    LANGUAGES C)

# package release version
# bump if re-releasing the same VERSION + patches
# set to 1 if releasing a new VERSION
set(PACKAGE_RELEASE_VERSION 1)

#-----------------------------------------------------------------------------
# Include cmake modules
#-----------------------------------------------------------------------------
include(GNUInstallDirs)
include(CheckIPOSupported)
include(CMakePackageConfigHelpers)
# default is "No tests"
option(BUILD_TESTING "" OFF)
include(CTest)
include(sources.cmake)

# The only direct cmake argument for now
option(BUILD_SHARED_LIBS "Build shared library and only the shared library if \"ON\", default is static" OFF)

#-----------------------------------------------------------------------------
# Compose CFLAGS
#-----------------------------------------------------------------------------

# Some information ported from makefile_include.mk


if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(TFM_C_FLAGS -Wall -Wsign-compare -Wextra -Wshadow
                -Wdeclaration-after-statement -Wbad-function-cast -Wcast-align
                -Wstrict-prototypes -Wpointer-arith -Wsystem-headers)
set(CMAKE_C_FLAGS_DEBUG "-g3")
set(CMAKE_C_FLAGS_RELEASE "-O3 -fomit-frame-pointer")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "-g3 -O2")
set(CMAKE_C_FLAGS_MINSIZEREL "-Os")

# What compiler do we have and what are their...uhm... peculiarities
if(CMAKE_C_COMPILER_ID MATCHES "(C|c?)lang")
    list(APPEND TFM_C_FLAGS -Wno-typedef-redefinition -Wno-tautological-compare -Wno-builtin-requires-header)
endif()
if(CMAKE_C_COMPILER MATCHES "mingw")
    list(APPEND TFM_C_FLAGS -Wno-shadow -Wno-expansion-to-defined -Wno-declaration-after-statement -Wno-bad-function-cast)
endif()
if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    list(APPEND TFM_C_FLAGS -Wno-nullability-completeness)
else()
    list(APPEND TFM_C_FLAGS  -funroll-loops)
endif()
if(CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
    list(APPEND TFM_C_FLAGS -no-undefined)
endif()

# TODO: coverage (lgcov)

# If the user set the environment variables at generate-time, append them
# in order to allow overriding our defaults.
# ${TFM_CFLAGS} means the user passed it via sth like:
# $ cmake -DTFM_CFLAGS="foo"
list(APPEND TFM_C_FLAGS ${TFM_CFLAGS})
list(APPEND TFM_LD_FLAGS ${TFM_LDFLAGS})

#-----------------------------------------------------------------------------
# library target
#-----------------------------------------------------------------------------
add_library(${PROJECT_NAME}
    ${SOURCES}
    ${PUBLIC_HEADERS}
    ${PRIVATE_HEADERS}
)

target_include_directories(${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/headers>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}>
)

target_compile_options(${PROJECT_NAME} BEFORE PRIVATE
    ${TFM_C_FLAGS}
)
target_link_options(${PROJECT_NAME} BEFORE PRIVATE
    ${TFM_LD_FLAGS}
)

set_target_properties(${PROJECT_NAME} PROPERTIES
    OUTPUT_NAME tfm
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
    PUBLIC_HEADER "${PUBLIC_HEADERS}"
)

option(COMPILE_LTO "Build with LTO enabled")
if(COMPILE_LTO)
    check_ipo_supported(RESULT COMPILER_SUPPORTS_LTO)
    if(COMPILER_SUPPORTS_LTO)
        set_property(TARGET ${PROJECT_NAME} PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    else()
        message(SEND_ERROR "This compiler does not support LTO. Reconfigure ${PROJECT_NAME} with -DCOMPILE_LTO=OFF.")
    endif()
endif()

#-----------------------------------------------------------------------------
# demo target
#-----------------------------------------------------------------------------

if(BUILD_TESTING)
    enable_testing()
    add_subdirectory(demo)
endif()

#-----------------------------------------------------------------------------
# Install/export targets and files
#-----------------------------------------------------------------------------
set(CONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
set(PROJECT_VERSION_FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake")
set(PROJECT_CONFIG_FILE "${PROJECT_NAME}-config.cmake")
set(TARGETS_EXPORT_NAME "${PROJECT_NAME}Targets")

install(TARGETS ${PROJECT_NAME}
    EXPORT ${TARGETS_EXPORT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT Libraries
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)

# Install tomsfastmath.pc for pkg-config if we build a shared library
if(BUILD_SHARED_LIBS)
    # Let the user override the default directory of the pkg-config file (usually this shouldn't be required to be changed)
    set(CMAKE_INSTALL_PKGCONFIGDIR "${CMAKE_INSTALL_LIBDIR}/pkgconfig" CACHE PATH "Folder where to install .pc files")

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}.pc.in
        ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
        @ONLY
    )

    install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}.pc
        DESTINATION ${CMAKE_INSTALL_PKGCONFIGDIR}
    )
endif()

# generate package version file
write_basic_package_version_file(
    ${PROJECT_VERSION_FILE}
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

# install version file
install(FILES ${PROJECT_VERSION_FILE}
        DESTINATION ${CONFIG_INSTALL_DIR}
)

# build directory package config
export(EXPORT ${TARGETS_EXPORT_NAME}
       FILE ${PROJECT_CONFIG_FILE}
)

# installed package config
install(EXPORT ${TARGETS_EXPORT_NAME}
        DESTINATION ${CONFIG_INSTALL_DIR}
        FILE ${PROJECT_CONFIG_FILE}
)

# add to CMake registry
export(PACKAGE ${PROJECT_NAME})

#---------------------------------------------------------------------------------------
# Create release packages
#---------------------------------------------------------------------------------------

# determine distribution and architecture
find_program(LSB_RELEASE lsb_release)
find_program(SYSCTL sysctl)
find_program(UNAME uname)

if(UNAME)
    execute_process(COMMAND uname -m OUTPUT_VARIABLE MACHINE_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
elseif(SYSCTL)
    execute_process(COMMAND sysctl -b hw.machine_arch OUTPUT_VARIABLE MACHINE_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
else()
    string(TOLOWER ${CMAKE_SYSTEM_NAME} MACHINE_ARCH)
endif()

if(LSB_RELEASE)
    execute_process(COMMAND lsb_release -si OUTPUT_VARIABLE LINUX_DISTRO OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND lsb_release -sc OUTPUT_VARIABLE LINUX_DISTRO_CODENAME OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND lsb_release -sr OUTPUT_VARIABLE LINUX_DISTRO_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

    string(TOLOWER ${LINUX_DISTRO} LINUX_DISTRO)
    if(LINUX_DISTRO_CODENAME STREQUAL "n/a")
        set(DISTRO_PACK_PATH ${LINUX_DISTRO}/${LINUX_DISTRO_VERSION}/)
    else()
        set(DISTRO_PACK_PATH ${LINUX_DISTRO}/${LINUX_DISTRO_CODENAME}/)
    endif()
else()
    set(DISTRO_PACK_PATH ${CMAKE_SYSTEM_NAME}/)
endif()

# make sure untagged versions get a different package name
execute_process(COMMAND git describe --exact-match --tags ERROR_QUIET RESULT_VARIABLE REPO_HAS_TAG)
if(REPO_HAS_TAG EQUAL 0)
    set(PACKAGE_NAME_SUFFIX "")
else()
    set(PACKAGE_NAME_SUFFIX "-git")
    message(STATUS "Use -git suffix")
endif()

# default CPack generators
set(CPACK_GENERATOR TGZ STGZ)

# extra CPack generators
if(LINUX_DISTRO STREQUAL "debian" OR LINUX_DISTRO STREQUAL "ubuntu" OR LINUX_DISTRO STREQUAL "linuxmint")
    list(APPEND CPACK_GENERATOR DEB)
elseif(LINUX_DISTRO STREQUAL "fedora" OR LINUX_DISTRO STREQUAL "opensuse" OR LINUX_DISTRO STREQUAL "centos")
    list(APPEND CPACK_GENERATOR RPM)
elseif(CMAKE_SYSTEM_NAME STREQUAL "FreeBSD")
    list(APPEND CPACK_GENERATOR FREEBSD)
endif()

set(TFM_DEBIAN_SHARED_PACKAGE_NAME "${LIBRARY_NAME}${PACKAGE_NAME_SUFFIX}${PROJECT_VERSION_MAJOR}")

# general CPack config
set(CPACK_PACKAGE_DIRECTORY ${CMAKE_BINARY_DIR}/packages/${DISTRO_PACK_PATH})
message(STATUS "CPack: packages will be generated under ${CPACK_PACKAGE_DIRECTORY}")
if(BUILD_SHARED_LIBS)
    set(CPACK_PACKAGE_NAME "${LIBRARY_NAME}${PROJECT_VERSION_MAJOR}")
    set(CPACK_DEBIAN_PACKAGE_NAME "${TFM_DEBIAN_SHARED_PACKAGE_NAME}")
else()
    set(CPACK_PACKAGE_NAME "${LIBRARY_NAME}-devel")
    set(CPACK_DEBIAN_PACKAGE_NAME "${LIBRARY_NAME}-dev")
    set(CPACK_DEBIAN_LIBRARIES_PACKAGE_NAME "${LIBRARY_NAME}${PACKAGE_NAME_SUFFIX}-dev")
endif()
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "TomsFastMath")
set(CPACK_PACKAGE_VENDOR "libtom projects")
set(CPACK_PACKAGE_CONTACT "libtom@googlegroups.com")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")
set(PACKAGE_NAME_TRAILER ${CPACK_PACKAGE_VERSION}-${PACKAGE_RELEASE_VERSION}_${MACHINE_ARCH})
set(CPACK_PACKAGE_FILE_NAME ${CPACK_PACKAGE_NAME}-${PACKAGE_NAME_TRAILER})

# deb specific CPack config
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_DEBIAN_DEBUGINFO_PACKAGE ON)
set(CPACK_DEBIAN_PACKAGE_RELEASE ${PACKAGE_RELEASE_VERSION})
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
if(BUILD_SHARED_LIBS)
    set(CPACK_DEBIAN_PACKAGE_SECTION "libs")
else()
    set(CPACK_DEBIAN_PACKAGE_SECTION "devel")
    set(CPACK_DEBIAN_PACKAGE_DEPENDS ${TFM_DEBIAN_SHARED_PACKAGE_NAME})
    set(CPACK_DEB_COMPONENT_INSTALL ON)
    set(CPACK_ARCHIVE_COMPONENT_INSTALL ON)
    set(CPACK_COMPONENTS_ALL Libraries)
endif()
## in case we build our `-git` version of the package, mark the packages as "provides:" and "conflicts:"
## so only this package will be installed instead of the regular one.
if(NOT REPO_HAS_TAG EQUAL 0)
    set(CPACK_DEBIAN_PACKAGE_CONFLICTS ${CPACK_PACKAGE_NAME})
    set(CPACK_DEBIAN_PACKAGE_PROVIDES ${CPACK_PACKAGE_NAME})
    set(CPACK_DEBIAN_LIBRARIES_PACKAGE_CONFLICTS ${CPACK_DEBIAN_PACKAGE_NAME})
    set(CPACK_DEBIAN_LIBRARIES_PACKAGE_PROVIDES ${CPACK_DEBIAN_PACKAGE_NAME})
endif()

# rpm specific CPack config
set(CPACK_RPM_PACKAGE_RELEASE ${PACKAGE_RELEASE_VERSION})
set(CPACK_RPM_PACKAGE_ARCHITECTURE ${MACHINE_ARCH})
set(CPACK_RPM_PACKAGE_NAME "${CPACK_PACKAGE_NAME}-${PROJECT_VERSION}")
set(CPACK_RPM_PACKAGE_LICENSE "The Unlicense")

# FreeBSD specific CPack config
set(CPACK_FREEBSD_PACKAGE_MAINTAINER "ports@FreeBSD.org")
set(CPACK_FREEBSD_PACKAGE_ORIGIN "math/tomsfastmath")
set(CPACK_FREEBSD_PACKAGE_CATEGORIES "math")

include(CPack)
