# - Try to find Eigen3 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen3 3.1.2)
# to require version 3.1.2 or newer of Eigen3.
#
# Once done this will define
#
#  EIGEN_FOUND - system has eigen lib with correct version
#  Eigen_INCLUDE_DIR - the eigen include directory
#  EIGEN_VERSION - eigen version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

if(NOT Eigen_FIND_VERSION)
  if(NOT Eigen_FIND_VERSION_MAJOR)
    set(Eigen_FIND_VERSION_MAJOR 2)
  endif(NOT Eigen_FIND_VERSION_MAJOR)
  if(NOT Eigen_FIND_VERSION_MINOR)
    set(Eigen_FIND_VERSION_MINOR 91)
  endif(NOT Eigen_FIND_VERSION_MINOR)
  if(NOT Eigen_FIND_VERSION_PATCH)
    set(Eigen_FIND_VERSION_PATCH 0)
  endif(NOT Eigen_FIND_VERSION_PATCH)

  set(Eigen_FIND_VERSION "${Eigen_FIND_VERSION_MAJOR}.${Eigen_FIND_VERSION_MINOR}.${Eigen_FIND_VERSION_PATCH}")
endif(NOT Eigen_FIND_VERSION)


# Construct consistent error messages for use below.
set(EIGEN_DIR_DESCRIPTION "directory containing the file 'Eigen/Core', i.e. the root of the build tree, or the PREFIX/include/ for an installation.")
set(EIGEN_DIR_MESSAGE "EIGEN not found.  Set the Eigen_INCLUDE_DIR cmake cache entry to the ${EIGEN_DIR_DESCRIPTION}")


if(NOT EIGEN_FOUND)
  LIST(APPEND EIGEN_CHECK_INCLUDE_DIRS
    /usr/local/include
    /usr/local/homebrew/include
    /opt/local/var/macports/software 
    /opt/local/include
    /usr/include
    ${CMAKE_INSTALL_PREFIX}/include/)
  # Additional suffixes to try appending to each search path.
  LIST(APPEND EIGEN_CHECK_PATH_SUFFIXES
    eigen3 # Default root directory for Eigen.
    Eigen/include/eigen3 ) 

  # Search supplied hint directories first if supplied.
  FIND_PATH(EIGEN_INCLUDE_DIR
    NAMES Eigen/Core
    PATHS ${EIGEN_INCLUDE_DIR_HINTS}
    ${EIGEN_CHECK_INCLUDE_DIRS}
    PATH_SUFFIXES ${EIGEN_CHECK_PATH_SUFFIXES})

  if(EIGEN_INCLUDE_DIR)
    if(EXISTS ${EIGEN_INCLUDE_DIR}/Eigen/Core)
      set(EIGEN_FOUND 1)
    else()
      set(EIGEN_INCLUDE_DIR "EIGEN_INCLUDE_DIR-NOTFOUND" CACHE PATH "The ${EIGEN_DIR_DESCRIPTION}" FORCE)
    endif()
  endif()

endif()



if (NOT EIGEN_FOUND)

  # Eigen not found, explain to the user how to specify its location.
  if(EIGEN_FIND_REQUIRED)
    message(FATAL_ERROR ${EIGEN_DIR_MESSAGE})
  else()
    if(NOT EIGEN_FIND_QUIETLY)
      message(STATUS ${EIGEN_DIR_MESSAGE})
    endif()
  endif()

endif()


