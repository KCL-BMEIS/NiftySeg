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
set(EIGEN_DIR_DESCRIPTION "directory containing the file 'Eigen/Core', i.e. the root of the build tree, or the PREFIX/include/eigen3 for an installation.")
set(EIGEN_DIR_MESSAGE "EIGEN not found.  Set the Eigen_INCLUDE_DIR cmake cache entry to the ${EIGEN_DIR_DESCRIPTION}")

if(NOT EIGEN_FOUND)

  # Look for signature_of_eigen3_matrix_library in build trees or under <prefix>/include/eigen3.
  find_path(Eigen_INCLUDE_DIR
    NAMES Eigen/Core
    PATH_SUFFIXES eigen3
    HINTS ENV EIGEN_DIR
    PATHS

    # Help the user find it if we cannot.
    DOC "The ${EIGEN_DIR_DESCRIPTION}"
    )

  if(Eigen_INCLUDE_DIR)
    if(EXISTS ${Eigen_INCLUDE_DIR}/Eigen/Core)
      set(EIGEN_FOUND 1)
    else()
      set(Eigen_INCLUDE_DIR "Eigen_INCLUDE_DIR-NOTFOUND" CACHE PATH "The ${EIGEN_DIR_DESCRIPTION}" FORCE)
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


