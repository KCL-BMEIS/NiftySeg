PROJECT(NiftySeg)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)
include(ExternalProject)


#-----------------------------------------------------------------------------
# Preflight Checks - CMake Version, InSource Build Type
#-----------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8.0)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.8\\.0$")
 mark_as_advanced(FORCE CMAKE_BACKWARDS_COMPATIBILITY)
else("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.8\\.0$")
 mark_as_advanced(CLEAR CMAKE_BACKWARDS_COMPATIBILITY)
endif("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.8\\.0$")

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
        message("In-source builds not allowed by niftyseg.")
        message("Please make a new directory (i.e. a build directory) and run CMake from there.")
        message(FATAL_ERROR "You may need to remove CMakeCache.txt.")
endif(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})


#-----------------------------------------------------------------------------
# Overall Options
#-----------------------------------------------------------------------------

OPTION(BUILD_ALL "To build the LoAd segmentation tool" ON)
OPTION(BUILD_SHARED "Build the libraries as shared." OFF)
OPTION(INSTALL_PRIORS "Install prior images. Turn this off if you are really short on disk space" ON)


#-----------------------------------------------------------------------------
# Set build details - C flags, Build Type, Library_type, OpenMP
#-----------------------------------------------------------------------------

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DFILE_OFFSETBITS=64 -DDIRENT=dirent64 -DREADDIR=readdir64")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DFILE_OFFSETBITS=64 -DDIRENT=dirent64 -DREADDIR=readdir64")

IF (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   add_definitions(/W1)
ELSE (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
   add_definitions(-Wall -pedantic)
ENDIF()

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(NOT cmake_build_type_tolower STREQUAL "debug"
         AND NOT cmake_build_type_tolower STREQUAL "release"
         AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
        message("Unknown build type \"${CMAKE_BUILD_TYPE}\".")
        message(FATAL_ERROR "Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

if(USE_OPENMP)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(USE_OPENMP)

IF(BUILD_SHARED)
  SET(LIBRARY_TYPE SHARED)
ELSE(BUILD_SHARED)
  SET(LIBRARY_TYPE STATIC)
ENDIF(BUILD_SHARED)

if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")
  # OpenMP is disable when the compiler is based on Clang
  find_package(OpenMP)
  if(OPENMP_FOUND)
    option(USE_OPENMP "To use openMP for multi-CPU processing" ON)
  endif(OPENMP_FOUND)
endif(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")


#-----------------------------------------------------------------------------
# Include Dirs and Module Paths
#-----------------------------------------------------------------------------

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/seg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nifti)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/zlib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/win32)

SET(CMAKE_MODULE_PATH 
${CMAKE_BINARY_DIR}
${CMAKE_SOURCE_DIR}/CMake
${CMAKE_MODULE_PATH}
)


#-----------------------------------------------------------------------------
# Optional installation of NiftyReg packages
#-----------------------------------------------------------------------------

OPTION(INSTALL_NIFTYREG "Install the NiftyReg package" OFF)
IF(INSTALL_NIFTYREG)
   OPTION(USE_SSE "Use SSE in NiftyReg package" OFF)
   ExternalProject_Add(NiftyReg
               SOURCE_DIR ${PROJECT_BINARY_DIR}/niftyreg/source
               PREFIX ${PROJECT_BINARY_DIR}/niftyreg
               STAMP_DIR ${PROJECT_BINARY_DIR}/niftyreg/stamps
               TMP_DIR ${PROJECT_BINARY_DIR}/niftyreg/tmp
               DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/niftyreg
               BINARY_DIR ${PROJECT_BINARY_DIR}/niftyreg/build
               GIT_REPOSITORY git://git.code.sf.net/p/niftyreg/git
               TIMEOUT 100
               UPDATE_COMMAND ""
               CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DUSE_SSE=${USE_SSE} -DUSE_OPENMP=${USE_OPENMP}
               BUILD_COMMAND ${CMAKE_BUILD_TOOL}
               INSTALL_COMMAND ${CMAKE_BUILD_TOOL} install
               )
ENDIF(INSTALL_NIFTYREG)


#-----------------------------------------------------------------------------
# Search for Eigen. Download and install if not found
#-----------------------------------------------------------------------------

find_package(EIGEN)
IF(EIGEN_FOUND)
    # if EIGEN is found, then just set its include dir
    set(EIGEN_INCLUDE_DIR  ${CMAKE_INSTALL_PREFIX}/include)
    MARK_AS_ADVANCED(EIGEN_INCLUDE_DIR)
ELSE(EIGEN_FOUND)
    # Else, download, install, and set the include dir
    set(eigen_url http://cmictig.cs.ucl.ac.uk/images/eigen-eigen-6b38706d90a9.tar.bz2)
    set(eigen_md5 ece1dbf64a49753218ce951624f4c487)
    ExternalProject_Add(Eigen
        SOURCE_DIR ${PROJECT_BINARY_DIR}/Eigen
        URL ${eigen_url}
        URL_MD5 ${eigen_md5}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ""
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy_directory "${PROJECT_BINARY_DIR}/Eigen/Eigen" "${CMAKE_INSTALL_PREFIX}/include/Eigen"
    )
    set(EIGEN_INCLUDE_DIR  ${CMAKE_INSTALL_PREFIX}/include)
    MARK_AS_ADVANCED(EIGEN_INCLUDE_DIR)
ENDIF(EIGEN_FOUND)


#-----------------------------------------------------------------------------
# Actually install NiftySeg
#-----------------------------------------------------------------------------

include_directories(${EIGEN_INCLUDE_DIR})
SET(ZLIB "z")
SUBDIRS(zlib)
SUBDIRS(nifti)
SUBDIRS(seg-lib)
SUBDIRS(seg-apps)

# Install scripts
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/LoAd_brainonly.sh ${CMAKE_BINARY_DIR}/LoAd_brainonly.sh @ONLY )
INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/LoAd_brainonly.sh DESTINATION ${CMAKE_INSTALL_PREFIX}/bin COMPONENT applications)

# Install priors
SET(INSTALL_PRIORS_DIRECTORY ${CMAKE_INSTALL_PREFIX}/priors)
MARK_AS_ADVANCED(INSTALL_PRIORS_DIRECTORY)
IF(INSTALL_PRIORS)
   FILE(GLOB PRIOR_FILES priors/*.nii.gz)
   INSTALL(PROGRAMS ${PRIOR_FILES} DESTINATION ${INSTALL_PRIORS_DIRECTORY} COMPONENT applications)
ENDIF(INSTALL_PRIORS)

#-----------------------------------------------------------------------------
