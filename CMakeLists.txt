PROJECT(nifty_seg)
cmake_minimum_required(VERSION 2.8)
cmake_policy(VERSION 2.8)
include(CPack)
include(ExternalProject)


SET(NiftySeg_VERSION_MAJOR 0)
SET(NiftySeg_VERSION_MINOR 9)
SET(CMAKE_CXX_FLAGS -Wall)
IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")
 MARK_AS_ADVANCED(FORCE CMAKE_BACKWARDS_COMPATIBILITY)
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")
 MARK_AS_ADVANCED(CLEAR CMAKE_BACKWARDS_COMPATIBILITY)
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")

#-----------------------------------------------------------------------------
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)
add_definitions(-Wall -pedantic)

#-----------------------------------------------------------------------------
SET(ZLIB "z")
#-----------------------------------------------------------------------------


OPTION(BUILD_ALL "To build the LoAd segmentation tool" ON)
OPTION(BUILD_SHARED "Build the libraries as shared." OFF)
OPTION(USE_OPENMP "To use openMP for multi-CPU processing" OFF)
OPTION(INSTALL_NIFTYREG "To download and install NiftyReg. Internet connection is required." OFF)
OPTION(INSTALL_PRIORS "Install prior images. Turn this off if you are really short on disk space" ON)
SET(INSTALL_PRIORS_DIRECTORY ${CMAKE_INSTALL_PREFIX}/priors CACHE PATH "Folder in which to install the prior images.")
#-----------------------------------------------------------------------------
IF(BUILD_SHARED)
  SET(LIBRARY_TYPE SHARED)
ELSE(BUILD_SHARED)
  SET(LIBRARY_TYPE STATIC)
ENDIF(BUILD_SHARED)

IF(USE_OPENMP)
SET(CMAKE_CXX_FLAGS -fopenmp)
ENDIF(USE_OPENMP)

IF(BUILD_SHARED)
    SET(LIBRARY_TYPE SHARED)
ELSE(BUILD_SHARED)
       SET(LIBRARY_TYPE STATIC)
ENDIF(BUILD_SHARED)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/seg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nifti)

#-----------------------------------------------------------------------------

SUBDIRS(zlib)
SUBDIRS(nifti)
SUBDIRS(seg-lib)
SUBDIRS(seg-apps)


# Install scripts
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/LoAd_brainonly.sh ${CMAKE_BINARY_DIR}/LoAd_brainonly.sh @ONLY )
INSTALL(PROGRAMS ${CMAKE_BINARY_DIR}/LoAd_brainonly.sh DESTINATION ${CMAKE_INSTALL_PREFIX}/bin COMPONENT applications)

# Install priors
IF(INSTALL_PRIORS)
  FILE(GLOB PRIOR_FILES priors/*.nii.gz)
  INSTALL(PROGRAMS ${PRIOR_FILES} DESTINATION ${INSTALL_PRIORS_DIRECTORY} COMPONENT applications)
ENDIF(INSTALL_PRIORS)

IF(USE_OPENMP)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
ENDIF(USE_OPENMP)

#-----------------------------------------------------------------------------

IF(INSTALL_NIFTYREG)
ExternalProject_Add(niftyreg
            SOURCE_DIR ${PROJECT_BINARY_DIR}/niftyreg
            PREFIX ${PROJECT_BINARY_DIR}/niftyreg
            STAMP_DIR ${PROJECT_BINARY_DIR}/niftyreg/stamps
            TMP_DIR ${PROJECT_BINARY_DIR}/niftyreg/stamps
            DOWNLOAD_DIR ${PROJECT_BINARY_DIR}/niftyreg
            BINARY_DIR ${PROJECT_BINARY_DIR}/niftyreg/build
            SVN_REPOSITORY https://niftyreg.svn.sourceforge.net/svnroot/niftyreg/trunk/nifty_reg
            TIMEOUT 100
            UPDATE_COMMAND ""
            CMAKE_ARGS -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX:PATH=${CMAKE_INSTALL_PREFIX} -DUSE_SSE=1
            BUILD_COMMAND ${CMAKE_BUILD_TOOL}
            INSTALL_COMMAND ${CMAKE_BUILD_TOOL} install
            )
ENDIF(INSTALL_NIFTYREG)
#-----------------------------------------------------------------------------

IF(DO_PACKAGE)
        IF(UNIX)
            IF(APPLE)
                SET(CPACK_GENERATOR "TZ;TGZ;PackageMaker")
            ELSE(APPLE)
                SET(CPACK_GENERATOR "TZ;TGZ;DEB")
            ENDIF(APPLE)
        ENDIF(UNIX)

        SET(CPACK_PACKAGE_NAME "NiftySeg")
        SET(CPACK_PACKAGE_VENDOR "CMIC-UCL")
        SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "NiftySeg - Software for Image Segmentation")
        SET(CPACK_PACKAGE_VERSION_MAJOR ${NiftySeg_VERSION_MAJOR})
        SET(CPACK_PACKAGE_VERSION_MINOR ${NiftySeg_VERSION_MINOR})
        SET(CPACK_PACKAGE_INSTALL_DIRECTORY "NiftySeg")
        SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
        SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE.txt")
        SET(CPACK_PACKAGE_ICON "${CMake_SOURCE_DIR}/images\\\\NiftySeg_logo.png")
        SET(CPACK_NSIS_MUI_ICON "${CMake_SOURCE_DIR}/images\\\\NiftySeg_logo.png")
        SET(CPACK_NSIS_MUI_UNIICON "${CMake_SOURCE_DIR}/images\\\\NiftySeg_logo.png")

        IF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
          SET (CMAKE_SYSTEM_PROCESSOR "x86")
        ENDIF (CMAKE_SYSTEM_PROCESSOR MATCHES "unknown")
        IF(NOT DEFINED CPACK_SYSTEM_NAME)
          SET(CPACK_SYSTEM_NAME ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})
        ENDIF(NOT DEFINED CPACK_SYSTEM_NAME)
        IF(${CPACK_SYSTEM_NAME} MATCHES Windows)
          IF(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win64-${CMAKE_SYSTEM_PROCESSOR})
          ELSE(CMAKE_CL_64)
            SET(CPACK_SYSTEM_NAME Win32-${CMAKE_SYSTEM_PROCESSOR})
          ENDIF(CMAKE_CL_64)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Windows)

        IF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)
          list(LENGTH CMAKE_OSX_ARCHITECTURES _length)
          IF(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-Universal)
          ELSE(_length GREATER 1)
            SET(CPACK_SYSTEM_NAME Darwin-${CMAKE_OSX_ARCHITECTURES})
          ENDIF(_length GREATER 1)
        ENDIF(${CPACK_SYSTEM_NAME} MATCHES Darwin AND CMAKE_OSX_ARCHITECTURES)

        #Deb
        SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "M. Jorge Cardoso, manuel.cardoso@ucl.ac.uk")
        SET(CPACK_DEBIAN_PACKAGE_DEPENDS "niftyreg (>= 1.2.1)")

        #Nsis
        SET(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\NiftySeg.exe")
        SET(CPACK_NSIS_DISPLAY_NAME "${CPACK_PACKAGE_INSTALL_DIRECTORY} NiftySeg")
        SET(CPACK_NSIS_HELP_LINK "http:\\\\\\\\cmic.cs.ucl.ac.uk")
        SET(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\cmic.cs.ucl.ac.uk")
        SET(CPACK_NSIS_CONTACT "manue.cardoso@ucl.ac.uk")
        SET(CPACK_NSIS_MODIFY_PATH ON)

        #Apple drag-and-drop
        IF(APPLE)
                SET(CPACK_BINARY_TBZ2 OFF)
                SET(CPACK_BINARY_DRAGNDROP ON)
                SET(CPACK_BINARY_PACKAGEMAKER ON)
                SET(CPACK_BINARY_STGZ OFF)
        ENDIF(APPLE)

        INCLUDE(CPack)
ENDIF(DO_PACKAGE)

#-----------------------------------------------------------------------------

