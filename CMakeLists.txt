PROJECT(nifty_seg)

cmake_minimum_required(VERSION 2.6.4)

IF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")
 MARK_AS_ADVANCED(FORCE CMAKE_BACKWARDS_COMPATIBILITY)
ELSE("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")
 MARK_AS_ADVANCED(CLEAR CMAKE_BACKWARDS_COMPATIBILITY)
ENDIF("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" MATCHES "^2\\.6\\.4$")

#-----------------------------------------------------------------------------
if(COMMAND cmake_policy)
      cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#-----------------------------------------------------------------------------
SET(ZLIB "z")
#-----------------------------------------------------------------------------


OPTION(BUILD_LOAD "To build the LoAd segmentation tool" ON)
OPTION(USE_SSE "to use SSE computation in some case" OFF)
OPTION(INSTALL_PRIORS "Install prior images. Turn this off if you are really short on disk space" ON)

SET(INSTALL_PRIORS_DIRECTORY ${CMAKE_INSTALL_PREFIX}/priors CACHE PATH "Folder in which to install the prior images.")
#-----------------------------------------------------------------------------

ADD_DEFINITIONS(-fPIC)

IF(USE_SSE)
    ADD_DEFINITIONS(-D_USE_SSE)
ENDIF(USE_SSE)

IF(BUILD_NR_DEV)
    ADD_DEFINITIONS(-D_NR_DEV)
ENDIF(BUILD_NR_DEV)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/reg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/seg-lib)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/nifti)

IF(USE_SSE)
	ADD_DEFINITIONS(-D_USE_SSE)
ENDIF(USE_SSE)

#-----------------------------------------------------------------------------

SUBDIRS(zlib)
SUBDIRS(nifti)
SUBDIRS(reg-lib)
SUBDIRS(reg-apps)
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

#-----------------------------------------------------------------------------
