#-----------------------------------------------------------------------------
ADD_DEFINITIONS(-DHAVE_ZLIB)
SET(NAME _seg_nifti)
ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} nifti1.h nifti1_io.c nifti1_io.h znzlib.c znzlib.h)
TARGET_LINK_LIBRARIES(${NAME} z)
INSTALL(TARGETS ${NAME}
	RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
	ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
INSTALL(FILES znzlib.h DESTINATION include)
INSTALL(FILES nifti1_io.h DESTINATION include)
INSTALL(FILES nifti1.h DESTINATION include)
#-----------------------------------------------------------------------------
