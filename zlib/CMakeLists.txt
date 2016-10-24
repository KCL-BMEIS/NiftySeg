#-----------------------------------------------------------------------------
SET(NAME z)

# On Windows none of this is exported, so leaving it on STATIC
# otherwise there will be no .lib file generated
if(WIN32)
  SET(LIBRARY_TYPE STATIC)
endif()

ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} adler32.c compress.c crc32.c gzio.c uncompr.c deflate.c trees.c zutil.c inflate.c infback.c inftrees.c inffast.c)
INSTALL(TARGETS ${NAME}
	RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
	LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
	ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
#add_dependencies(${NAME} Eigen)

#-----------------------------------------------------------------------------
