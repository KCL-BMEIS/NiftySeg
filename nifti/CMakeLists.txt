#-----------------------------------------------------------------------------
option(USE_NII_NAN "To enable NaN support with nifti images" ON)
mark_as_advanced(FORCE USE_NII_NAN)
if(USE_NII_NAN)
    add_definitions(-DUSE_NII_NAN)
endif(USE_NII_NAN)
add_definitions(-DHAVE_ZLIB)

set(NAME _seg_nifti)

# On Windows none of this is exported, so leaving it on STATIC
# otherwise there will be no .lib file generated
if(WIN32)
  SET(LIBRARY_TYPE STATIC)
endif()
 
add_library(${NAME}  ${LIBRARY_TYPE} nifti1.h nifti1_io.c nifti1_io.h znzlib.c znzlib.h)
target_link_libraries(${NAME} z)
install(TARGETS ${NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(FILES nifti1_io.h znzlib.h nifti1.h DESTINATION include COMPONENT Development)
#add_dependencies(${NAME} Eigen)
#-----------------------------------------------------------------------------
