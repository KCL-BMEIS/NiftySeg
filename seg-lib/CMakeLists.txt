#-----------------------------------------------------------------------------

    SET(NAME _seg_tools)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_matrix.h)
    TARGET_LINK_LIBRARIES(${NAME} _seg_FMM _seg_nifti)
    INSTALL(TARGETS ${NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h _seg_matrix.h DESTINATION include)
    ##add_dependencies(_seg_tools Eigen)

    SET(NAME _seg_FMM)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
    INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h DESTINATION include)
    #add_dependencies(_seg_FMM Eigen)

    SET(NAME _seg_LoAd)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
    TARGET_LINK_LIBRARIES(${NAME} _seg_nifti _seg_tools _seg_FMM)
    INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h DESTINATION include)
    #add_dependencies(_seg_LoAd Eigen)

    SET(NAME _seg_EM)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
    TARGET_LINK_LIBRARIES(${NAME} _seg_nifti _seg_tools _seg_FMM)
    INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h DESTINATION include)

    SET(NAME _seg_LabFusion)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_tools)
    TARGET_LINK_LIBRARIES(${NAME} _seg_nifti _seg_tools _seg_FMM)
    INSTALL(TARGETS ${NAME}
            RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
            LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
            ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h DESTINATION include)
    
    set(NAME _seg_PatchMatch)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_matrix.h _seg_PatchMatchResult.h _seg_PatchMatchResult.cpp)
    TARGET_LINK_LIBRARIES(${NAME} _seg_tools _seg_FMM _seg_nifti)
    INSTALL(TARGETS ${NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    )
    INSTALL(FILES ${NAME}.h _seg_matrix.h DESTINATION include)
    
    set(NAME _seg_fill_lesions)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_matrix.h)
    TARGET_LINK_LIBRARIES(${NAME} _seg_tools _seg_FMM _seg_nifti)
    INSTALL(FILES ${NAME}.h _seg_matrix.h DESTINATION include)
    
    set(NAME _seg_fill_lesions_other)
    ADD_LIBRARY(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_matrix.h)
    TARGET_LINK_LIBRARIES(${NAME} _seg_tools _seg_FMM _seg_nifti)
    INSTALL(FILES ${NAME}.h _seg_matrix.h DESTINATION include)

INSTALL(FILES _seg_common.h NiftySegWinExportHeader.h DESTINATION include)



