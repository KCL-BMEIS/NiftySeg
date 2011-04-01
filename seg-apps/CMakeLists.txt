#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_STAPLE seg_STAPLE.cpp)
        TARGET_LINK_LIBRARIES(seg_STAPLE _seg_STAPLE _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_STAPLE)

#-----------------------------------------------------------------------------
