#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_Split seg_Split.cpp)
        TARGET_LINK_LIBRARIES(seg_Split seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_Split)

    ADD_EXECUTABLE(seg_ProbLabFusion seg_ProbLabFusion.cpp)
        TARGET_LINK_LIBRARIES(seg_ProbLabFusion _seg_ProbLabFusion _seg_Topo _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_ProbLabFusion)

    ADD_EXECUTABLE(seg_MultiLabFusion seg_MultiLabFusion.cpp)
        TARGET_LINK_LIBRARIES(seg_MultiLabFusion _seg_MultiLabFusion _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_MultiLabFusion)

#-----------------------------------------------------------------------------
