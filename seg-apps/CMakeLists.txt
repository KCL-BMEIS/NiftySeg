#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)

    ADD_EXECUTABLE(seg_Stats seg_Stats.cpp)
        TARGET_LINK_LIBRARIES(seg_Stats _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_Stats)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_LabFusion seg_LabFusion.cpp)
        TARGET_LINK_LIBRARIES(seg_LabFusion _seg_LabFusion _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LabFusion)

    ADD_EXECUTABLE(seg_maths seg_maths.cpp)
        TARGET_LINK_LIBRARIES(seg_maths _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_maths)


#-----------------------------------------------------------------------------
