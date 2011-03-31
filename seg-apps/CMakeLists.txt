#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)

    ADD_EXECUTABLE(seg_Test seg_Test.cpp)
        TARGET_LINK_LIBRARIES(seg_Test _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_Test)

    ADD_EXECUTABLE(seg_copy_header seg_copy_header.cpp)
        TARGET_LINK_LIBRARIES(seg_copy_header seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_copy_header)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_STAPLE seg_STAPLE.cpp)
        TARGET_LINK_LIBRARIES(seg_STAPLE _seg_STAPLE _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_STAPLE)

    ADD_EXECUTABLE(seg_neonate_PV seg_neonate_PV.cpp)
        TARGET_LINK_LIBRARIES(seg_neonate_PV _seg_EM _seg_tools _seg_Topo _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_neonate_PV)

    ADD_EXECUTABLE(seg_Split seg_Split.cpp)
        TARGET_LINK_LIBRARIES(seg_Split _seg_LoAd _seg_tools _seg_MRF _seg_FMM _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_Split)

    ADD_EXECUTABLE(seg_CT_Lung seg_CT_Lung.cpp)
        TARGET_LINK_LIBRARIES(seg_CT_Lung _seg_CT_lungs _seg_Topo _seg_MRF _seg_tools _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_CT_Lung)

    ADD_EXECUTABLE(seg_CT_Colon seg_CT_Colon.cpp)
        TARGET_LINK_LIBRARIES(seg_CT_Colon _seg_CT_colon _seg_Topo _seg_tools _seg_MRF _seg_BiasCorrection seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_CT_Colon)

#-----------------------------------------------------------------------------
