#-----------------------------------------------------------------------------

    ADD_EXECUTABLE(seg_LoAd seg_LoAd.cpp)
        TARGET_LINK_LIBRARIES(seg_LoAd _seg_LoAd _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LoAd)
   #add_dependencies(seg_LoAd Eigen)

    ADD_EXECUTABLE(seg_stats seg_Stats.cpp)
        TARGET_LINK_LIBRARIES(seg_stats _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_stats)

    ADD_EXECUTABLE(seg_EM seg_EM.cpp)
        TARGET_LINK_LIBRARIES(seg_EM _seg_EM _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_EM)

    ADD_EXECUTABLE(seg_LabFusion seg_LabFusion.cpp)
        TARGET_LINK_LIBRARIES(seg_LabFusion _seg_LabFusion _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_LabFusion)

    ADD_EXECUTABLE(seg_maths seg_maths.cpp)
        TARGET_LINK_LIBRARIES(seg_maths _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_maths)
   #add_dependencies(seg_maths Eigen)

    ADD_EXECUTABLE(seg_CalcTopNCC seg_CalcTopNCC.cpp)
        TARGET_LINK_LIBRARIES(seg_CalcTopNCC _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_CalcTopNCC)

    ADD_EXECUTABLE(seg_PatchMatch seg_PatchMatch.cpp)
	TARGET_LINK_LIBRARIES(seg_PatchMatch _seg_PatchMatch _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_PatchMatch)

    ADD_EXECUTABLE(seg_FillLesions seg_FillLesions.cpp)
	TARGET_LINK_LIBRARIES(seg_FillLesions _seg_fill_lesions _seg_fill_lesions_other _seg_tools _seg_FMM _seg_nifti ${ZLIB})
    INSTALL_TARGETS(/bin seg_FillLesions)


#-----------------------------------------------------------------------------
set(MODULE_LIST
  seg_LoAd
  seg_stats
  seg_EM
  seg_LabFusion
  seg_maths
  seg_CalcTopNCC
  seg_PatchMatch
  seg_FillLesions
)

foreach(MODULE_NAME ${MODULE_LIST})
  install(TARGETS ${MODULE_NAME}
    RUNTIME DESTINATION bin COMPONENT Runtime
    LIBRARY DESTINATION lib COMPONENT Runtime
    ARCHIVE DESTINATION lib COMPONENT Runtime
  )
endforeach(MODULE_NAME)
