execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzvf ${TEST_DIR}/schwarzschild.tgz WORKING_DIRECTORY ${TEST_DIR} RESULT_VARIABLE tar_failed)

if(tar_failed)
	message("File Extract Failed")
endif(tar_failed)

execute_process(COMMAND ${EXE_DIR}/ray_trace_ellipse --glellipse 0,128,0,128 --usesurfacemassdensity schwarzschild_mass.nc --parameters schwarzschild_parameters.xml --newlens test.lens OUTPUT_FILE ray_trace_test.out WORKING_DIRECTORY ${TEST_DIR}/schwarzschild RESULT_VARIABLE test_failed)

if(test_failed)
  message("Ray Trace Test Failed")
endif(test_failed)

