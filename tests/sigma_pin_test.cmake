include("${PRISM_SOURCE_DIR}/tests/commands.cmake")
set(copy "${PRISM_BINARY_DIR}/sigma-pin-check/source")
file(MAKE_DIRECTORY "${copy}/include/sigma")
foreach(header core.h rng.h)
    configure_file("${PRISM_SIGMA_SOURCE_DIR}/include/sigma/${header}"
        "${copy}/include/sigma/${header}" COPYONLY)
endforeach()
configure_file("${PRISM_SIGMA_SOURCE_DIR}/LICENSE" "${copy}/LICENSE" COPYONLY)
file(APPEND "${copy}/include/sigma/rng.h" "\n/* changed source: pin must reject */\n")
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${PRISM_SOURCE_DIR}"
    -B "${PRISM_BINARY_DIR}/sigma-pin-check/build" "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}"
    -DPRISM_ENABLE_FFTW=OFF -DPRISM_BUILD_TESTS=OFF -DPRISM_BUILD_EXAMPLES=OFF
    -DPRISM_BUILD_BENCHMARKS=OFF -DPRISM_BUILD_EXPERIMENTS=OFF
    -DPRISM_ENABLE_SIGMA=ON "-DPRISM_SIGMA_SOURCE_DIR=${copy}"
    RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(status EQUAL 0 OR NOT error MATCHES "Pinned SIGMA hash mismatch")
    message(FATAL_ERROR "Expected changed-header rejection: ${output}\n${error}")
endif()
