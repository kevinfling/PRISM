include("${PRISM_SOURCE_DIR}/tests/commands.cmake")
set(consumer "${PRISM_BINARY_DIR}/subproject-consumer")
check_command("${CMAKE_COMMAND}" -S "${PRISM_SOURCE_DIR}/tests/subproject" -B "${consumer}"
    "-DPRISM_SOURCE_DIR=${PRISM_SOURCE_DIR}" "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}"
    "-DPRISM_ENABLE_FFTW=${PRISM_FFTW}")
check_command("${CMAKE_COMMAND}" --build "${consumer}" --config "${PRISM_CONFIG}")
check_command("${CMAKE_CTEST_COMMAND}" --test-dir "${consumer}" -C "${PRISM_CONFIG}" --output-on-failure)
