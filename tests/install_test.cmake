include("${PRISM_SOURCE_DIR}/tests/commands.cmake")
set(prefix "${PRISM_BINARY_DIR}/consumer-prefix")
set(consumer "${PRISM_BINARY_DIR}/consumer-build")
check_command("${CMAKE_COMMAND}" --install "${PRISM_BINARY_DIR}" --prefix "${prefix}" --config "${PRISM_CONFIG}")
check_command("${CMAKE_COMMAND}" -S "${PRISM_SOURCE_DIR}/tests/consumer" -B "${consumer}"
    "-DCMAKE_PREFIX_PATH=${prefix}" "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}")
check_command("${CMAKE_COMMAND}" --build "${consumer}" --config "${PRISM_CONFIG}")
check_command("${CMAKE_CTEST_COMMAND}" --test-dir "${consumer}" -C "${PRISM_CONFIG}" --output-on-failure)
# Compile and run the complete C example printed in the README against the installed library.
file(READ "${PRISM_SOURCE_DIR}/README.md" readme)
string(FIND "${readme}" "```c\n" begin)
if(begin EQUAL -1)
    message(FATAL_ERROR "README has no C example")
endif()
math(EXPR begin "${begin}+5")
string(SUBSTRING "${readme}" ${begin} -1 rest)
string(FIND "${rest}" "```" end)
if(end EQUAL -1)
    message(FATAL_ERROR "README C example has no closing fence")
endif()
string(SUBSTRING "${rest}" 0 ${end} code)
set(documented "${PRISM_BINARY_DIR}/documented-example")
file(MAKE_DIRECTORY "${documented}")
file(WRITE "${documented}/main.c" "${code}")
check_command("${CMAKE_COMMAND}" -S "${PRISM_SOURCE_DIR}/tests/consumer" -B "${documented}/build"
    "-DCMAKE_PREFIX_PATH=${prefix}" "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}"
    "-DPRISM_EXAMPLE_SOURCE=${documented}/main.c")
check_command("${CMAKE_COMMAND}" --build "${documented}/build" --config "${PRISM_CONFIG}")
check_command("${CMAKE_CTEST_COMMAND}" --test-dir "${documented}/build" -C "${PRISM_CONFIG}" --output-on-failure)
# The 0.1 package must not advertise compatibility with the 0.2 API.
execute_process(COMMAND "${CMAKE_COMMAND}" -S "${PRISM_SOURCE_DIR}/tests/consumer"
    -B "${PRISM_BINARY_DIR}/incompatible-consumer" "-DCMAKE_PREFIX_PATH=${prefix}"
    "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}" -DPRISM_REQUESTED_VERSION=0.2
    RESULT_VARIABLE incompatible OUTPUT_VARIABLE output ERROR_VARIABLE error)
if(incompatible EQUAL 0 OR NOT error MATCHES "version")
    message(FATAL_ERROR "Expected incompatible package-version rejection: ${output}\n${error}")
endif()
