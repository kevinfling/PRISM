include("${PRISM_SOURCE_DIR}/tests/commands.cmake")
execute_process(COMMAND "${PRISM_GIT}" -C "${PRISM_SOURCE_DIR}" ls-files --cached --others --exclude-standard
    RESULT_VARIABLE status OUTPUT_VARIABLE files ERROR_VARIABLE error)
if(NOT status EQUAL 0)
    message(FATAL_ERROR "Cannot enumerate reviewable source files: ${error}")
endif()
set(snapshot "${PRISM_BINARY_DIR}/clean-source-check")
file(REMOVE_RECURSE "${snapshot}")
string(REPLACE "\n" ";" files "${files}")
foreach(path IN LISTS files)
    if(NOT path STREQUAL "")
        get_filename_component(directory "${path}" DIRECTORY)
        file(COPY "${PRISM_SOURCE_DIR}/${path}" DESTINATION "${snapshot}/source/${directory}")
    endif()
endforeach()
check_command("${CMAKE_COMMAND}" -S "${snapshot}/source" -B "${snapshot}/build"
    "-DCMAKE_C_COMPILER=${PRISM_C_COMPILER}" -DCMAKE_BUILD_TYPE=Release
    "-DPRISM_ENABLE_FFTW=${PRISM_FFTW}" -DPRISM_WARNINGS_AS_ERRORS=ON)
check_command("${CMAKE_COMMAND}" --build "${snapshot}/build" --parallel 2)
check_command("${CMAKE_CTEST_COMMAND}" --test-dir "${snapshot}/build"
    -R "^(contracts|properties|documentation|installed_consumer)$" --output-on-failure)
