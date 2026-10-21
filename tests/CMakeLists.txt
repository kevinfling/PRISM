add_executable(test_core test_core.c)
add_executable(test_contracts test_contracts.c)
foreach(t test_core test_contracts)
    target_link_libraries(${t} PRIVATE prism::prism)
    prism_checks(${t})
endforeach()
add_test(NAME reference COMMAND test_core reference)
add_test(NAME contracts COMMAND test_contracts)
set_tests_properties(reference contracts PROPERTIES LABELS prism TIMEOUT 240)
if(PRISM_ENABLE_FFTW)
    add_test(NAME fftw COMMAND test_core fftw)
    set_tests_properties(fftw PROPERTIES LABELS prism TIMEOUT 120)
endif()
add_executable(test_properties test_properties.c)
target_link_libraries(test_properties PRIVATE prism::prism)
prism_checks(test_properties)
add_test(NAME properties COMMAND test_properties)
set_tests_properties(properties PROPERTIES LABELS prism TIMEOUT 120)
if(UNIX)
    find_package(Threads REQUIRED)
    add_executable(test_concurrency test_concurrency.c)
    target_link_libraries(test_concurrency PRIVATE prism::prism Threads::Threads)
    prism_checks(test_concurrency)
    add_test(NAME concurrency COMMAND test_concurrency)
    set_tests_properties(concurrency PROPERTIES LABELS prism TIMEOUT 120)
endif()
add_test(NAME installed_consumer COMMAND ${CMAKE_COMMAND}
    -DPRISM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPRISM_BINARY_DIR=${PROJECT_BINARY_DIR}
    -DPRISM_C_COMPILER=${CMAKE_C_COMPILER}
    -DPRISM_CONFIG=$<CONFIG>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/install_test.cmake)
set_tests_properties(installed_consumer PROPERTIES LABELS prism TIMEOUT 120 RUN_SERIAL TRUE)
add_executable(test_math test_math.c)
target_link_libraries(test_math PRIVATE prism::prism)
prism_checks(test_math)
add_test(NAME mathematical_properties COMMAND test_math)
set_tests_properties(mathematical_properties PROPERTIES LABELS prism TIMEOUT 120)
add_test(NAME subproject_consumer COMMAND ${CMAKE_COMMAND}
    -DPRISM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
    -DPRISM_BINARY_DIR=${PROJECT_BINARY_DIR}
    -DPRISM_C_COMPILER=${CMAKE_C_COMPILER}
    -DPRISM_FFTW=${PRISM_ENABLE_FFTW}
    -DPRISM_SIGMA=${PRISM_ENABLE_SIGMA}
    -DPRISM_SIGMA_SOURCE_DIR=${PRISM_SIGMA_SOURCE_DIR}
    -DPRISM_CONFIG=$<CONFIG>
    -P ${CMAKE_CURRENT_SOURCE_DIR}/subproject_test.cmake)
set_tests_properties(subproject_consumer PROPERTIES LABELS prism TIMEOUT 120 RUN_SERIAL TRUE)
find_package(Python3 QUIET COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(NAME documentation COMMAND ${Python3_EXECUTABLE} ${PROJECT_SOURCE_DIR}/tools/check_docs.py)
    set_tests_properties(documentation PROPERTIES LABELS prism TIMEOUT 30)
endif()
find_package(Git QUIET)
if(Git_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    add_test(NAME clean_source COMMAND ${CMAKE_COMMAND}
        -DPRISM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DPRISM_BINARY_DIR=${PROJECT_BINARY_DIR}
        -DPRISM_C_COMPILER=${CMAKE_C_COMPILER}
        -DPRISM_FFTW=${PRISM_ENABLE_FFTW}
        -DPRISM_SIGMA=${PRISM_ENABLE_SIGMA}
        -DPRISM_SIGMA_SOURCE_DIR=${PRISM_SIGMA_SOURCE_DIR}
        -DPRISM_GIT=${GIT_EXECUTABLE}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/clean_source.cmake)
    set_tests_properties(clean_source PROPERTIES LABELS prism TIMEOUT 180 RUN_SERIAL TRUE)
endif()

add_executable(test_noise test_noise.c)
target_link_libraries(test_noise PRIVATE prism::prism)
prism_checks(test_noise)
add_test(NAME stochastic_covariance COMMAND test_noise)
set_tests_properties(stochastic_covariance PROPERTIES LABELS prism TIMEOUT 120)

if(UNIX)
    add_executable(test_noise_concurrency test_noise_concurrency.c)
    target_link_libraries(test_noise_concurrency PRIVATE prism::prism Threads::Threads)
    prism_checks(test_noise_concurrency)
    add_test(NAME stochastic_concurrency COMMAND test_noise_concurrency)
    set_tests_properties(stochastic_concurrency PROPERTIES LABELS prism TIMEOUT 120)
endif()
add_executable(test_innovation test_innovation.c)
target_link_libraries(test_innovation PRIVATE prism::prism)
prism_checks(test_innovation)
add_test(NAME innovation_providers COMMAND test_innovation)
set_tests_properties(innovation_providers PROPERTIES LABELS prism TIMEOUT 120)
add_executable(test_noise_extended test_noise_extended.c)
target_link_libraries(test_noise_extended PRIVATE prism::prism)
prism_checks(test_noise_extended)
add_test(NAME covariance_awkward_lengths COMMAND test_noise_extended)
set_tests_properties(covariance_awkward_lengths PROPERTIES LABELS prism TIMEOUT 240)
if(PRISM_ENABLE_SIGMA)
    add_executable(test_sigma_upstream test_sigma_upstream.c)
    target_include_directories(test_sigma_upstream SYSTEM PRIVATE "${PRISM_SIGMA_SOURCE_DIR}/include")
    target_link_libraries(test_sigma_upstream PRIVATE prism::prism)
    prism_checks(test_sigma_upstream)
    add_test(NAME sigma_upstream_audit COMMAND test_sigma_upstream)
    set_tests_properties(sigma_upstream_audit PROPERTIES LABELS prism TIMEOUT 30)
    add_test(NAME sigma_dependency_pin COMMAND ${CMAKE_COMMAND}
        -DPRISM_SOURCE_DIR=${PROJECT_SOURCE_DIR}
        -DPRISM_BINARY_DIR=${PROJECT_BINARY_DIR}
        -DPRISM_C_COMPILER=${CMAKE_C_COMPILER}
        -DPRISM_SIGMA_SOURCE_DIR=${PRISM_SIGMA_SOURCE_DIR}
        -P ${CMAKE_CURRENT_SOURCE_DIR}/sigma_pin_test.cmake)
    set_tests_properties(sigma_dependency_pin PROPERTIES LABELS prism TIMEOUT 60)
endif()
