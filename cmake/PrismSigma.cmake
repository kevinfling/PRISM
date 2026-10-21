# Only these two pinned upstream headers are compiled. Avoid importing SIGMA's
# directory-wide optimization/fast-math/ISA flags or optional executables.
set(PRISM_SIGMA_REVISION 6e9a1ef526365474738d80d4eb71fd945a13abd6)
if(NOT PRISM_SIGMA_SOURCE_DIR)
    message(FATAL_ERROR "PRISM_ENABLE_SIGMA requires PRISM_SIGMA_SOURCE_DIR; see docs/SIGMA.md")
endif()
set(sigma_files include/sigma/core.h include/sigma/rng.h LICENSE)
set(sigma_hashes
    84860d2b17deedd77a13878ea5bfc99ebd8a06429c01c4640658cac1f81e44f4
    8dd7b30eb892e4dce51338e040f9dc56b250d936ff21254c26467e8dadba3f2a
    8d9f5e03c2637169e3a0c6b7d45cf0a9122991884d7218879bec5b1fe6a8dd6f)
foreach(i RANGE 0 2)
    list(GET sigma_files ${i} path)
    list(GET sigma_hashes ${i} expected)
    if(NOT EXISTS "${PRISM_SIGMA_SOURCE_DIR}/${path}")
        message(FATAL_ERROR "Pinned SIGMA file missing: ${path}")
    endif()
    file(SHA256 "${PRISM_SIGMA_SOURCE_DIR}/${path}" actual)
    if(NOT actual STREQUAL expected)
        message(FATAL_ERROR "Pinned SIGMA hash mismatch: ${path} (expected revision ${PRISM_SIGMA_REVISION})")
    endif()
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${PRISM_SIGMA_SOURCE_DIR}/${path}")
endforeach()
target_include_directories(prism SYSTEM PRIVATE "${PRISM_SIGMA_SOURCE_DIR}/include")
target_compile_definitions(prism PRIVATE PRISM_HAVE_SIGMA=1 PRISM_SIGMA_REVISION="${PRISM_SIGMA_REVISION}")
install(FILES "${PRISM_SIGMA_SOURCE_DIR}/LICENSE"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/licenses/PRISM RENAME SIGMA-MIT.txt)
