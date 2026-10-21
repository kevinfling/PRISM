add_executable(prism_example decompose.c)
target_link_libraries(prism_example PRIVATE prism::prism)
prism_checks(prism_example)
