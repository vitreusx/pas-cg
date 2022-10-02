cmake_minimum_required(VERSION 3.17)
enable_language(Fortran)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cg-cxx
    COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:cg>" ${CMAKE_CURRENT_BINARY_DIR}/cg-cxx
    DEPENDS cg
)

add_custom_target(copy-cg-cxx DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/cg-cxx")

add_custom_target(parity)
add_dependencies(parity copy-cg-cxx copy-tests-dir)

set(FortranVariants pho5 pho6)
foreach (Variant IN LISTS FortranVariants)
    add_executable(cg-fort-${Variant})
    target_sources(cg-fort-${Variant} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/cg_${Variant}.f)

    target_compile_options(cg-fort-${Variant}
        PRIVATE
        -std=legacy;-mcmodel=large
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Release>:-Ofast;-march=native>)

    add_dependencies(parity cg-fort-${Variant})
endforeach ()

file(GLOB Tests CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/*)

set(TestNames)
foreach (Test IN LISTS Tests)
    get_filename_component(TestName "${Test}" NAME)
    list(APPEND TestNames "${TestName}")
endforeach ()

file(GLOB_RECURSE TestFiles CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/**)

add_custom_target(remove-tests-dir
    COMMAND ${CMAKE_COMMAND} -E rm -rf
    "${CMAKE_CURRENT_BINARY_DIR}/tests")

add_custom_target(copy-tests-dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    "${CMAKE_CURRENT_BINARY_DIR}/tests"
    DEPENDS ${TestFiles} remove-tests-dir)

foreach (TestName IN LISTS TestNames)
    add_test(
        NAME "${TestName}"
        COMMAND python3 "${CMAKE_CURRENT_SOURCE_DIR}/parity.py"
        --bin "${CMAKE_CURRENT_BINARY_DIR}"
        --test "${CMAKE_CURRENT_BINARY_DIR}/tests/${TestName}"
        --source "${CMAKE_CURRENT_SOURCE_DIR}/tests/${TestName}"
    )
    set_tests_properties("${TestName}" PROPERTIES DEPENDS parity)
endforeach ()
