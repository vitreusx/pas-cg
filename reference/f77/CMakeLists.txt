cmake_minimum_required(VERSION 3.14)
project(cg-f77)
enable_language(Fortran)

file(GLOB_RECURSE DataFiles CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/data/**")

#add_custom_target(remove-data-dir
#    COMMAND ${CMAKE_COMMAND} -E rm -rf
#    "${CMAKE_CURRENT_BINARY_DIR}/data")

add_custom_target(copy-data-dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_BINARY_DIR}/data"
    DEPENDS ${DataFiles})

foreach (target cg_pho5)
    add_executable(${target})
    target_sources(${target} PRIVATE ${target}.f)

    target_compile_options(${target}
        PRIVATE
        -std=legacy;-mcmodel=large
        $<$<CONFIG:Debug>:-g>
        $<$<CONFIG:Release>:-Ofast;-march=native>)

    #find_package(OpenMP REQUIRED)
    #target_link_libraries(cg-f77 PRIVATE $<$<CONFIG:Release>:OpenMP::OpenMP_Fortran>)

    add_dependencies(${target} copy-data-dir)

    add_custom_command(
        TARGET ${target}
        POST_BUILD
        COMMAND objdump -drlwC -Mintel ${target} >${target}.s
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endforeach ()