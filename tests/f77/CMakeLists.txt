cmake_minimum_required(VERSION 3.14)
project(cg-f77)
enable_language(Fortran)

add_executable(cg-f77)
target_sources(cg-f77 PRIVATE cg.f)

target_compile_options(cg-f77
    PRIVATE
    -std=legacy
    $<$<CONFIG:Debug>:-g;-O1>
    $<$<CONFIG:Release>:-Ofast;-march=native>)

#find_package(OpenMP REQUIRED)
#target_link_libraries(cg-f77 PRIVATE $<$<CONFIG:Release>:OpenMP::OpenMP_Fortran>)

file(COPY data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_custom_command(
    TARGET cg-f77
    POST_BUILD
    COMMAND objdump -drlwC -Mintel cg-f77 >cg-f77.s
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})