option(LEGACY_MODE "Keep the results as close to Fortran version(s)
as possible." ON)

option(MIXED_PRECISION "Whether to use floats for computing the
forces and doubles for integration, instead of using doubles
throughout." OFF)

function(message)
    if (NOT MESSAGE_QUIET)
        _message(${ARGN})
    endif ()
endfunction()

add_executable(cg)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/src/cg/simul/default.yml.h
    COMMAND xxd -i default.yml >src/cg/simul/default.yml.h
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/default.yml
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

add_custom_target(defaults ALL
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/cg/simul/default.yml.h)

add_dependencies(cg defaults)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp src/*.hpp src/*.h src/*.inl src/*.cu src/*.cuh)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl include/*.cuh)

target_sources(cg
    PRIVATE ${PrivateSources}
    PUBLIC ${PublicSources})

target_compile_features(cg
    PUBLIC cxx_std_17)

target_include_directories(cg
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/cg>)

target_compile_definitions(cg
    PRIVATE
    LEGACY_MODE=$<BOOL:${LEGACY_MODE}>;
    $<$<NOT:$<CONFIG:Debug>>:NDEBUG>;
    MIXED_PRECISION=$<BOOL:${MIXED_PRECISION}>;)

add_subdirectory(ext/cparse EXCLUDE_FROM_ALL)
add_subdirectory(ext/boost EXCLUDE_FROM_ALL)
#add_subdirectory(ext/vcl-2.01.04 EXCLUDE_FROM_ALL)
add_subdirectory(ext/version2-2.02.00 EXCLUDE_FROM_ALL)
add_subdirectory(ext/yaml-cpp-0.7.0 EXCLUDE_FROM_ALL)

find_package(Eigen3 3.4 EXACT QUIET)
if (NOT "${Eigen3_FOUND}")
    set(MESSAGE_QUIET ON)
    set(BUILD_TESTING OFF)
    add_subdirectory(ext/eigen-3.4.0 EXCLUDE_FROM_ALL)
    unset(MESSAGE_QUIET)
endif ()

find_package(OpenMP REQUIRED)
target_link_libraries(cg PRIVATE OpenMP::OpenMP_CXX)

target_link_libraries(cg
    PUBLIC vcl yaml-cpp cparse Boost_ Eigen3::Eigen)

target_compile_options(cg
    PUBLIC
    -Wall;-Wextra;-Wpedantic;-Wno-unused-parameter;-Wno-c99-extensions;
    $<$<CONFIG:Release,Staging>:-Ofast;-march=native>;
    $<$<CONFIG:Debug,Staging>:-gdwarf-4>;)

set_target_properties(cg
    PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION TRUE)

file(GLOB_RECURSE DataFiles CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/data/**")

#add_custom_target(remove-data-dir
#    COMMAND ${CMAKE_COMMAND} -E rm -rf
#    "${CMAKE_CURRENT_BINARY_DIR}/data")

add_custom_target(copy-data-dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_BINARY_DIR}/data"
    DEPENDS ${DataFiles})

file(GLOB_RECURSE Scripts CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/scripts/**")

add_custom_target(copy-scripts-dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/scripts"
    "${CMAKE_CURRENT_BINARY_DIR}/scripts"
    DEPENDS ${Scripts})

add_custom_command(
    TARGET cg
    POST_BUILD
    COMMAND objdump -drlwC -Mintel cg >cg.s
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_dependencies(cg copy-data-dir copy-scripts-dir)



