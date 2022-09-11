option(USE_OMP "Use OpenMP" ON)
option(SINGLE_FILE "Concatenate all source files" OFF)
option(USE_MIXED_PRECISION "Use floats for computing the forces and doubles for integration, instead of using doubles throughout" OFF)
option(USE_VECTORIZED_IMPLS "Use vectorized implementations of the algorithms" ON)

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

if (SINGLE_FILE)
    file(GLOB_RECURSE SourceFiles CONFIGURE_DEPENDS
        src/*.cpp)

    add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/.fused.cpp
        COMMAND awk '{print}' ${SourceFiles} >.fused.cpp
        DEPENDS ${SourceFiles}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_custom_target(fused_sources
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/.fused.cpp)

    add_dependencies(cg fused_sources)

    target_sources(cg
        PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/.fused.cpp)
else ()
    file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
        src/*.cpp src/*.hpp src/*.h src/*.inl)

    file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
        include/*.hpp include/*.h include/*.inl)

    target_sources(cg
        PRIVATE ${PrivateSources}
        PUBLIC ${PublicSources})
endif ()

target_compile_features(cg
    PUBLIC cxx_std_17)

target_include_directories(cg
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/cg>)

target_compile_options(cg
    PRIVATE
    -Wall;-Wextra;-Wpedantic;-Wno-unused-but-set-parameter;-Wno-c99-extensions;
    $<$<CONFIG:Deployment>:-Ofast>;
    $<$<CONFIG:Release>:-Ofast;-march=native>;
    $<$<CONFIG:Staging>:-Ofast;-march=native;-g>;
    $<$<CONFIG:Debug>:-g>;
    $<$<CONFIG:Compat>:-g>)

target_compile_definitions(cg
    PRIVATE
    COMPAT_MODE
    $<$<NOT:$<CONFIG:Debug>>:NDEBUG>
    USE_MIXED_PRECISION=$<BOOL:${USE_MIXED_PRECISION}>
    USE_VECTORIZED_IMPLS=$<BOOL:${USE_VECTORIZED_IMPLS}>)

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

if (${USE_OMP})
    find_package(OpenMP REQUIRED)
    target_link_libraries(cg PRIVATE OpenMP::OpenMP_CXX)
endif ()

target_link_libraries(cg
    PRIVATE vcl yaml-cpp cparse Boost Eigen3::Eigen)

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

add_custom_command(
    TARGET cg
    POST_BUILD
    COMMAND objdump -drlwC -Mintel cg >cg.s
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

add_dependencies(cg copy-data-dir)
