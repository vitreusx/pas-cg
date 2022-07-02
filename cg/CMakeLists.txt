option(USE_OMP "Use OpenMP" ON)
option(SINGLE_FILE "Concatenate all source files" OFF)

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
    -Wall;-Wextra;-Wpedantic;-Wno-unused-but-set-parameter;
    $<$<CONFIG:Deployment>:-Ofast>;
    $<$<CONFIG:Release>:-Ofast;-march=native>;
    $<$<CONFIG:Staging>:-Ofast;-march=native;-g>;
    $<$<CONFIG:Debug>:-g>;
    $<$<CONFIG:Compat>:-g>)

target_compile_definitions(cg
    PRIVATE
    COMPAT_MODE
    $<$<NOT:$<CONFIG:Debug>>:NDEBUG>)

add_subdirectory(ext/cparse EXCLUDE_FROM_ALL)
add_subdirectory(ext/boost EXCLUDE_FROM_ALL)
add_subdirectory(ext/vcl-2.01.04 EXCLUDE_FROM_ALL)
add_subdirectory(ext/yaml-cpp-0.7.0 EXCLUDE_FROM_ALL)

find_package(Eigen3 3.4 QUIET)
if (NOT "${Eigen3_FOUND}")
    add_subdirectory(ext/eigen-3.4.0 EXCLUDE_FROM_ALL)
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

add_custom_target(copy_data_dir
    COMMAND ${CMAKE_COMMAND} -E copy_directory
    "${CMAKE_CURRENT_SOURCE_DIR}/data"
    "${CMAKE_CURRENT_BINARY_DIR}/data")

add_dependencies(cg copy_data_dir)
