add_executable(cg)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp src/*.hpp src/*.h src/*.inl)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

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
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/cg>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/cg>)

target_compile_options(cg
    PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:
    -Wall;
    $<$<CONFIG:Release>:-Ofast;-march=native>;
    $<$<CONFIG:Debug>:-g>>)

add_subdirectory(ext/ioxx/ioxx EXCLUDE_FROM_ALL)
add_subdirectory(ext/nitro/nitro EXCLUDE_FROM_ALL)
add_subdirectory(ext/cparse EXCLUDE_FROM_ALL)

find_package(OpenMP REQUIRED)

find_package(Eigen3 3.4 QUIET)
if (NOT ${Eigen3_FOUND})
    add_subdirectory(ext/eigen-3.4.0 EXCLUDE_FROM_ALL)
endif ()

target_link_libraries(cg
    PRIVATE ioxx nitro OpenMP::OpenMP_CXX cparse Eigen3::Eigen)

set_target_properties(cg
    PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION TRUE)

file(COPY data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")