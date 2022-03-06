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

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/ext/eigen-3.4.0/cmake")

add_subdirectory(ext/ioxx/ioxx)
add_subdirectory(ext/nitro/nitro)
add_subdirectory(ext/cparse)
find_package(OpenMP REQUIRED)
find_package(Eigen3 3.4 REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

target_link_libraries(cg
    PRIVATE ioxx nitro OpenMP::OpenMP_CXX cparse Eigen3::Eigen Boost::program_options)

set_target_properties(cg
    PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION TRUE)

file(COPY data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

#add_custom_command(
#    TARGET cg
#    POST_BUILD
#    COMMAND /usr/bin/sh ${CMAKE_CURRENT_SOURCE_DIR}/scripts/gen_asm
#    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})