add_executable(pas-cg)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp src/*.hpp src/*.h src/*.inl)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

target_sources(pas-cg
    PRIVATE ${PrivateSources}
    PUBLIC ${PublicSources})

target_compile_features(pas-cg
    PUBLIC cxx_std_17)

target_include_directories(pas-cg
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/cg>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/cg>)

target_compile_options(pas-cg
    PRIVATE
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:-Ofast;-march=native>>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

add_subdirectory(ext/ioxx/ioxx)
add_subdirectory(ext/nitro/nitro)
add_subdirectory(ext/cparse)
find_package(OpenMP REQUIRED)

target_link_libraries(pas-cg PRIVATE ioxx nitro OpenMP::OpenMP_CXX cparse)

file(COPY data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")