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
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/pas-cg>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/ioxx>)

target_compile_options(pas-cg
    PRIVATE
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:-Ofast;-march=native>>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

file(COPY data DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")