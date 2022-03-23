add_library(nitro INTERFACE)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp src/*.hpp src/*.h src/*.inl)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

target_sources(nitro
    INTERFACE ${PrivateSources} ${PublicSources})

target_compile_features(nitro
    INTERFACE cxx_std_17)

target_include_directories(nitro
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

add_subdirectory(ext/vcl-2.01.04)
target_link_libraries(nitro INTERFACE vcl)