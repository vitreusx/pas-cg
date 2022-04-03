add_library(Boost INTERFACE)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

target_sources(Boost
    INTERFACE ${PublicSources})

target_compile_features(Boost
    INTERFACE cxx_std_17)

target_include_directories(Boost
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)