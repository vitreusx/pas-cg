add_library(Boost_)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp)

target_sources(Boost_
    PUBLIC ${PublicSources}
    PRIVATE ${PrivateSources})

target_compile_options(Boost_
    PUBLIC
    -Wno-deprecated-declarations)

target_compile_features(Boost_
    PUBLIC cxx_std_17)

target_include_directories(Boost_
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)