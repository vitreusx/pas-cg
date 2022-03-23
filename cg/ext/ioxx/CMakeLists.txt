add_library(ioxx)

file(GLOB_RECURSE PrivateSources CONFIGURE_DEPENDS
    src/*.cpp src/*.hpp src/*.h src/*.inl)

file(GLOB_RECURSE PublicSources CONFIGURE_DEPENDS
    include/*.hpp include/*.h include/*.inl)

target_sources(ioxx
    PRIVATE ${PrivateSources}
    PUBLIC ${PublicSources})

target_compile_features(ioxx
    PUBLIC cxx_std_17)

target_include_directories(ioxx
    PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/ioxx>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
    $<INSTALL_INTERFACE:include/ioxx>
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

target_compile_options(ioxx
    PRIVATE
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:Release>:$<$<COMPILE_LANGUAGE:CXX>:-Ofast;-march=native>>
    $<$<COMPILE_LANGUAGE:CXX>:-Wall;-Wextra;-Wpedantic>)

find_package(yaml-cpp QUIET)
if (NOT ${yaml-cpp_FOUND})
    add_subdirectory(ext/yaml-cpp-0.7.0)
endif()

target_link_libraries(ioxx PUBLIC yaml-cpp)