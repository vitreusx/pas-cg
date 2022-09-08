add_library(vcl INTERFACE)

file(GLOB_RECURSE VCL_HEADERS include/vcl/*.h)
target_sources(vcl INTERFACE "${VCL_HEADERS}")

target_include_directories(vcl
    INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

