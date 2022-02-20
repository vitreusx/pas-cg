add_library(cparse)

target_sources(cparse
    PRIVATE
    src/builtin-features.cpp src/shunting-yard.cpp src/packToken.cpp
    src/functions.cpp src/containers.cpp)

target_compile_features(cparse
    PUBLIC cxx_std_11)

target_compile_options(cparse
    PRIVATE
    $<$<CONFIG:Debug>:-g>)

function(add_include_directories target dir)
    target_include_directories(${target}
        PRIVATE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/${dir}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include/${dir}>
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>)
endfunction()

add_include_directories(cparse cparse)

find_package(Catch2 2 REQUIRED)

add_executable(test_shunting_yard
    src/test-shunting-yard.cpp src/catch.cpp)

target_link_libraries(test_shunting_yard
    PRIVATE Catch2::Catch2 cparse)

add_include_directories(test_shunting_yard cparse)
