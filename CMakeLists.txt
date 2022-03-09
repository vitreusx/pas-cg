cmake_minimum_required(VERSION 3.14)
project(pas-cg LANGUAGES CXX)

if (NOT DEFINED CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

add_subdirectory(cg)
add_subdirectory(docs EXCLUDE_FROM_ALL)