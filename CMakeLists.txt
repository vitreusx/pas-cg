cmake_minimum_required(VERSION 3.14)
project(pas-cg LANGUAGES CXX)

add_subdirectory(cg)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/docs)
    add_subdirectory(docs)
endif ()