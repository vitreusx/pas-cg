cmake_minimum_required(VERSION 3.14)
project(pas-cg LANGUAGES CXX)
enable_testing()

add_subdirectory(cg)

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/parity)
    add_subdirectory(parity)
endif ()

if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/docs)
    add_subdirectory(docs)
endif ()