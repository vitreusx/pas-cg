cmake_minimum_required(VERSION 3.14)
project(pas-cg LANGUAGES CXX)

add_subdirectory(cg)

if (EXISTS docs)
    add_subdirectory(docs EXCLUDE_FROM_ALL)
endif ()