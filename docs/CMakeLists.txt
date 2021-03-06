cmake_minimum_required(VERSION 3.14)
project(cg-docs)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(Doxygen)
find_package(Sphinx)

if ("${Doxygen_FOUND}" AND "${Sphinx_FOUND}")
    set(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/cg)
    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)

    set(CG_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../cg)
    file(GLOB_RECURSE CG_SOURCES CONFIGURE_DEPENDS
        ${CG_DIR}/include/* ${CG_DIR}/src/*)

    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
    add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
        DEPENDS ${CG_SOURCES}
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
        MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
        COMMENT "Generating Doxygen docs")

    add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/docs)
    set(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

    file(GLOB_RECURSE SPHINX_FILES CONFIGURE_DEPENDS ${SPHINX_SOURCE}/*)

    add_custom_command(OUTPUT ${SPHINX_INDEX_FILE}
        COMMAND ${SPHINX_EXECUTABLE} -b html -Dbreathe_projects.cg=${DOXYGEN_OUTPUT_DIR}/xml ${SPHINX_SOURCE} ${SPHINX_BUILD}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS ${SPHINX_FILES} ${DOXYGEN_INDEX_FILE}
        MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
        COMMENT "Generating Sphinx docs")

    add_custom_target(docs ALL DEPENDS ${SPHINX_INDEX_FILE})
endif ()