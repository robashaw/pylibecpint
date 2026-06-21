cmake_minimum_required(VERSION 3.15)
project(pylibecpint LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(LIBECPINT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LIBECPINT_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(libecpint)

pybind11_add_module(_core bindings/bindings.cpp)
target_link_libraries(_core PRIVATE ecpint)
target_include_directories(_core PRIVATE
    ${CMAKE_SOURCE_DIR}/libecpint/include
    ${CMAKE_SOURCE_DIR}/libecpint/include/libecpint
    ${CMAKE_BINARY_DIR}/libecpint/include/libecpint
)

install(TARGETS _core DESTINATION pylibecpint)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/libecpint/share/libecpint/xml
        DESTINATION pylibecpint/data)
