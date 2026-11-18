cmake_minimum_required(VERSION 3.15)
project(pylibecpint LANGUAGES CXX)

find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

set(LIBECPINT_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(LIBECPINT_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
add_subdirectory(libecpint)
get_directory_property(LIBECPINT_VERSION DIRECTORY libecpint DEFINITION libecpint_VERSION)

pybind11_add_module(_core bindings/bindings.cpp)
# ecpint's include directories and C++17 requirement are propagated publicly
target_link_libraries(_core PRIVATE ecpint)
target_compile_definitions(_core PRIVATE LIBECPINT_VERSION_STR="${LIBECPINT_VERSION}")

# Only the "python" component is installed into the wheel (see install.components in
# pyproject.toml), so libecpint's own install rules (headers, static library, share/) are skipped
install(TARGETS _core DESTINATION pylibecpint COMPONENT python)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/libecpint/share/libecpint/raw
        DESTINATION pylibecpint/data COMPONENT python)
