cmake_minimum_required(VERSION 3.16)
project(pyopatra)

option(BUILD_CPP "Build the C++ Library" OFF)
option(BUILD_PYTHON "Build the Python Library" OFF)
option(BUILD_TEST "Build the test suite" OFF)
option(CODE_COVERAGE "Enable code coverage report" OFF)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(MPI REQUIRED COMPONENTS CXX)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Eigen3 REQUIRED NO_MODULE)

find_package(Python REQUIRED COMPONENTS Interpreter Development)
set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
find_package(pybind11 REQUIRED)
set(PYBIND11_CPP_STANDARD -std=c++14)

if (${BUILD_TEST})
    find_package(Catch2 2 REQUIRED)
endif ()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

if (${BUILD_PYTHON})
    message(STATUS "Building Python library...")
    pybind11_add_module(pyopatra_pybind
            src/PYOPATRA/cpp/bindings.cpp
            src/PYOPATRA/cpp/particle.h
            src/PYOPATRA/cpp/particle.cpp
            src/PYOPATRA/cpp/mesh/mesh_element.h
            src/PYOPATRA/cpp/coordinate.h
            src/PYOPATRA/cpp/illnode.h
            src/PYOPATRA/cpp/list.h
            src/PYOPATRA/cpp/mesh/mesh.cpp
            src/PYOPATRA/cpp/mesh/mesh.h
            src/PYOPATRA/cpp/mesh/mesh_vertex.inl
            src/PYOPATRA/cpp/mesh/mesh_vertex.h
            src/PYOPATRA/cpp/mesh/mesh_water_column.inl
            src/PYOPATRA/cpp/mesh/mesh_water_column.h
            src/PYOPATRA/cpp/mesh/mesh_vertex.cpp
            src/PYOPATRA/cpp/mesh/mesh_element.cpp
            src/PYOPATRA/cpp/particle_list.h
            src/PYOPATRA/cpp/inversion_tools/objective_functions.h
            src/PYOPATRA/cpp/util.h
            src/PYOPATRA/cpp/solver.h
            )
    target_include_directories(pyopatra_pybind PUBLIC SYSTEM ${MPI_CXX_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(pyopatra_pybind PUBLIC ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen)
    target_compile_definitions(pyopatra_pybind PUBLIC ${HDF5_DEFINITIONS})
    target_compile_options(pyopatra_pybind PUBLIC -O3 -Wall -Wextra -pedantic -Werror)
endif()

if (${BUILD_CPP} OR ${BUILD_TEST})
    add_library(pyopatra_lib
            src/PYOPATRA/cpp/particle.h
            src/PYOPATRA/cpp/particle.cpp
            src/PYOPATRA/cpp/mesh/mesh_element.h
            src/PYOPATRA/cpp/coordinate.h
            src/PYOPATRA/cpp/illnode.h
            src/PYOPATRA/cpp/list.h
            src/PYOPATRA/cpp/mesh/mesh.cpp
            src/PYOPATRA/cpp/mesh/mesh.h
            src/PYOPATRA/cpp/mesh/mesh_vertex.inl
            src/PYOPATRA/cpp/mesh/mesh_vertex.h
            src/PYOPATRA/cpp/mesh/mesh_water_column.inl
            src/PYOPATRA/cpp/mesh/mesh_water_column.h
            src/PYOPATRA/cpp/mesh/mesh_vertex.cpp
            src/PYOPATRA/cpp/mesh/mesh_element.cpp
            src/PYOPATRA/cpp/particle_list.h
            src/PYOPATRA/cpp/inversion_tools/objective_functions.h
            src/PYOPATRA/cpp/util.h
            src/PYOPATRA/cpp/solver.h)
    target_include_directories(pyopatra_lib PUBLIC SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(pyopatra_lib PUBLIC ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen)
    target_compile_definitions(pyopatra_lib PUBLIC ${HDF5_DEFINITIONS})
    target_compile_options(pyopatra_lib PUBLIC -g -Wall -Wextra -pedantic -Werror)
endif ()

if (${BUILD_TEST})
    if (${CODE_COVERAGE})
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
        set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    endif ()
    add_executable(pyopatra_test
            test/cpp_tests/test_base.cpp
            test/cpp_tests/test_particle.cpp
            test/cpp_tests/test_mesh_vertex.cpp
            test/cpp_tests/test_mesh.cpp
            test/cpp_tests/test_ill.cpp
            test/cpp_tests/test_mesh_element.cpp
            test/cpp_tests/test_coordinate.cpp
            test/cpp_tests/test_obj_function.cpp)
    target_link_libraries(pyopatra_test pyopatra_lib Catch2::Catch2)
    target_include_directories(pyopatra_test PUBLIC src )
endif ()