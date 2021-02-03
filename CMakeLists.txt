cmake_minimum_required(VERSION 3.16)
project(lptcpp)

if(NOT EXISTS ${BUILD_CPP})
    set(BUILD_CPP FALSE)
endif ()

if(NOT EXISTS ${BUILD_PYTHON})
    set(BUILD_PYTHON FALSE)
endif ()

if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(DEBUG TRUE)
endif ()

find_package(MPI REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS CXX)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (${BUILD_PYTHON})
    find_package(Python REQUIRED COMPONENTS Interpreter Development)
    set(PYTHON_EXECUTABLE ${Python_EXECUTABLE})
    find_package(pybind11 REQUIRED)
    set(PYBIND11_CPP_STANDARD -std=c++14)
endif ()

if (${DEBUG})
    find_package(Catch2 2 REQUIRED)
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -W -Wall -Wextra -pedantic -march=native")
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

if (${BUILD_PYTHON})
    pybind11_add_module(pylpt_cpp src/pyLPT/cpp/library.cpp src/pyLPT/cpp/library.h)
    target_include_directories(pylpt_cpp PUBLIC SYSTEM ${MPI_CXX_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(pylpt_cpp PUBLIC ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen)
    target_compile_definitions(pylpt_cpp PUBLIC ${HDF5_DEFINITIONS})
endif()

if (${BUILD_CPP} OR ${DEBUG})
    add_library(lptcpp_lib src/pyLPT/cpp/library.cpp src/pyLPT/cpp/library.h)
    target_include_directories(lptcpp_lib PUBLIC SYSTEM ${MPI_CXX_INCLUDE_PATH} ${HDF5_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(lptcpp_lib PUBLIC ${MPI_CXX_LIBRARIES} ${HDF5_LIBRARIES} Eigen3::Eigen)
    target_compile_definitions(lptcpp_lib PUBLIC ${HDF5_DEFINITIONS})
endif ()

if (${DEBUG})
    add_executable(lptcpp_test test/test_base.cpp)
endif ()