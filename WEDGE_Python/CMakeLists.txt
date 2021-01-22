cmake_minimum_required(VERSION  3.16.3)
project(SNMF_inner)
add_compile_options(-std=c++11)
set(pybind11_DIR ${PROJECT_SOURCE_DIR}/env/lib/python3.8/site-packages/pybind11/share/cmake/pybind11) #change by the path of pybind11

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g")
set(CMAKE_BUILD_TYPE "Release")
set(PYBIND11_CPP_STANDARD -std=c++11)

include_directories("./")
include_directories("./env/lib/python3.8/site-packages/pybind11/include") #change by the path of pybind11
include_directories("./eigen-3.3.8")

find_package(pybind11 REQUIRED)
 
add_library(NNLS MODULE 
            NNLS.cpp
    )
target_link_libraries(NNLS
            pybind11::module
  )

# message("${PROJECT_SOURCE_DIR}")
#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/../bin)       #set the out location of executable file

SET(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})
