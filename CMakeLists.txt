cmake_minimum_required(VERSION 3.19)
project(Heisenberg_Chain_1D_C_)

set(CMAKE_CXX_STANDARD 20)

if (WIN32)
    include_directories("C:\\Users\\masch\\Code\\eigen-3.4.0")
endif (WIN32)
if (UNIX)
    include_directories("/home/mmaschke/BA Code/Libraries/eigen-3.4.0")
endif (UNIX)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wa,-mbig-obj")

add_executable(Heisenberg1D main.cpp main.h helperFunctions.cpp helperFunctions.h constants.h)
