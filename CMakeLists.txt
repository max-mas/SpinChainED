cmake_minimum_required(VERSION 3.21)
project(SpinChainED)

set(CMAKE_CXX_STANDARD 20)

if (WIN32)
    include_directories("D:\\Code\\C++\\eigen-3.4.0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wa,-mbig-obj")
endif (WIN32)
if (UNIX)
    include_directories("/home/mmaschke/BA_Code/Libraries/eigen-3.4.0")
endif (UNIX)

find_package(MKL REQUIRED)
find_package(OpenMP REQUIRED)

set(MKLROOT /opt/intel/oneapi/mkl/latest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL")

add_executable(SpinChainED main.cpp main.h helperFunctions.cpp helperFunctions.h constants.h hamiltonianBuilders.cpp hamiltonianBuilders.h thermodynamics.cpp thermodynamics.h dataGenerators.cpp dataGenerators.h)

target_link_libraries(SpinChainED PUBLIC OpenMP::OpenMP_CXX -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl)
