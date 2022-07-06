cmake_minimum_required(VERSION 3.21)
project(SpinChainED)

set(CMAKE_CXX_STANDARD 20)

find_package(MKL REQUIRED)
find_package(OpenMP REQUIRED)

if (WIN32)
    include_directories("D:/Code/C++/eigen-3.4.0")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wa,-mbig-obj")
    set(MKLROOT "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -I"%MKLROOT%"\\include -DEIGEN_USE_MKL_ALL")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -I"%MKLROOT%"\\include -DEIGEN_USE_MKL_ALL")
endif (WIN32)
if (UNIX)
    include_directories("/home/mmaschke/BA_Code/Libraries/eigen-3.4.0")
    set(MKLROOT /opt/intel/oneapi/mkl/latest)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -pg -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -m64 -I${MKLROOT}/include -DEIGEN_USE_MKL_ALL")
endif (UNIX)

add_executable(SpinChain)
add_subdirectory(src)

if (WIN32)
    target_link_libraries(SpinChain PUBLIC OpenMP::OpenMP_CXX -Wl,--start-group
            ${MKLROOT}/lib/intel64/mkl_intel_lp64.lib ${MKLROOT}/lib/intel64/mkl_intel_thread.lib
            ${MKLROOT}/lib/intel64/mkl_core.lib -Wl,--end-group -lgomp -lpthread -lm)
endif(WIN32)
if (UNIX)
    target_link_libraries(SpinChain PUBLIC OpenMP::OpenMP_CXX -Wl,--start-group
            ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a
            ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl)
endif(UNIX)