cmake_minimum_required(VERSION 3.8.2)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)

if ( DO_TRACES )
    add_definitions(-DTRACES)
    message( STATUS "Traces are activated" )
endif ( DO_TRACES )

if ( USE_CUDA )
    find_package(CUDA REQUIRED)
    enable_language(CUDA)
    set(CUDA_NVCC_FLAGS "-std=c++11 -arch=compute_60 -lineinfo")
    cuda_compile(cuda_aevol Algorithms.cu ThreefryGPU.cu)
endif ( USE_CUDA )

if ( USE_OMP )
    find_package(OpenMP REQUIRED)
    if ( OPENMP_FOUND )
        add_definitions(-DUSE_OMP)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ( USE_OMP )

add_library(micro_aevol ExpManager.cpp ExpManager.h AeTime.cpp AeTime.h DnaMutator.h Organism.cpp Organism.h RNA.h Protein.h Stats.cpp Stats.h Gaussian.h Threefry.cpp Threefry.h Timetracer.h bits.h bitset.h)

find_package(ZLIB REQUIRED)
if ( ZLIB_FOUND )
    target_include_directories(micro_aevol PUBLIC ${ZLIB_INCLUDE_DIRS})
    target_link_libraries(micro_aevol PUBLIC ${ZLIB_LIBRARIES})
endif ( ZLIB_FOUND )

if ( CUDA_FOUND )
    add_definitions(-DUSE_CUDA)
    cuda_add_executable(pdc_micro_aevol_gpu main.cpp ${cuda_aevol})
    # nvToolsExt for enhanced profiling (ad-hoc chunks)
    target_link_libraries(pdc_micro_aevol_gpu micro_aevol nvToolsExt)
else ( CUDA_FOUND )
    add_executable(pdc_micro_aevol main.cpp)
    target_link_libraries(pdc_micro_aevol micro_aevol)
endif ( CUDA_FOUND )
