# Root/CMakeLists.txt

# Specify the minimum version of CMake required
cmake_minimum_required(VERSION 3.10)

# Define the project name and the languages used (CXX for main code, C for argtable3)
project(DYNAMIC_SET_COVER_BENCHMARK CXX C)

# Require C++17 standard for modern features
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set optimization flags for Release builds
# Default to Release mode if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add -O3 optimization for Release builds
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3")
endif()

# Supercomputer / AMD EPYC specific optimizations
# Detect if building on supercomputer or if SUPERCOMPUTER_OPTIMIZE is set
message(STATUS "Checking for SUPERCOMPUTER_OPTIMIZE: ${CMAKE_SYSTEM_PROCESSOR}")
if(DEFINED ENV{SUPERCOMPUTER_OPTIMIZE})
    message(STATUS "SUPERCOMPUTER_OPTIMIZE is defined: $ENV{SUPERCOMPUTER_OPTIMIZE}")
endif()

if(DEFINED ENV{SUPERCOMPUTER_OPTIMIZE} OR CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64")
    # AMD EPYC 7763 (Zen3) architecture flags
    # Use -march=znver3 for EPYC 7763, or -march=native for auto-detection
    if(DEFINED ENV{SUPERCOMPUTER_OPTIMIZE})
        message(STATUS "Adding AMD EPYC optimizations: -march=znver3 -mtune=znver3")
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=znver3 -mtune=znver3" CACHE STRING "Flags used by the compiler during release builds" FORCE)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=znver3 -mtune=znver3" CACHE STRING "Flags used by the compiler during release builds" FORCE)
    else()
        # For other systems, use native if explicitly requested
        option(USE_NATIVE_ARCH "Use native CPU architecture optimizations" OFF)
        if(USE_NATIVE_ARCH)
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native" CACHE STRING "Flags used by the compiler during release builds" FORCE)
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native" CACHE STRING "Flags used by the compiler during release builds" FORCE)
        endif()
    endif()
    
    # Additional optimization flags for AMD EPYC
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        message(STATUS "Adding vectorization flags: -mavx2 -mfma")
        # Vectorization support (AVX2, FMA)
        set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx2 -mfma" CACHE STRING "Flags used by the compiler during release builds" FORCE)
        set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -mavx2 -mfma" CACHE STRING "Flags used by the compiler during release builds" FORCE)
        
        # Link-time optimization (can improve performance)
        option(ENABLE_LTO "Enable Link Time Optimization" ON)
        if(ENABLE_LTO)
            message(STATUS "Adding LTO flag: -flto")
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -flto" CACHE STRING "Flags used by the compiler during release builds" FORCE)
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -flto" CACHE STRING "Flags used by the compiler during release builds" FORCE)
            set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto" CACHE STRING "Flags used by the linker during release builds" FORCE)
        endif()

    endif()
endif()

# --- Argtable3 Library Definition ---
add_library(argtable3 STATIC
    extern/argtable3.c
)

# Compile argtable3.c as C (not C++) to avoid void* cast issues
set_source_files_properties(extern/argtable3.c PROPERTIES LANGUAGE C)

target_include_directories(argtable3
    # Add the 'extern' directory as an include path for Argtable3
    # Both PRIVATE (for argtable3.c's own includes) and PUBLIC (for dynsc.cpp's includes)
    PUBLIC
        extern
)
# --- End Argtable3 ---

# Add the 'lib' and 'app' directories to the build, CMake will look for a CMakeLists.txt file in each of these folders
add_subdirectory(lib)
add_subdirectory(app)