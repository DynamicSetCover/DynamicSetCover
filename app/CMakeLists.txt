# app/CMakeLists.txt

# 1. Define the Argtable3 source file location
set(ARGTABLE3_SOURCE ${CMAKE_SOURCE_DIR}/extern/argtable3.c)

# 2. Create the executable. Includes main app source and the Argtable3 source.
add_executable(dynsc_benchmark 
    dynsc.cpp
)

# 3. Link the executable to the required libraries.
target_link_libraries(dynsc_benchmark 
    PRIVATE 
        dynsc_lib       # The static library built in Step 2
        m                   # Standard math library (for sin/cos/pow/etc. on Linux)
        argtable3
)

# 4. Ensure the compiler can find the Argtable3 header
target_include_directories(dynsc_benchmark
    PRIVATE
        ${CMAKE_SOURCE_DIR}/extern # Directory containing argtable3.h
)
