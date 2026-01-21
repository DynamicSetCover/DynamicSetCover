# lib/CMakeLists.txt

# Gather all .cpp files recursively from the current directory (lib/) 
# and all its subdirectories (algorithms/, io/)
file(GLOB_RECURSE LIB_SOURCES "*.cpp")

# Create a static library target named 'set_cover_lib'
add_library(dynsc_lib STATIC ${LIB_SOURCES})

# Specify where the compiler should look for headers.
# We must include the subdirectories explicitly so that files inside the 'lib/algorithms' folder can include files from 'lib/io'.
target_include_directories(dynsc_lib PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}             # For files directly in lib/
    ${CMAKE_CURRENT_SOURCE_DIR}/algorithms  # For headers in lib/algorithms
    ${CMAKE_CURRENT_SOURCE_DIR}/io          # For headers in lib/io
)