
add_library(ndarray_library ndArray.cpp ndArray.hpp)
target_compile_features(ndarray_library PUBLIC cxx_std_17)
set_target_properties(ndarray_library
    PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic"
               LINKER_FLAGS "-g -m64 -arch x86_64 -msse")