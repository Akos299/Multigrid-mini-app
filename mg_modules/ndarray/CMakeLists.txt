
add_library(ndarray_library ndArray.cpp ndArray.hpp)
target_compile_features(ndarray_library PUBLIC cxx_std_17)