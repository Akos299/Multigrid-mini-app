add_executable(poisson poisson.cpp main.cpp poisson.hpp)

target_link_libraries(poisson PUBLIC multigrid_library ndarray_library)
target_compile_features(poisson PUBLIC cxx_std_17)