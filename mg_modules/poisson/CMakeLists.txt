add_library(poisson_library poisson.cpp)
target_include_directories(poisson_library PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ../include)

add_executable(poisson main.cpp)
target_link_libraries(poisson PUBLIC multigrid_library poisson_library)
target_compile_features(poisson PUBLIC cxx_std_17)
target_include_directories(poisson PUBLIC ../include )
set_target_properties(poisson
    PROPERTIES COMPILER_FLAGS "-g -m64 -arch x86_64 -msse -Wall -pedantic -v"
               LINKER_FLAGS "-g -m64 -arch x86_64 -msse -v")