# project name
project(multigrid)      

#Configure cmake build
cmake_minimum_required(VERSION 3.7...3.31)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

#set C++ version to be used
set(CMAKE_CXX_STANDARD 17) 
# set PATH to cmake config files       
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "cmake") 
# 
set(CMAKE_VERBOSE_MAKEFILE false)


# Let's nicely support folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

add_subdirectory(poisson)

# Here are sources code for 
add_subdirectory(src)

# the  ndArray lib is here
add_subdirectory(ndarray)



