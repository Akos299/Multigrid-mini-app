
project(multigrid VERSION 1.0
                DESCRIPTION "a mini app for mcell-centering geometric multigrid"
                LANGUAGES CXX)

#Configure cmake build
cmake_minimum_required(VERSION 3.7...3.31)

if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

#set source file directory
set(source_directory .)

#decide wether or not to build library
set(build_library true)
IF(${build_library})

    file(GLOB hearders ${source_directory}/*.hpp ${source_directory}/*.h)
    file(GLOB sources 
    ${source_directory}/multigrid_base.cpp
    ${source_directory}/mgrid_stack.cpp
    )

    #set up library
    add_library(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory})


    #install commands
    install(TARGETS multigrid 
            LIBRARY DESTINATION /usr/local/lib
            ARCHIVE DESTINATION /usr/local/lib 
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE)
    install(FILES ${headers}
            DESTINATION /usr/local/include/multigrid
            PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
            GROUP_READ GROUP_EXECUTE)
ELSE(${build_library}
    # get sources
    file(GLOB sources ${source_directory}/*.cpp)

    # set up executable
    add_executable(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory})
    )      


# add_library(ndarray PUBLIC ndArray.h ndArray.hpp)



