project(multigrid)      

# Configure cmake build
#Configure cmake build
cmake_minimum_required(VERSION 3.7...3.31)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
    cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()
set(CMAKE_CXX_STANDARD 17)        
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "cmake") 
set(CMAKE_VERBOSE_MAKEFILE false)

# Set source file directory  
set(source_directory .)   

# Decide what to build
set(build_library false) 
IF(${build_library}) 
    # Get headers & sources but not main.cpp  
    file(GLOB headers ${source_directory}/*.hpp ${source_directory}/*.h )
    file(GLOB sources 
        ${source_directory}/mgrid_base.cpp
        ${source_directory}/mgrid_stack.cpp
)    
    
    # Set up library    
    add_library(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory})

    # Install commands
    install(TARGETS multigrid
        LIBRARY DESTINATION /usr/local/lib
        ARCHIVE DESTINATION /usr/local/lib
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE)    
    install(FILES ${headers}
        DESTINATION /usr/local/include/multigrid
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE
                    GROUP_READ GROUP_EXECUTE)
ELSE(${build_library}) 
    # Get sources
    file(GLOB sources ${source_directory}/*.cpp)
    
    # Set up executable 
    add_executable(${PROJECT_NAME} ${sources})
    include_directories(${INCLUDES} ${source_directory} 
       ) 
ENDIF(${build_library})          