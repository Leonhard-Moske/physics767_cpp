## Build the project
#cd build
#cmake ..
#cmake --build .

cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME "test")

# project name
project(${PROJECT_NAME} LANGUAGES CXX)

# set the C++ standard
set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -O4")



include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/mul    
)

set(SOURCES
    main.cpp
)

# Add executable
add_executable(${PROJECT_NAME} ${SOURCES})

# Link libraries
# target_link_libraries(
# )
