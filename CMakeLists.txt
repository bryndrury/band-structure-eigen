cmake_minimum_required(VERSION 3.5.0)
project(minimalcpp VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Xpreprocessor -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /openmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(minimalcpp main.cpp)

# Include Eigen headers
include_directories("${CMAKE_SOURCE_DIR}/dependencies/eigen")

# Link the custom library
include_directories("module")

# Add the OpenMP library
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(minimalcpp OpenMP::OpenMP_CXX)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)