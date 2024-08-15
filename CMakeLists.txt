cmake_minimum_required(VERSION 3.5.0)
project(minimalcpp VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Xpreprocessor -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include")

if (APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Xpreprocessor -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include")
elseif(MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 /fopenmp")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(minimalcpp main.cpp)

# link the Eigen library
include_directories("${CMAKE_SOURCE_DIR}/dependencies/eigen")

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/dependencies")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/dependencies")

find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    target_link_libraries(minimalcpp Eigen3::Eigen)
endif()

# link the custom library
include_directories("module")

# Add the OpenMP library
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(minimalcpp OpenMP::OpenMP_CXX)
endif()


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
