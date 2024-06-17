cmake_minimum_required(VERSION 3.5.0)
project(minimalcpp VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Xpreprocessor -fopenmp -lomp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(minimalcpp main.cpp)

# link the Eigen library
include_directories("dependencies/eigen")
find_package(Eigen3 REQUIRED)
if(Eigen3_FOUND)
    target_link_libraries(minimalcpp Eigen3::Eigen)
endif()

# Add the gsl library
# include_directories(/opt/homebrew/Cellar/gsl/2.7/include)
# link_directories(/opt/homebrew/Cellar/gsl/2.7/lib)
# find_package(GSL REQUIRED)
# if(GSL_FOUND)
#     target_link_libraries(minimalcpp GSL::gsl GSL::gslcblas)
# endif()

# Add the OpenMP library
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(minimalcpp OpenMP::OpenMP_CXX)
endif()


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
