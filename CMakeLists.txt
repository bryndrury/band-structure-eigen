cmake_minimum_required(VERSION 3.5.0)
project(minimalcpp VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

if (APPLE)
    # include_directories("/opt/homebrew/opt/libomp/lib")
    # link_directories("/opt/homebrew/opt/libomp/lib")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Xpreprocessor -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lomp")

    find_package(Eigen3 REQUIRED)
    if (Eigen3_FOUND)
        include_directories(${EIGEN3_INCLUDE_DIR})
    endif()

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
include_directories("matrix_module")

# Add the OpenMP library
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(minimalcpp OpenMP::OpenMP_CXX)
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

# add_test(NAME 2_25_6 COMMAND minimalcpp 2 25 6)
add_test(NAME 2_250_6 COMMAND minimalcpp 2 250 6)

# add_test(NAME 3_25_6 COMMAND minimalcpp 3 25 6)
add_test(NAME 3_250_6 COMMAND minimalcpp 3 250 6)

# add_test(NAME 4_25_6 COMMAND minimalcpp 4 25 6)
add_test(NAME 4_250_6 COMMAND minimalcpp 4 250 6)

# add_test(NAME 5_25_6 COMMAND minimalcpp 5 25 6)
add_test(NAME 5_250_6 COMMAND minimalcpp 5 250 6)

# add_test(NAME 6_25_6 COMMAND minimalcpp 6 25 6)
add_test(NAME 6_250_6 COMMAND minimalcpp 6 250 6)

# add_test(NAME 7_25_6 COMMAND minimalcpp 7 25 6)
add_test(NAME 7_250_6 COMMAND minimalcpp 7 250 6)

# add_test(NAME 8_25_6 COMMAND minimalcpp 8 25 6)
add_test(NAME 8_250_6 COMMAND minimalcpp 8 250 6)

# add_test(NAME 9_25_6 COMMAND minimalcpp 9 25 6)
add_test(NAME 9_250_6 COMMAND minimalcpp 9 250 6)

# add_test(NAME 10_25_6 COMMAND minimalcpp 10 25 6)
add_test(NAME 10_250_6 COMMAND minimalcpp 10 250 6)