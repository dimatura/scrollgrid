cmake_minimum_required(VERSION 2.8.3)
project(scrollgrid)

# Load catkin and all dependencies required for this package
find_package(catkin REQUIRED COMPONENTS roscpp pcl_util geom_cast cmake_modules)
include_directories(${catkin_INCLUDE_DIRS})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

find_package(Eigen REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#common commands for building c++ executables and libraries
#add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)

#target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES})
#add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})

# catkin_package parameters: http://ros.org/doc/groovy/api/catkin/html/dev_guide/generated_cmake_api.html#catkin-package
# TODO: fill in what other packages will need to use this package
catkin_package(
    DEPENDS Eigen
    CATKIN_DEPENDS roscpp pcl_util geom_cast
    INCLUDE_DIRS include
    LIBRARIES 
)

set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")

#if (CATKIN_ENABLE_TESTING)
catkin_add_gtest(${PROJECT_NAME}_scrollgrid3_tests tests/scrollgrid3.cpp)
target_link_libraries(${PROJECT_NAME}_scrollgrid3_tests ${catkin_LIBRARIES})
#endif()

include_directories(include)
add_executable(${PROJECT_NAME}_examples src/examples.cpp)
target_link_libraries(${PROJECT_NAME}_examples ${catkin_LIBRARIES})
