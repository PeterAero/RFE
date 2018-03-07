project(RFE)
cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem thread )

find_package(GDAL)

if(GDAL_INCLUDE_DIR AND GDAL_LIBRARY)
   SET(GDAL_FOUND TRUE)
endif(GDAL_INCLUDE_DIR AND GDAL_LIBRARY)

if(GDAL_FOUND)
   MESSAGE(STATUS "Found GDAL: ${GDAL_LIBRARY}")
else(GDAL_FOUND)
   MESSAGE(GDAL_INCLUDE_DIR=${GDAL_INCLUDE_DIR})
   MESSAGE(GDAL_LIBRARY=${GDAL_LIBRARY})
   MESSAGE(FATAL_ERROR "Could not find GDAL")
endif(GDAL_FOUND)

set(HEADER_FILES "./c++/data_management/datamanager.h"
                 "./c++/image_proc/imageproc.h"
                 "./c++/parallel_processing/parallelprocessing.h"
                 "./c++/OSM/osm.h")

set(SRC_FILES "./c++/data_management/datamanager.cpp"
              "./c++/image_proc/imageproc.cpp"
              "./c++/parallel_processing/parallelprocessing.cpp"
              "./c++/OSM/osm.cpp")

aux_source_directory(. SRC_LIST)

set(INCLUDE_DIRS "./c++/data_management"
                 "./c++/image_proc"
                 "./c++/parallel_processing"
                 "./c++/OSM")


if(Boost_FOUND)
    MESSAGE (STATUS "Project Name: " ${PROJECT_NAME})
    include_directories(${Boost_INCLUDE_DIRS} ${GDAL_INCLUDE_DIR} ${INCLUDE_DIRS})
    add_executable(${PROJECT_NAME}  main.cpp ${SRC_FILES} ${HEADER_FILES})
    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${GDAL_LIBRARY})
endif()





