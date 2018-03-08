project(RFE)
cmake_minimum_required(VERSION 3.4)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)


set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem thread )



find_library(MYGDAL NAMES gdal PATHS "/home/fisc_p0/Desktop/local_libs/gdal/lib" NO_DEFAULT_PATH)

if(MYGDAL)
    set(GDAL_LIBRARY ${MYGDAL})
    MESSAGE(STATUS "Found GDAL LIB in non-system eco-system.")
    MESSAGE(STATUS "Used GDAL library: " ${GDAL_LIBRARY})

    set(MY_GDAL_INCLUDES "/home/fisc_p0/Desktop/local_libs/gdal/.libs"
    "/home/fisc_p0/Desktop/local_libs/gdal/gcore"
    "/home/fisc_p0/Desktop/local_libs/gdal/port"
    "/home/fisc_p0/Desktop/local_libs/gdal/ogr"
    "/home/fisc_p0/Desktop/local_libs/gdal/ogr/ogrsf_frmts"
    "/home/fisc_p0/Desktop/local_libs/gdal/data"
    "/home/fisc_p0/Desktop/local_libs/gdal/alg" )
    MESSAGE(STATUS "Used GDAL library: " ${MY_GDAL_INCLUDES})
    set(GDAL_INCLUDE_DIR ${MY_GDAL_INCLUDES})

else(MYGDAL)

    MESSAGE(STATUS "No Lib Found: " ${MYGDAL})
    find_package(GDAL REQUIRED )

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

endif(MYGDAL)







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
#    include_directories(${Boost_INCLUDE_DIRS} ${INCLUDE_DIRS} "/home/fisc_p0/Desktop/local_libs/gdal/.libs"
# "/home/fisc_p0/Desktop/local_libs/gdal/gcore"
# "/home/fisc_p0/Desktop/local_libs/gdal/port"
# "/home/fisc_p0/Desktop/local_libs/gdal/ogr"
# "/home/fisc_p0/Desktop/local_libs/gdal/ogr/ogrsf_frmts"
# "/home/fisc_p0/Desktop/local_libs/gdal/data"
# "/home/fisc_p0/Desktop/local_libs/gdal/alg")

    add_executable(${PROJECT_NAME}  main.cpp ${SRC_FILES} ${HEADER_FILES})

    target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} ${GDAL_LIBRARY})
    #target_link_libraries(${PROJECT_NAME} ${Boost_LIBRARIES} "/home/fisc_p0/Desktop/local_libs/gdal/lib/libgdal.so")

endif()





