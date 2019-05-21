CMAKE_MINIMUM_REQUIRED(VERSION 3.0.0)

PROJECT(firesheds)

SET(CMAKE_CXX_STANDARD 11)

FIND_PACKAGE(GDAL REQUIRED)

INCLUDE(FindOpenMP)
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    IF(MSVC)
        SET(CMAKE_CXX_FLAGS "${cmake_cxx_flags} ${OpenMP_CXX_FLAGS} /EHsc")
    ELSE(MSVC) # Not using MSVC
        SET(CMAKE_CXX_FLAGS "${cmake_cxx_flags} ${OpenMP_CXX_FLAGS}")
    ENDIF(MSVC)
ELSE(OPENMP_FOUND)
    MESSAGE("OpenMP could not be found on your system, disabling support")
ENDIF(OPENMP_FOUND)

INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src
					${GDAL_INCLUDE_DIRS})

ADD_EXECUTABLE(${PROJECT_NAME} "src/main.cpp"
				"src/sqlite3.c"
				"src/WfipsGrid.cpp"
				"src/WfipsGridData.cpp")

TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${GDAL_LIBRARY})
