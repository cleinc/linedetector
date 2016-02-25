cmake_minimum_required(VERSION 2.8)
enable_testing()

# make release version
set(CMAKE_BUILD_TYPE RelWithDebInfo)
#set(CMAKE_BUILD_TYPE release)
set(CMAKE_VERBOSE_MAKEFILE yes)
#set(CMAKE_CXX_FLAGS "-w -g -DNDEBUG")
set(CMAKE_CXX_FLAGS "-O3 -DNDEBUG")
set(CMAKE_BUILD_TYPE Release)

find_package(OpenCV PATHS "/usr/local/share/OpenCV")
include_directories(${OpenCV_INCLUDE_DIRS})

enable_testing()
if (APPLE)
  add_definitions(-DAPPLE)
endif (APPLE)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    "${PROJECT_SOURCE_DIR}/cmake_modules/")

# Default locations to search for on various platforms.
LIST(APPEND SEARCH_LIBS /usr/lib)
LIST(APPEND SEARCH_LIBS /usr/local/lib)
LIST(APPEND SEARCH_LIBS /opt/local/lib)

LIST(APPEND SEARCH_HEADERS /usr/include)
LIST(APPEND SEARCH_HEADERS /usr/local/include)
LIST(APPEND SEARCH_HEADERS /opt/local/include)

message("-- Check for Google Flags")
find_library(GFLAGS_LIB NAMES gflags PATHS ${SEARCH_LIBS})
find_path(GFLAGS_INCLUDE NAMES gflags/gflags.h PATHS ${SEARCH_HEADERS})

set(GOOGLE_LIBRARIES ${GFLAGS_LIB})

include_directories(${linedetection_SOURCE_DIR}/include)
include_directories(${GFLAGS_INCLUDE})
include_directories(.)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})

add_executable(linedetection linedetection.cpp  linedetector.cpp ) 
target_link_libraries(linedetection ${GOOGLE_LIBRARIES} ${OpenCV_LIBS})
