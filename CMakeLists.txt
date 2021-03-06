cmake_minimum_required(VERSION 3.10)

PROJECT(hungarian)

########################
# BUILD SETTING

if(WIN32)
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
        set(DEFAILT_BUILD_TYPE "Release")
        message(${CMAKE_BUILD_TYPE})
        if(${CMAKE_BUILD_TYPE} MATCHES "Debug")
            message(${CMAKE_BUILD_TYPE})
            #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ")
        endif()
    endif(MSVC)
endif(WIN32)

IF(NOT MSVC)
    set(DEFAILT_BUILD_TYPE "Release")
    message(${CMAKE_BUILD_TYPE})
    if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
        message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' since none specified")
        set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING "Choose the type of build" FORCE)
        set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
    endif()
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
ADD_DEFINITIONS(-D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}/${CMAKE_BUILD_TYPE})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib/${CMAKE_SYSTEM_NAME}/${CMAKE_SYSTEM_PROCESSOR}/${CMAKE_BUILD_TYPE})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/)

include_directories(
    ${CMAKE_SOURCE_DIR}/src
)

add_subdirectory(example)