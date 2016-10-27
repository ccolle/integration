
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.10)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_SOURCE_DIR}/cmake/Modules/")
MESSAGE("cmake_module_path is ${CMAKE_MODULE_PATH}")

FIND_PACKAGE(GSL REQUIRED)

set(CMAKE_CXX_FLAGS "-std=c++11 -Wall -Wextra -Wshadow -lm")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

MESSAGE("-- [Info] CMAKE_CXX_FLAGS  are ${CMAKE_CXX_FLAGS}")

include_directories("${PROJECT_SOURCE_DIR}/inc" 
                    "${PROJECT_SOURCE_DIR}/Cuba"
                    "${PROJECT_SOURCE_DIR}/cubature")


add_subdirectory("${PROJECT_SOURCE_DIR}/src")
add_subdirectory("${PROJECT_SOURCE_DIR}/cubature")

set(Cuba_dir "${PROJECT_SOURCE_DIR}/Cuba")

#execute_process(COMMAND readlink -f Cuba
#                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                OUTPUT_VARIABLE Cuba_dir
#                OUTPUT_STRIP_TRAILING_WHITESPACE
#                )

include(ExternalProject)
ExternalProject_Add(
                Cuba
                SOURCE_DIR         "${Cuba_dir}"
                CONFIGURE_COMMAND  ${Cuba_dir}/configure --prefix=${Cuba_dir} --libdir=${PROJECT_SOURCE_DIR}/lib
                BINARY_DIR         "${Cuba_dir}"
                BUILD_COMMAND      make
                INSTALL_COMMAND    make install
                )
