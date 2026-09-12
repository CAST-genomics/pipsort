cmake_minimum_required(VERSION 3.16)

include(ExternalProject)

project(PIPSORT 
    VERSION 1.0.1 
    LANGUAGES CXX
)
configure_file(PIPSORTConfig.h.in PIPSORTConfig.h)


list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

########## dependencies #########

ExternalProject_Add(pgsl
        PREFIX ${CMAKE_BINARY_DIR}/thirdparty/gsl
        GIT_REPOSITORY "https://github.com/ampl/gsl.git"
        GIT_TAG "v2.5.0"
        UPDATE_COMMAND ""
        BUILD_IN_SOURCE 1
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/thirdparty/gsl -DBUILD_SHARED_LIBS=OFF
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        LOG_DOWNLOAD 1
        )

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_BINARY_DIR}/thirdparty/gsl/include)

find_package(OpenMP REQUIRED)
find_package(BLAS REQUIRED)

########## PIPSORT code #########

# add the executable
add_executable(${PROJECT_NAME} pipsort.cpp postcal.cpp sss_postcal.cpp util.cpp)

target_include_directories(${PROJECT_NAME} PUBLIC "${PROJECT_BINARY_DIR}")
target_include_directories(${PROJECT_NAME} PUBLIC 
    ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/armadillo/include
)

target_link_libraries(${PROJECT_NAME} PRIVATE gsl)
target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(${PROJECT_NAME} PRIVATE ${BLAS_LIBRARIES})
target_link_directories(${PROJECT_NAME} PUBLIC ${CMAKE_BINARY_DIR}/thirdparty/gsl/lib)

install(TARGETS ${PROJECT_NAME} DESTINATION bin)

add_dependencies(${PROJECT_NAME} pgsl)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "") #required

set(CPACK_PACKAGE_VERSION "${PIPSORT_VERSION_MAJOR}.${PIPSORT_VERSION_MINOR}.${PIPSORT_VERSION_PATCH}")
set(CPACK_PACKAGE_VERSION_MAJOR "${PIPSORT_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${PIPSORT_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${PIPSORT_VERSION_PATCH}")

set(CPACK_RESOURCE_FILE_README "${CMAKE_SOURCE_DIR}/README.md")

include(CPack)