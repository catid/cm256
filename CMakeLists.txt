cmake_minimum_required(VERSION 3.5)
project(cm256)

set(CMAKE_CXX_STANDARD 11)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "" FORCE)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

option(MARCH_NATIVE "Use -march=native option" ON)

set(LIB_SOURCE_FILES
    src/cm256.cpp
    include/cm256.h
    src/gf256.cpp
    include/gf256.h
)

set(UNIT_TEST_SOURCE_FILES
    unit_test/main.cpp
    unit_test/SiameseTools.cpp
    unit_test/SiameseTools.h
)

if(MSVC)
else()
    set(CMAKE_CXX_FLAGS "-Wall -Wextra")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

if(MARCH_NATIVE)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -march=native")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native")
endif()

add_library(cm256 ${LIB_SOURCE_FILES})
set_target_properties(cm256 PROPERTIES VERSION 2)
set_target_properties(cm256 PROPERTIES SOVERSION 2)
target_include_directories(cm256 PUBLIC ${PROJECT_SOURCE_DIR}/include)

add_executable(unit_test ${UNIT_TEST_SOURCE_FILES})
target_link_libraries(unit_test cm256)

install(TARGETS cm256
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY include DESTINATION ${CMAKE_INSTALL_PREFIX})
