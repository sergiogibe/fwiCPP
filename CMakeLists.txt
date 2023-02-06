cmake_minimum_required(VERSION 3.16)
project(fwiCPP)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-O2")
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Armadillo REQUIRED)

set(SOURCE_FILES
    main.cpp
    Device.cpp
    Device.h
    MaterialModel.cpp
    MaterialModel.h
    Mesh.cpp
    Mesh.h
    Problem.cpp
    Problem.h
)

source_group(
    TREE "${CMAKE_CURRENT_SOURCE_DIR}"
    PREFIX "Source Files"
    FILES ${SOURCE_FILES}
)

add_executable(
    fwiCPP
    ${SOURCE_FILES}
)

target_link_libraries(
    fwiCPP
    ${ARMADILLO_LIBRARIES}
)

target_include_directories(
    fwiCPP
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/nlohmann/json/
)
