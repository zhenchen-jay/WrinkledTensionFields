# - Try to find the JSON library
# Once done this will define
#
#  JSON_FOUND - system has JSON
#  JSON_INCLUDE_DIR - **the** JSON include directory
find_path(JSON_INCLUDE_DIR json/json.h
    PATHS
        ${CMAKE_SOURCE_DIR}/../..
        ${CMAKE_SOURCE_DIR}/..
        ${CMAKE_SOURCE_DIR}
        ${CMAKE_SOURCE_DIR}/../tools/json
        ${CMAKE_SOURCE_DIR}/json
        ${CMAKE_SOURCE_DIR}/../json
        ${CMAKE_SOURCE_DIR}/../../json
        ${CMAKE_SOURCE_DIR}/external/json
        ${CMAKE_SOURCE_DIR}/../external/json
        /usr
        /usr/local
        /usr/local/igl/json
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JSON
    "\nJSON not found --- You can download it using:\n\tgit clone --recursivehttps://github.com/nlohmann/json.git ${CMAKE_SOURCE_DIR}/../JSON"
    JSON_INCLUDE_DIR)
mark_as_advanced(JSON_INCLUDE_DIR)
