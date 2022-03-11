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
        ${CMAKE_SOURCE_DIR}/../tools/Jsoncpp
        ${CMAKE_SOURCE_DIR}/Jsoncpp
        ${CMAKE_SOURCE_DIR}/../Jsoncpp
        ${CMAKE_SOURCE_DIR}/../../Jsoncpp
        ${CMAKE_SOURCE_DIR}/external/Jsoncpp
        ${CMAKE_SOURCE_DIR}/../external/Jsoncpp
        /usr
        /usr/local
        /usr/local/igl/Jsoncpp
    PATH_SUFFIXES include
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(JSON
    "\nJSON not found --- You can download it using:\n\tgit clone --recursive https://github.com/open-source-parsers/jsoncpp.git ${CMAKE_SOURCE_DIR}/../JSON"
    JSON_INCLUDE_DIR)
mark_as_advanced(JSON_INCLUDE_DIR)

list(APPEND CMAKE_MODULE_PATH "${JSON_INCLUDE_DIR}/../cmake")
include(JSON)
