
set(HEADERS
    include/json/autolink.h
    include/json/config.h
    include/json/features.h
    include/json/forwards.h
    include/json/json.h
    include/json/reader.h
    include/json/value.h
    include/json/writer.h

    src/lib_json/json_batchallocator.h
)

set(INLINES
    src/lib_json/json_internalarray.inl
    src/lib_json/json_internalmap.inl
    src/lib_json/json_valueiterator.inl

)

set(SOURCES
    src/lib_json/json_reader.cpp
    src/lib_json/json_value.cpp
    src/lib_json/json_writer.cpp
)

# Set include directories.
include_directories(
    include/
)

# Disable stupid MSVC compiler warnings.
add_definitions(-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS)

set(CMAKE_DEBUG_POSTFIX "d")

add_library(json STATIC
    ${HEADERS}
    ${INLINES}
    ${SOURCES}
)

# propnet has the following dependencies.
#add_dependencies(propnet
#)

# Set link libraries.
#target_link_libraries(propnet
#    ${Boost_LIBRARIES}
#)

# install(TARGETS propnet
#         RUNTIME DESTINATION bin
#         LIBRARY DESTINATION lib
#         ARCHIVE DESTINATION lib)
#
# #install(DIRECTORY include
#         DESTINATION .
#         PATTERN "*.hpp")
