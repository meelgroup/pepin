# Copyright (c) 2021, Mate Soos
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${GMP_INCLUDE_DIRS})
include_directories(${GMPXX_INCLUDE_DIR})
include_directories(${MPFR_INCLUDE_DIR})

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/GitSHA1.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp" @ONLY)

add_library (pepin
    pepin-int.cpp
    pepin.cpp
    ${CMAKE_CURRENT_BINARY_DIR}/GitSHA1.cpp
)
target_link_libraries (pepin
  LINK_PUBLIC ${GMPXX_LIBRARY}
  LINK_PUBLIC ${GMP_LIBRARY}
)

add_executable (pepin-bin
    pepin-main.cpp
)

add_executable (pepin-test-bin
    test.cpp
)

set(pepin_exec_link_libs
    ${GMPXX_LIBRARY}
    ${GMP_LIBRARY}
)
IF (ZLIB_FOUND)
    SET(pepin_exec_link_libs ${pepin_exec_link_libs} ${ZLIB_LIBRARY})
ENDIF()

target_link_libraries(pepin-bin
  pepin
  ${pepin_exec_link_libs}
)
target_link_libraries(pepin-test-bin pepin
)

set_target_properties(pepin PROPERTIES
    PUBLIC_HEADER "${pepin_public_headers}"
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    SOVERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
)

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    set_target_properties(pepin-bin PROPERTIES
        OUTPUT_NAME pepin
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE
        LINK_FLAGS "-s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s EXPORTED_RUNTIME_METHODS='[\"callMain\", \"ccall\", \"cwrap\", \"FS\", \"print\"]' -s FORCE_FILESYSTEM=1"
    )
else()
    set_target_properties(pepin-bin PROPERTIES
        OUTPUT_NAME pepin
        RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
        INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

# set_target_properties(pepin-test-bin PROPERTIES
#     OUTPUT_NAME test
#     RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
# )

pepin_add_public_header(
    pepin ${CMAKE_CURRENT_SOURCE_DIR}/pepin.h )

# -----------------------------------------------------------------------------
# Copy public headers into build directory include directory.
# The pepinConfig.cmake we generate in the build directory depends on
# this.
# -----------------------------------------------------------------------------
set(HEADER_DEST "${PROJECT_BINARY_DIR}/include/pepin")
add_custom_target(CopyPublicHeaders ALL)
get_target_property(pepin_public_headers pepin PUBLIC_HEADER)
foreach(public_header ${pepin_public_headers})
    get_filename_component(HEADER_NAME ${public_header} NAME)
    add_custom_command(TARGET CopyPublicHeaders PRE_BUILD
                       COMMAND ${CMAKE_COMMAND} -E make_directory
                               "${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E echo
                       "Copying ${HEADER_NAME} to ${HEADER_DEST}"
                       COMMAND ${CMAKE_COMMAND} -E
                           copy_if_different
                           ${public_header}
                           "${HEADER_DEST}"
                      )
endforeach()


install(TARGETS pepin
    EXPORT ${PEPIN_EXPORT_NAME}
    LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/pepin"
)

install(TARGETS pepin-bin
    EXPORT ${PEPIN_EXPORT_NAME}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
