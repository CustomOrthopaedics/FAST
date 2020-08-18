# Download and set up json

include(cmake/Externals.cmake)

ExternalProject_Add(json
        PREFIX ${FAST_EXTERNAL_BUILD_DIR}/json
        BINARY_DIR ${FAST_EXTERNAL_BUILD_DIR}/json
        #GIT_REPOSITORY "https://github.com/nlohmann/json.git"
        #GIT_TAG "db78ac1d7716f56fc9f1b030b715f872f93964e4"
        URL "https://github.com/nlohmann/json/archive/v3.9.1.tar.gz"
        INSTALL_DIR ${FAST_EXTERNAL_INSTALL_DIR}
        CMAKE_CACHE_ARGS
            -DCMAKE_BUILD_TYPE:STRING=Release
            -DCMAKE_VERBOSE_MAKEFILE:BOOL=OFF
            -DCMAKE_INSTALL_MESSAGE:BOOL=LAZY
            -DCMAKE_INSTALL_PREFIX:STRING=${FAST_EXTERNAL_INSTALL_DIR}
            -DBUILD_TESTING:BOOL=OFF
)

list(APPEND FAST_INCLUDE_DIRS ${FAST_EXTERNAL_INSTALL_DIR}/include/json/)
list(APPEND FAST_EXTERNAL_DEPENDENCIES json)
