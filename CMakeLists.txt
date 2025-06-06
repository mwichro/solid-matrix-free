##
#  CMake script for the step-<XYZ> tutorial program:
##

OPTION(WITH_BREAKDOWN "Hack breackdown measurements of MF vmult() with LIKWID" OFF)

# Set the name of the project and target:
SET(TARGET "main")
SET(TARGETLIB "mf_elasticity")

INCLUDE_DIRECTORIES(
 ${CMAKE_SOURCE_DIR}/include
)

# Declare all source files the target consists of:
SET(TARGETLIB_SRC
  src/mf_elasticity.cc
  ${CMAKE_CURRENT_BINARY_DIR}/version.cc
)

SET(CLEAN_UP_FILES
  # a custom list of globs, e.g. *.log *.vtk
  *.vtk
)

# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 9.2 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

# Run macro to get Git info:
DEAL_II_QUERY_GIT_INFORMATION()
# If we could not get the tag from Git, get it from VERSION file
IF(NOT GIT_TAG)
  FILE(STRINGS "${CMAKE_SOURCE_DIR}/VERSION" GIT_TAG LIMIT_COUNT 1)
ENDIF()
MESSAGE("-- Code version ${GIT_TAG}")
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/src/version.cc.in
               ${CMAKE_CURRENT_BINARY_DIR}/version.cc)

PROJECT(${TARGET})

#
# RPATHS: https://cmake.org/Wiki/CMake_RPATH_handling#Mac_OS_X_and_the_RPATH
#
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  ## cmake  --help-policy CMP0042
  set(CMAKE_MACOSX_RPATH ON)
ENDIF()


# Specify the directory where the model header files are located
set(MODEL_DIR "${CMAKE_SOURCE_DIR}/include/models")

# Find all the header files in the models directory
file(GLOB MODEL_FILES "${MODEL_DIR}/*.h")

# Create a list of filenames for users to choose from
list(TRANSFORM MODEL_FILES PREPEND "models/")
set(MODEL_OPTIONS "")
foreach(file ${MODEL_FILES})
    get_filename_component(file_name ${file} NAME_WE)
    list(APPEND MODEL_OPTIONS ${file_name})
endforeach()

message(STATUS "Available model files: ${MODEL_OPTIONS}")
message(STATUS "Use -DSELECTED_MODEL to select different one")

if(DEFINED SELECTED_MODEL)
    # Print the selected model for verification
    message(STATUS "Selected model: ${SELECTED_MODEL}")
    # Check if SELECTED_MODEL is in the MODEL_OPTIONS list
    list(FIND MODEL_OPTIONS ${SELECTED_MODEL} SELECTED_INDEX)

    if(${SELECTED_INDEX} EQUAL -1)
      message(FATAL_ERROR "Selected model '${SELECTED_MODEL}' is not in the list of available options.")
    else()
      message(STATUS "Model ${SELECTED_MODEL} is available")
    endif()

else()
    list(GET MODEL_OPTIONS 0 SELECTED_MODEL)
    set(SELECTED_MODEL ${MODEL_OPTIONS} CACHE STRING "Select the model header to include")
endif()

# Print the selected model for verification


add_definitions(-DMODEL_INCLUDE=${SELECTED_MODEL})

# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# make sure RPATHs are used on macOS
SET(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

ADD_LIBRARY(${TARGETLIB} SHARED
    ${TARGETLIB_SRC}
  )

ADD_EXECUTABLE(${TARGET}
    src/main.cc
  )

DEAL_II_SETUP_TARGET(${TARGETLIB})
DEAL_II_SETUP_TARGET(${TARGET})

TARGET_LINK_LIBRARIES(${TARGET}
    ${TARGETLIB}
  )

#
# Likwid, configure with:
# -DWITH_LIKWID=TRUE -DLIKWID_DIR=/path/to/prefix
#
OPTION(WITH_LIKWID "Enable LIKWID" OFF)


IF (WITH_LIKWID)
set(LIKWID_DIR "/usr/local/" CACHE PATH "Path to LIKWID installation directory")
  FIND_LIBRARY(LIKWID_LIBRARY
    NAMES likwid
    HINTS ${LIKWID_DIR}/lib
    NO_DEFAULT_PATH
  )
  IF(LIKWID_LIBRARY STREQUAL "LIKWID_LIBRARY-NOTFOUND")
    MESSAGE("-- Build without Likwid")
    SET(WITH_LIKWID "FALSE")
  ELSE()
    MESSAGE("-- Use Likwid from ${LIKWID_DIR}")
    TARGET_LINK_LIBRARIES(${TARGETLIB}
      "${LIKWID_LIBRARY}"
    )
    # TODO: use FIND_PATH() to look for the headers
    INCLUDE_DIRECTORIES ("${LIKWID_DIR}/include")

    # Finally, add necessary flags:
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DLIKWID_PERFMON")
  ENDIF()
ELSE()
  MESSAGE("-- Build without Likwid")
ENDIF ()

INCLUDE_DIRECTORIES (${CMAKE_CURRENT_BINARY_DIR}/include)

#
# Create a config file
#
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/include/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/include/config.h)

#
# Custom "debug" and "release" make targets:
#
ADD_CUSTOM_TARGET(debug
COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
)

ADD_CUSTOM_TARGET(release
COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
COMMENT "Switch CMAKE_BUILD_TYPE to Release"
)

# Benchmarks
ADD_SUBDIRECTORY(benchmark)

# Enable testing and descent into tests/ subdirectory:
ENABLE_TESTING()
INCLUDE(CTest)
ADD_SUBDIRECTORY(tests)
