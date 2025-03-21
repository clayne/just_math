cmake_minimum_required(VERSION 2.8)
set(PROJNAME math_3ddda)
Project(${PROJNAME})
Message(STATUS "-------------------------------")
Message(STATUS "Processing Project ${PROJNAME}:")

#####################################################################################
# LIBMIN Installed Path
#
set ( LIBMIN_PATH "d:/codes/build/libmin" CACHE STRING "Location of Libmin")
list( APPEND CMAKE_MODULE_PATH "${LIBMIN_PATH}/cmake" )
list( APPEND CMAKE_PREFIX_PATH "${LIBMIN_PATH}/cmake" )

#####################################################################################
# Include LIBMIN
#
find_package(Libmin QUIET)

if (NOT LIBMIN_FOUND)

  Message ( FATAL_ERROR "
  This project requires Libmin. Build Libmin first.
  Set the LIBMIN_PATH to the installed location (not the libmin source).
  Libmin also provides the cmakes used to build this app. 
  " )

else()
  add_definitions(-DUSE_LIBMIN)  
  include_directories(${LIBMIN_INC_DIR})
  include_directories(${LIBRARIES_INC_DIR})  

  if (DEFINED ${BUILD_LIBMIN_STATIC})
    add_definitions(-DLIBMIN_STATIC) 
    file(GLOB LIBMIN_SRC "${LIBMIN_SRC_DIR}/*.cpp" )
    file(GLOB LIBMIN_INC "${LIBMIN_INC_DIR}/*.h" )
    LIST( APPEND LIBMIN_SOURCE_FILES ${LIBMIN_SRC} ${LIBMIN_INC} )
    message ( STATUS "  ---> Using LIBMIN (static)")
  else()    
    LIST( APPEND LIBRARIES_OPTIMIZED "${LIBMIN_LIB_DIR}/${LIBMIN_REL}")
    LIST( APPEND LIBRARIES_DEBUG "${LIBMIN_LIB_DIR}/${LIBMIN_DEBUG}")	    
    _EXPANDLIST( OUTPUT PACKAGE_DLLS SOURCE ${LIBMIN_LIB_DIR} FILES ${LIBMIN_DLLS} )
    message ( STATUS "  ---> Using LIBMIN")
  endif() 
endif()

#####################################################################################
# Options

set( USE_NVTX true CACHE BOOL "")

_REQUIRE_MAIN()

_REQUIRE_GL()

_REQUIRE_GLEW()

#--- symbols in release mode
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Zi" CACHE STRING "" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG /OPT:REF /OPT:ICF" CACHE STRING "" FORCE)

######################
# CMAKE_INSTALL_PREFIX -- path where library will be installed to

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
   if ( WIN32 )
      get_filename_component ( _instpath "${CMAKE_CURRENT_BINARY_DIR}" REALPATH )
   else()
      get_filename_component ( _instpath "/usr/local/shapes" REALPATH )
   endif()
   set ( CMAKE_INSTALL_PREFIX ${_instpath} CACHE PATH "default install path" FORCE)   
endif()

#####################################################################################
# Asset Path
#
if ( NOT DEFINED ASSET_PATH ) 
   get_filename_component ( _assets "${BASE_DIRECTORY}/assets" REALPATH )
   set ( ASSET_PATH ${_assets} CACHE PATH "Full path to /assets" )
   add_definitions( -DASSET_PATH="${ASSET_PATH}/" )
endif()
add_definitions(-DASSET_PATH="${ASSET_PATH}/")

#####################################################################################
# Executable
#
file(GLOB MAIN_FILES *.cpp *.c *.h )

unset ( ALL_SOURCE_FILES )

list( APPEND ALL_SOURCE_FILES ${MAIN_FILES} )
list( APPEND ALL_SOURCE_FILES ${COMMON_SOURCE_FILES} )
list( APPEND ALL_SOURCE_FILES ${PACKAGE_SOURCE_FILES} )
list( APPEND ALL_SOURCE_FILES ${UTIL_SOURCE_FILES} )

if ( NOT DEFINED WIN32 )
    set(libdeps GL GLEW X11)
  LIST(APPEND LIBRARIES_OPTIMIZED ${libdeps})
  LIST(APPEND LIBRARIES_DEBUG ${libdeps})
ENDIF()
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")    

add_executable (${PROJNAME} ${ALL_SOURCE_FILES} ${GLSL_FILES} )

set_property ( TARGET ${PROJNAME} APPEND PROPERTY DEPENDS )

#--- debug and release exe
set ( CMAKE_DEBUG_POSTFIX "d" CACHE STRING "" )
set_target_properties( ${PROJNAME} PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

#####################################################################################
# Additional Libraries
#
_LINK ( PROJECT ${PROJNAME} OPT ${LIBRARIES_OPTIMIZED} DEBUG ${LIBRARIES_DEBUG} PLATFORM ${PLATFORM_LIBRARIES} )

#####################################################################################
# Windows specific
#
_MSVC_PROPERTIES()
source_group("Source Files" FILES ${MAIN_FILES} ${COMMON_SOURCE_FILES} ${PACKAGE_SOURCE_FILES})

#####################################################################################
# Install Binaries
#
_DEFAULT_INSTALL_PATH()

# *NOTE*: file COPY is at cmake-time, not compile-time. Need to replace with add_custom_command -E copy_directory (my own _COPY)

file (COPY "${CMAKE_CURRENT_SOURCE_DIR}/assets" DESTINATION ${CMAKE_INSTALL_PREFIX} )	# assets folder
_INSTALL ( FILES ${GLEW_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX} )
_INSTALL ( FILES ${SHADERS} DESTINATION "${CMAKE_INSTALL_PREFIX}/assets" )		# shaders
_INSTALL ( FILES ${PACKAGE_DLLS} DESTINATION ${CMAKE_INSTALL_PREFIX} )			# DLLs
# _INSTALL_PTX ( FILES ${CUDA_PTX_PATHS} DESTINATION ${CMAKE_INSTALL_PREFIX} OUTPUT INSTALL_LIST )
install ( FILES $<TARGET_PDB_FILE:${PROJNAME}> DESTINATION ${CMAKE_INSTALL_PREFIX} OPTIONAL )		# PDB

install ( FILES ${INSTALL_LIST} DESTINATION ${EXECUTABLE_OUTPUT_PATH} )		# exe, pdb

###########################
# Done
message ( STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}" )
message ( STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}" )
message ( STATUS "LIBS OPTIMIZED: ${LIBRARIES_OPTIMIZED}" )
message ( STATUS "LIBS DEBUG:     ${LIBRARIES_DEBUG}" )
message ( STATUS "PACKAGED DLLS:  ${PACKAGE_DLLS}" )
message ( STATUS "INSTALLED PATH: ${CMAKE_INSTALL_PREFIX}" )
