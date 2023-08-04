#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "HYPRE::HYPRE" for configuration "Release"
set_property(TARGET HYPRE::HYPRE APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(HYPRE::HYPRE PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "C"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libHYPRE.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS HYPRE::HYPRE )
list(APPEND _IMPORT_CHECK_FILES_FOR_HYPRE::HYPRE "${_IMPORT_PREFIX}/lib64/libHYPRE.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
