#---------------------------------------------------------------------------
# Project related configuration options
#---------------------------------------------------------------------------

DOXYFILE_ENCODING      = UTF-8
PROJECT_NAME           = @CMAKE_PROJECT_NAME@
PROJECT_NUMBER         = @AKANTU_VERSION@
OUTPUT_DIRECTORY       = .
OUTPUT_LANGUAGE        = English
RECURSIVE              = YES
FULL_PATH_NAMES        = YES
STRIP_FROM_PATH        = @CMAKE_SOURCE_DIR@
STRIP_FROM_INC_PATH    = @CMAKE_SOURCE_DIR@
TAB_SIZE               = 4
BUILTIN_STL_SUPPORT    = NO

#---------------------------------------------------------------------------
# configuration options related to warning and progress messages
#---------------------------------------------------------------------------

QUIET                  = @DOXYGEN_QUIET@
WARNINGS               = @DOXYGEN_WARNINGS@
WARN_IF_UNDOCUMENTED   = @DOXYGEN_WARNINGS@
WARN_IF_DOC_ERROR      = @DOXYGEN_WARNINGS@

#---------------------------------------------------------------------------
# configuration options related to the input files
#---------------------------------------------------------------------------

INPUT                  = @CMAKE_SOURCE_DIR@/src
EXAMPLE_PATH           = @CMAKE_SOURCE_DIR@/test

#---------------------------------------------------------------------------
# configuration options related to source browsing
#---------------------------------------------------------------------------

SOURCE_BROWSER         = YES
INLINE_SOURCES         = NO

#---------------------------------------------------------------------------
# Configuration options related to the preprocessor
#---------------------------------------------------------------------------

ENABLE_PREPROCESSING   = YES
MACRO_EXPANSION        = YES
EXPAND_ONLY_PREDEF     = NO
SEARCH_INCLUDES        = YES
INCLUDE_PATH           = @CMAKE_SOURCE_DIR@/src/common \
		         @CMAKE_SOURCE_DIR@/src/synchronizer \
		         @CMAKE_SOURCE_DIR@/src/solver \
			 @CMAKE_SOURCE_DIR@/src/fem \
			 @CMAKE_SOURCE_DIR@/src/fem/element_classes \
			 @CMAKE_SOURCE_DIR@/src/mesh_utils \
			 @CMAKE_SOURCE_DIR@/src/mesh_utils/mesh_io \
			 @CMAKE_SOURCE_DIR@/src/mesh_utils/mesh_partition \
			 @CMAKE_SOURCE_DIR@/src/model \
			 @CMAKE_SOURCE_DIR@/src/model/integration_scheme \
			 @CMAKE_SOURCE_DIR@/src/model/solid_mechanics \
			 @CMAKE_SOURCE_DIR@/src/model/solid_mechanics/materials \
			 @CMAKE_SOURCE_DIR@/src/model/solid_mechanics/contact

INCLUDE_FILE_PATTERNS  = *.hh
PREDEFINED             = @AKANTU_DOXYGEN_DEFINTIONS@
EXPAND_AS_DEFINED      = __BEGIN_AKANTU__ \
                         __END_AKANTU__ \
                         AKANTU_GET_MACRO \
                         AKANTU_SET_MACRO \
                         AKANTU_GET_MACRO_NOT_CONST \
                         AKANTU_GET_MACRO_BY_ELEMENT_TYPE

SKIP_FUNCTION_MACROS   = YES

#---------------------------------------------------------------------------
# Configuration options related to the dot tool
#---------------------------------------------------------------------------

CLASS_DIAGRAMS         = YES
HAVE_DOT               = YES
CLASS_GRAPH            = YES
COLLABORATION_GRAPH    = YES
TEMPLATE_RELATIONS     = YES
CALL_GRAPH             = YES
CALLER_GRAPH           = YES
DOT_PATH               = @DOXYGEN_DOT_PATH@
