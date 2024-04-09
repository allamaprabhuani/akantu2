find_package(Git)

set(CURRENT_WORKDIR @WORKING_DIRECTORY@)

set(AKANTU_DIRTY_PATCH)
if (Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} diff @
    WORKING_DIRECTORY ${CURRENT_WORKDIR}
    OUTPUT_FILE ${CURRENT_WORKDIR}/git_patch
    ERROR_QUIET
  )
  file(STRINGS ${CURRENT_WORKDIR}/git_patch _res)
  set(_res_escaped)
  foreach(str ${_res})
    list(APPEND _res_escaped "R\"(${str})\"")
  endforeach()
  string(REPLACE ";" ", " _res_escaped "${_res_escaped}")
  set(AKANTU_DIRTY_PATCH "${_res_escaped}")
endif()

set(AKANTU_VERSION @AKANTU_VERSION@)

configure_file(common/aka_config.cc.in
  "${CURRENT_WORKDIR}/aka_config.cc" @ONLY)
