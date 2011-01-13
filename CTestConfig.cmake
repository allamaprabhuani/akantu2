## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "Akantu")
set(CTEST_NIGHTLY_START_TIME "06:10:00 EST")
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "lsmssrv1.epfl.ch")
set(CTEST_DROP_LOCATION "/CDash/submit.php?project=Akantu")
set(CTEST_DROP_SITE_CDASH TRUE)
