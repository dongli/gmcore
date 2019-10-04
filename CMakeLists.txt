cmake_minimum_required(VERSION 3.0)

project(gmcore LANGUAGES Fortran)

option(STAGGER_V_ON_POLE "Stagger V on the two Poles" OFF)

if (STAGGER_V_ON_POLE)
  add_definitions(-DSTAGGER_V_ON_POLE)
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffpe-summary=none -fdefault-real-8 -ffree-line-length-none -Ofast -fopenmp")
endif ()

find_package(MPI)
include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

add_subdirectory(lib/container)
add_subdirectory(lib/datetime)
add_subdirectory(lib/fiona)
add_subdirectory(lib/flogger)

set(source_files
  src/utils/const_mod.F90
  src/utils/namelist_mod.F90
  src/utils/time_mod.F90
  src/utils/math_mod.F90
  src/utils/history_mod.F90 
  src/utils/time_mod.F90
  src/utils/sphere_geometry_mod.F90
  src/utils/allocator_mod.F90
  src/utils/parallel_mod.F90
  src/utils/debug_mod.F90
  src/core/mesh_mod.F90
  src/core/state_mod.F90
  src/core/static_mod.F90
  src/core/tend_mod.F90
  src/core/pv_mod.F90
  src/core/ke_mod.F90
  src/core/operators_mod.F90
  src/core/reduce_mod.F90
  src/core/gmcore_mod.F90
  src/tests/swm/rossby_haurwitz_wave_test_mod.F90
  src/tests/swm/mountain_zonal_flow_test_mod.F90
  src/tests/swm/jet_zonal_flow_test_mod.F90
  src/tests/swm/steady_geostrophic_flow_test_mod.F90
  lib/quadpack.f90
)

add_library(gmcore ${source_files})
target_link_libraries(gmcore netcdff fortran_container fortran_datetime fiona flogger ${MPI_Fortran_LIBRARIES})

add_executable(gmcore_swm_driver.exe src/drivers/gmcore_swm_driver.F90)
target_link_libraries(gmcore_swm_driver.exe gmcore)
