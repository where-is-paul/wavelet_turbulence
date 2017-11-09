cmake_minimum_required( VERSION 2.8 )

message(STATUS "*** Wavelet Turbulence ***")

find_package( ZLIB REQUIRED )
find_package( PNG REQUIRED )

include_directories( . ${ZLIB_INCLUDE_DIRS} ${PNG_INCLUDE_DIRS} )

set( WTURB_HEADERS
    FLUID_3D.h
	EIGENVALUE_HELPER.h
	IMAGE.h
	INTERPOLATE.h
	MERSENNETWISTER.h
	OBSTACLE.h
	SPHERE.h
	VEC3.h
	WAVELET_NOISE.h
	WTURBULENCE.h
	GABOR_NOISE.h
)

set( WTURB_SOURCES
    FLUID_3D.cpp
	FLUID_3D_SOLVERS.cpp
	FLUID_3D_STATIC.cpp
	main.cpp
	#noiseFFT.cpp
	SPHERE.cpp
	WTURBULENCE.cpp
)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable( FLUID_3D ${WTURB_SOURCES} ${WTURB_HEADERS} )
target_link_libraries(FLUID_3D ${ZLIB_LIBRARIES} ${PNG_LIBRARIES})
