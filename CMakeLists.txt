#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8.11)

#Name your project here
project(NetGen)

#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")

#find the STIM library
find_package(STIM)

#if the STIM library isn't found, download it
if(NOT STIM_INCLUDE_DIRS)
	file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/stimlib)	#remove the stimlib directory if it exists
	set(STIM_GIT "https://git.stim.ee.uh.edu/codebase/stimlib.git")
	execute_process(COMMAND git clone --depth 1 ${STIM_GIT} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
	set(STIM_INCLUDE_DIRS "${CMAKE_BINARY_DIR}/stimlib" CACHE TYPE PATH)
endif(NOT STIM_INCLUDE_DIRS)

#find BOOST
find_package(Boost REQUIRED)

#find cuda
find_package(CUDA REQUIRED)

#find the GLUT library for visualization
find_package(OpenGL REQUIRED)

find_package(GLUT REQUIRED)

#find OpenCV
find_package(OpenCV REQUIRED)
add_definitions(-DUSING_OPENCV)

#find the pthreads package
find_package(Threads)

#find the X11 package
find_package(X11)

#find the Approximate Nearest Neighbor Library
#find_package(ANN REQUIRED)

#build the executable in the binary directory on MS Visual Studio
if ( MSVC )
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
endif ( MSVC )

include_directories(
					${OPENGL_INCLUDE_DIRS}
					${GLUT_INCLUDE_DIR}
					${STIM_INCLUDE_DIRS}
          			#${ANN_INCLUDE_DIR}
					${Boost_INCLUDE_DIR}
					${CUDA_INCLUDE_DIRS}
					${OpenCV_INCLUDE_DIRS}
					)

#Assign source files to the appropriate variables
#file(GLOB SRC_CU "*.cu")
#file(GLOB SRC_H "*.h")

#create an executable file
cuda_add_executable(NetGen main.cpp)

#set the link libraries
target_link_libraries(NetGen
#						${OpenGL_LIBRARIES}
						${OPENGL_gl_LIBRARY}
						${OPENGL_glu_LIBRARY}
						${GLUT_LIBRARIES}
						${CMAKE_THREAD_LIBS_INIT}
            			#${ANN_LIBRARY}
						${X11_LIBRARIES}
						#${CUDA_cublas_LIBRARY}
						${OpenCV_LIBS}
		  )
