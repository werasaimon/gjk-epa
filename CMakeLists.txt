
project(gjk-epa)
cmake_minimum_required(VERSION 2.8.12)
#set(CMAKE_PREFIX_PATH "/Users/andrea/Developer/Qt5/5.2.1/clang_64")

## Qt library
#find_package(Qt5Widgets)
#find_package(Qt5OpenGL)
#set(CMAKE_AUTOMOC ON)

# OpenGL & GLEW library
set(QT_USE_QTOPENGL TRUE)
find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

# executable
add_executable( gjk-epa main.cpp
                freeglut/freeglut_xinput.c
                freeglut/freeglut_window.c
                freeglut/freeglut_videoresize.c
                freeglut/freeglut_teapot.c
                freeglut/freeglut_structure.c
                freeglut/freeglut_stroke_roman.c
                freeglut/freeglut_stroke_mono_roman.c
                freeglut/freeglut_state.c
                freeglut/freeglut_spaceball.c
                freeglut/freeglut_overlay.c
                freeglut/freeglut_misc.c
                freeglut/freeglut_menu.c
                freeglut/freeglut_main.c
                freeglut/freeglut_joystick.c
                freeglut/freeglut_input_devices.c
                freeglut/freeglut_init.c
                freeglut/freeglut_glutfont_definitions.c
                freeglut/freeglut_geometry.c
                freeglut/freeglut_gamemode.c
                freeglut/freeglut_font_data.c
                freeglut/freeglut_font.c
                freeglut/freeglut_ext.c
                freeglut/freeglut_display.c
                freeglut/freeglut_cursor.c
                freeglut/freeglut_callbacks.c
                gjk_epa.hpp)




#qt5_use_modules(gjk-epa Widgets OpenGL)
target_link_libraries(gjk-epa ${QT_LIBRARIES} ${OPENGL_LIBRARIES})

include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
set(OPENGL_LIBS ${OPENGL_LIBRARIES} ${GLUT_LIBRARY})

target_compile_options( gjk-epa PUBLIC -std=c++11 -fpermissive -w -Wall )

target_link_libraries( gjk-epa ${OpenCV_LIBS}  -lX11 )
