/** 
 *  @file
 *  @author Fabian Bösch
 *  opengl includes
 */

#ifndef LB_OPENGL_HPP_INCLUDED
#define LB_OPENGL_HPP_INCLUDED

// OpenGL Graphics includes
#include <GL/glew.h>
#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#ifdef WIN32
     bool IsOpenGLAvailable(const char *appName) { return true; }
#else
  #if (defined(__APPLE__) || defined(MACOSX))
     bool IsOpenGLAvailable(const char *appName) { return true; }
  #else
     // check if this is a linux machine
     #include <X11/Xlib.h>
     bool IsOpenGLAvailable(const char *appName)
     {
        Display *Xdisplay = XOpenDisplay(NULL);
        if (Xdisplay == NULL) {
           return false;
        } else {
           XCloseDisplay(Xdisplay);
           return true;
        }
     }
  #endif
#endif

#if defined(__APPLE__) || defined(MACOSX)
#include <GLUT/glut.h>
#else
#include <GL/freeglut.h>
#endif

#include <stdexcept>

#define GL_CHECK_AND_PRINT_ERROR                                                    \
{                                                                                   \
	GLenum err = glGetError();                                                      \
	if (err != GL_NO_ERROR)                                                         \
	{                                                                               \
		std::cerr << "OPENGL_ERROR: " << __FILE__ << ", line " << __LINE__ << ": "; \
		std::cerr << gluErrorString(err) << " (code " << err << ")" << std::endl;   \
	}                                                                               \
}

#define GL_CHECK_ERROR_AND_THROW                                                    \
{                                                                                   \
	GLenum err = glGetError();                                                      \
	if (err != GL_NO_ERROR)                                                         \
	{                                                                               \
		std::cerr << "OPENGL_ERROR: " << __FILE__ << ", line " << __LINE__ << ": "; \
		std::cerr << gluErrorString(err) << " (code " << err << ")" << std::endl;   \
		throw std::runtime_error("opengl error occured");                           \
	}                                                                               \
}

static GLuint make_buffer(GLenum target, const void *buffer_data, GLsizei buffer_size, GLenum usage = GL_STATIC_DRAW) 
{
    GLuint buffer;
    glGenBuffers(1, &buffer);
    glBindBuffer(target, buffer);
    glBufferData(target, buffer_size, buffer_data, usage);
    return buffer;
}

#endif // LB_OPENGL_HPP_INCLUDED
