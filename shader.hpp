/** 
 *  @file
 *  @author Fabian Bösch
 *  opengl shader
 */

#ifndef LB_SHADER_HPP_INCLUDED
#define LB_SHADER_HPP_INCLUDED

#include "opengl.hpp"
#include <string>
#include <fstream>
#include <vector>

template<typename T>
struct gl_uniform_traits {};

template<>
struct gl_uniform_traits<GLfloat>
{ 
	static void gl_uniform(GLint location, GLfloat v0) { glUniform1f(location,v0); }
	static void gl_uniform(GLint location, GLfloat v0, GLfloat v1) { glUniform2f(location,v0,v1); }
	static void gl_uniform(GLint location, GLfloat v0, GLfloat v1, GLfloat v2) { glUniform3f(location,v0,v1,v2); }
	static void gl_uniform(GLint location, GLfloat v0, GLfloat v1, GLfloat v2, GLfloat v3) { glUniform4f(location,v0,v1,v2,v3); }
	static void gl_uniform1v(GLint location, GLsizei count, const GLfloat *value) { glUniform1fv(location,count,value); }
	static void gl_uniform2v(GLint location, GLsizei count, const GLfloat *value) { glUniform2fv(location,count,value); }
	static void gl_uniform3v(GLint location, GLsizei count, const GLfloat *value) { glUniform3fv(location,count,value); }
	static void gl_uniform4v(GLint location, GLsizei count, const GLfloat *value) { glUniform4fv(location,count,value); }
};

template<>
struct gl_uniform_traits<GLint>
{ 
	static void gl_uniform(GLint location, GLint v0) { glUniform1i(location,v0); }
	static void gl_uniform(GLint location, GLint v0, GLint v1) { glUniform2i(location,v0,v1); }
	static void gl_uniform(GLint location, GLint v0, GLint v1, GLint v2) { glUniform3i(location,v0,v1,v2); }
	static void gl_uniform(GLint location, GLint v0, GLint v1, GLint v2, GLint v3) { glUniform4i(location,v0,v1,v2,v3); }
	static void gl_uniform1v(GLint location, GLsizei count, const GLint *value) { glUniform1iv(location,count,value); }
	static void gl_uniform2v(GLint location, GLsizei count, const GLint *value) { glUniform2iv(location,count,value); }
	static void gl_uniform3v(GLint location, GLsizei count, const GLint *value) { glUniform3iv(location,count,value); }
	static void gl_uniform4v(GLint location, GLsizei count, const GLint *value) { glUniform4iv(location,count,value); }
};

template<>
struct gl_uniform_traits<GLuint>
{ 
	static void gl_uniform(GLint location, GLuint v0) { glUniform1ui(location,v0); }
	static void gl_uniform(GLint location, GLuint v0, GLuint v1) { glUniform2ui(location,v0,v1); }
	static void gl_uniform(GLint location, GLuint v0, GLuint v1, GLuint v2) { glUniform3ui(location,v0,v1,v2); }
	static void gl_uniform(GLint location, GLuint v0, GLuint v1, GLuint v2, GLuint v3) { glUniform4ui(location,v0,v1,v2,v3); }
	static void gl_uniform1v(GLint location, GLsizei count, const GLuint *value) { glUniform1uiv(location,count,value); }
	static void gl_uniform2v(GLint location, GLsizei count, const GLuint *value) { glUniform2uiv(location,count,value); }
	static void gl_uniform3v(GLint location, GLsizei count, const GLuint *value) { glUniform3uiv(location,count,value); }
	static void gl_uniform4v(GLint location, GLsizei count, const GLuint *value) { glUniform4uiv(location,count,value); }
};


class shader
{
public: // ctors and assignement

	shader()
	: vertex_shader(0), fragment_shader(0), shader_program(0)
	{}
	
	~shader()
	{
		if (shader_program != 0) glDeleteProgram(shader_program);
		if (vertex_shader != 0) glDeleteShader(vertex_shader);
		if (fragment_shader != 0) glDeleteShader(vertex_shader);
	}
	
	shader(const shader&) = delete;
	shader& operator=(const shader&) = delete;
	
public: // load shaders from file

	void load(std::string vshader_file_name, std::string fshader_file_name)
	{
		load(std::vector<std::string>(1,vshader_file_name),std::vector<std::string>(1,fshader_file_name));
	}
	
	void load(std::vector<std::string> vshader_file_names, std::string fshader_file_name)
	{
		load(vshader_file_names,std::vector<std::string>(1,fshader_file_name));
	}
	
	void load(std::string vshader_file_name, std::vector<std::string> fshader_file_names)
	{
		load(std::vector<std::string>(1,vshader_file_name),fshader_file_names);
	}
	
	void load(std::vector<std::string> vshader_file_names, std::vector<std::string> fshader_file_names)
	{
		// delete shaders and program if existing
		if (shader_program != 0) glDeleteProgram(shader_program);
		if (vertex_shader != 0) glDeleteShader(vertex_shader);
		if (fragment_shader != 0) glDeleteShader(vertex_shader);
		
		vertex_shader_file_names = vshader_file_names;
		fragment_shader_file_names = fshader_file_names;
		
		// create both shader objects
		vertex_shader = glCreateShader(GL_VERTEX_SHADER);
		fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);

		// load source code from file an add to shader objects
		std::string source_code;
		for (std::string str : vertex_shader_file_names) source_code += read_file(str);
		if(source_code.size() == 0) throw std::runtime_error("vertex shader object could not be created!");
		GLint length = source_code.size();
		const GLchar* source = &source_code[0];
		glShaderSource(vertex_shader, 1, &source, &length);

		source_code = "";
		for (std::string str : fragment_shader_file_names) source_code += read_file(str);
		if(source_code.size() == 0) throw std::runtime_error("fragment shader object could not be created!");
		length = source_code.size();
		source = &source_code[0];
		glShaderSource(fragment_shader, 1, &source, &length);

		// create Program object
		shader_program = glCreateProgram();

		// attach both shaders
		glAttachShader(shader_program, vertex_shader);
		glAttachShader(shader_program, fragment_shader);
		
		GLint status = 0;
		
		// compile both shaders
		glCompileShader(vertex_shader);
		glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &status);
		if (status == GL_FALSE)
		{
			show_info_log(vertex_shader, glGetShaderiv, glGetShaderInfoLog);
			glDeleteShader(vertex_shader);
			throw std::runtime_error("failed to compile vertex shader!");
		}

		glCompileShader(fragment_shader);
		glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &status);
		if(status == GL_FALSE)
		{
			show_info_log(fragment_shader, glGetShaderiv, glGetShaderInfoLog);
			glDeleteShader(vertex_shader);
			glDeleteShader(fragment_shader);
			throw std::runtime_error("failed to compile fragment shader!");
		}

		// Link Program
		glLinkProgram(shader_program);
		glGetProgramiv(shader_program, GL_LINK_STATUS, &status);
		if(status == GL_FALSE)
		{
			show_info_log(shader_program, glGetProgramiv, glGetProgramInfoLog);
			glDeleteShader(vertex_shader);
			glDeleteShader(fragment_shader);
			glDeleteProgram(shader_program);
			throw std::runtime_error("failed to link shader program!");
		}
	}

public: // activate

	void activate()
	{
		GLint status = 0;
		glGetProgramiv(shader_program, GL_LINK_STATUS, &status);
		if(status == GL_FALSE) throw std::runtime_error("shaders can only be used after successful linking!");
		glUseProgram(shader_program);
	}

	void deactivate()
	{
		// reactivate fixed function pipeline
		// (note: ffp is deprecated in opengl3.0 and not implmented in opengles 2.0
		glUseProgram(0);
	}
	
public: // get location of a uniform or attribute

	GLint get_uniform_location(std::string name) const
	{
		GLint return_value = glGetUniformLocation(shader_program, name.c_str());
		if (return_value == -1) { std::cout << name << std::endl; throw std::runtime_error("uniform does not exist"); }
		GL_CHECK_ERROR_AND_THROW
		return return_value;
	}
	
	GLint get_attribute_location(std::string name) const
	{
		GLint return_value = glGetAttribLocation(shader_program, name.c_str());
		if (return_value == -1) throw std::runtime_error("attribute does not exist");
		GL_CHECK_ERROR_AND_THROW
		return return_value;
	}
	
public: // set and get uniforms

	template <typename T>
	void set_uniform(GLint location, T v0) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(location,v0);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(std::string name, T v0) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(get_uniform_location(name),v0);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(GLint location, T v0, T v1) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(location,v0,v1);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(std::string name, T v0, T v1) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(get_uniform_location(name),v0,v1);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(GLint location, T v0, T v1, T v2) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(location,v0,v1,v2);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(std::string name, T v0, T v1, T v2) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(get_uniform_location(name),v0,v1,v2);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(GLint location, T v0, T v1, T v2, T v3) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(location,v0,v1,v2,v3);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform(std::string name, T v0, T v1, T v2, T v3) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform(get_uniform_location(name),v0,v1,v2,v3);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform1v(GLint location, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform1v(location,count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform1v(std::string name, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform1v(get_uniform_location(name),count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform2v(GLint location, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform2v(location,count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform2v(std::string name, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform2v(get_uniform_location(name),count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform3v(GLint location, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform3v(location,count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform3v(std::string name, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform3v(get_uniform_location(name),count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform4v(GLint location, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform4v(location,count,value);
		GL_CHECK_ERROR_AND_THROW
	}
	
	template <typename T>
	void set_uniform4v(std::string name, GLsizei count, const T* value) 
	{
		program_context p_ctxt(shader_program);
		gl_uniform_traits<T>::gl_uniform4v(get_uniform_location(name),count,value);
		GL_CHECK_ERROR_AND_THROW
	}

private: // impl

	std::string read_file(std::string file_name)
	{
		// open file
		std::ifstream ifs(file_name, std::ios::binary);
		if (!ifs)
		{
			std::cerr << "File not found " << file_name << std::endl;
			return "";
		}
		
		// get length of file
		ifs.seekg(0, std::ios::end);
		std::size_t file_length = ifs.tellg();
		ifs.seekg(0, std::ios::beg);
		
		// allocate string buffer and read string from file
		char* str_ptr = new char[file_length + 1];
		ifs.read(str_ptr, file_length);
		str_ptr[file_length] = 0;
		
		// convert to STL string
		std::string str(str_ptr);
		delete[] str_ptr;
		return str;
	}
	
	void show_info_log(GLuint object, PFNGLGETSHADERIVPROC glGet__iv, PFNGLGETSHADERINFOLOGPROC glGet__InfoLog) const
	{
		GLint log_length;
		glGet__iv(object, GL_INFO_LOG_LENGTH, &log_length);
		char *log = new char[log_length];
		glGet__InfoLog(object, log_length, NULL, log);
		std::cerr << log << std::endl;
		delete[] log;
	}
	
private: // impl

	// Uses a program as long as this object exists
	class program_context
	{
	public:
		program_context(GLint programID) 
		{
			glGetIntegerv(GL_CURRENT_PROGRAM, &m_oldPogramID);
			glUseProgram(programID);
		}

		~program_context() 
		{
			glUseProgram(m_oldPogramID);
		}

	private:
		GLint m_oldPogramID;
	};
	
	
private: // members
	
	GLuint vertex_shader;
	GLuint fragment_shader;
	GLuint shader_program;
	
	std::vector<std::string> vertex_shader_file_names;
	std::vector<std::string> fragment_shader_file_names;
	
};


#endif // LB_SHADER_HPP_INCLUDED
