/** 
 *  @file
 *  @author Fabian Bösch
 *  opengl visualization
 */

#ifndef LB_VISUALIZATION_HPP_INCLUDED
#define LB_VISUALIZATION_HPP_INCLUDED

#include "simulation.hpp"
#include "opengl.hpp"
#include "shader.hpp"
#include <memory> // unique_ptr
#include <typeinfo> // typeid
#include <cmath>
#include <vector>
#include <list>
#include <cstdlib> // rand
#include <sstream>
#include "third_party/EasyBMP/EasyBMP.h"

namespace lb {

template <typename T>
inline T scale(T value, T min_value, T max_value)
{
	return (value - min_value) / (max_value - min_value);
}

class visualization // singleton
{
public: // ctors
	
	visualization(const visualization& other) = delete;
	visualization& operator=(const visualization& other) = delete;
	
private: // ctors

	visualization(simulation* _sim, int argc, char *argv[])
	: // window properties
	  panel_size(170),
	  width(((800-panel_size)/_sim->l.nx)*(_sim->l.nx)+panel_size),
	  height(((800-panel_size)/_sim->l.nx)*(_sim->l.ny)),
	  //width(std::min(static_cast<int>(600+panel_size), static_cast<int>(_sim->l.nx+panel_size))), 
	  //height(std::min(static_cast<int>(600),static_cast<int>( _sim->l.ny))),
	  // pointer to simulation
	  sim(_sim),
	  // simulation properties
	  dim_x(sim->l.nx),
	  dim_y(sim->l.ny),
	  buffer_size(sim->l.buffer_size),
	  real_dim_x(dim_x + 2*buffer_size),
	  real_dim_y(dim_y + 2*buffer_size),
	  periodic_x(sim->l.periodic_x),
	  periodic_y(sim->l.periodic_y),
	  // if simulation does not use floats -> cast every time
	  float_cast(((typeid(float_type)==typeid(float)) ? false : true)),
	  rho_data_cast((float_cast ? real_dim_x*real_dim_y : 0),0),
	  u_data_cast((float_cast ? real_dim_x*real_dim_y : 0),0),
	  v_data_cast((float_cast ? real_dim_x*real_dim_y : 0),0),
	  // wall container
	  wall(real_dim_x*real_dim_y,0),
	  // data pointers
	  rho_data((float_cast ? &rho_data_cast[0] : reinterpret_cast<float*>(&(sim->l.rho[0])))),
	  u_data((float_cast ? &u_data_cast[0] : reinterpret_cast<float*>(&(sim->l.u[0])))),
	  v_data((float_cast ? &v_data_cast[0] : reinterpret_cast<float*>(&(sim->l.v[0])))),
	  wall_data(&wall[0]),
	  // opacity value for drawing
	  alpha0(0.8),
	  // ibfv params
	  ibfv_alpha(0.12*255),//(0.12*255),
	  ibfv_texture_dim(64),
	  ibfv_sample_rate_x(std::max(dim_x/100, 1)),
	  ibfv_sample_rate_y(std::max(dim_y/100, 1)),
	  ibfv_additional_pts(5),
	  ibfv_grid_dim_x(dim_x/ibfv_sample_rate_x+2*ibfv_additional_pts),
	  ibfv_grid_dim_y(dim_y/ibfv_sample_rate_y+2*ibfv_additional_pts),
	  ibfv_scale(0.1),
	  ibfv_num_pattern(32),
	  frame_index(0),
	  ibfv_u_data(ibfv_grid_dim_x*ibfv_grid_dim_y,0),
	  ibfv_v_data(ibfv_grid_dim_x*ibfv_grid_dim_y,0),
	  ibfv_x_index_lookup(ibfv_grid_dim_x),
	  ibfv_y_index_lookup(ibfv_grid_dim_y),
	  // opengl buffers (location and sizes)
	  simulation_grid_vbo_size(real_dim_x*real_dim_y*4),
	  simulation_grid_ibo_size((real_dim_x-1)*(real_dim_y-1)*4),
	  field_vbo_size(real_dim_x*real_dim_y),
	  // timings
	  fps_hist(30),
	  fps_hist_index(0),
	  rho_max_hist(30,1.001),
	  rho_min_hist(30,0.999),
	  u_max_hist(30,sim->Vmax),
	  u_min_hist(30,-sim->Vmax),
	  v_max_hist(30,sim->Vmax),
	  v_min_hist(30,-sim->Vmax),
	  vel_mag2_max_hist(30,(sim->Vmax)*(sim->Vmax)),
	  vel_mag2_min_hist(30,0),
	  data_hist_index(0),
	  // limits for simulation data
	  min_rho(rho_min_hist[0]), max_rho(rho_max_hist[0]), 
	  min_u(u_min_hist[0]), max_u(u_max_hist[0]), 
	  min_v(v_min_hist[0]), max_v(v_max_hist[0]),
	  min_vel2(vel_mag2_min_hist[0]), max_vel2(vel_mag2_max_hist[0]),
	  // sate flags
	  running(0),
	  num_modes(8),
	  mode(0),
	  start_particles(true),
	  start_timeline(false),
	  picture_index(0),
	  wall_drawing_mode(false)
	{
		glutInit(&argc, argv);
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);
		glutInitWindowPosition(100,100);
		glutInitWindowSize(width,height);
		glutCreateWindow("2D Lattice Boltzmann");

		glutDisplayFunc(&visualization::display_function);
		glutIdleFunc(&visualization::idle_function);
		glutReshapeFunc(&visualization::reshape_function);
		glutKeyboardFunc(&visualization::keybord_function);
		glutMouseFunc(&visualization::mouse_function);
		glutMotionFunc(&visualization::motion_function);

		//	glEnable(GL_DEPTH_TEST);
		//	glClearColor(1.0,1.0,1.0,1.0);
		//	glEnable(GL_CULL_FACE);

		glewInit();
		
		if (glewIsSupported("GL_VERSION_2_0")) std::cout << "Ready for OpenGL 2.0\n";
		else throw std::runtime_error("OpenGL 2.0 not supported");
		
		// load the shaders
		vel_mag_shader_ptr = std::unique_ptr<shader>( new shader() );
		rho_shader_ptr  = std::unique_ptr<shader>( new shader() );
		u_shader_ptr  = std::unique_ptr<shader>( new shader() );
		v_shader_ptr  = std::unique_ptr<shader>( new shader() );
		wall_shader_ptr = std::unique_ptr<shader>( new shader() );
		
		vel_mag_shader_ptr->load(std::vector<std::string>({"shaders/jet_color_map.glsl", "shaders/velocity_magnitude.vert"}), "shaders/trivial.frag");
		rho_shader_ptr->load(std::vector<std::string>({"shaders/jet_color_map.glsl", "shaders/scalar.vert"}), "shaders/trivial.frag");
		u_shader_ptr->load(std::vector<std::string>({"shaders/jet_color_map.glsl", "shaders/scalar.vert"}), "shaders/trivial.frag");
		v_shader_ptr->load(std::vector<std::string>({"shaders/jet_color_map.glsl", "shaders/scalar.vert"}), "shaders/trivial.frag");
		wall_shader_ptr->load("shaders/wall.vert", "shaders/trivial.frag");

		GLint location;
		
		// create simulation grid geometry
		// -------------------------------
		// vertex buffer object
		float dx = 2.0/(dim_x-1);
		float dy = 2.0/(dim_y-1);
		unsigned int k(0);
		std::vector<float> vertices(simulation_grid_vbo_size);
		for (int j=-buffer_size; j<dim_y+buffer_size; ++j)
		{
			for (int i=-buffer_size; i<dim_x+buffer_size; ++i)
			{
				vertices[k++] = -1.0 + i*dx;
				vertices[k++] = -1.0 + j*dy;
				vertices[k++] = 0.0;
				vertices[k++] = 1.0;
			}
		}
		// make buffer, bind it and fill with values
		simulation_grid_vbo = make_buffer(GL_ARRAY_BUFFER, &vertices[0], simulation_grid_vbo_size*sizeof(float), GL_STATIC_DRAW);
		location = vel_mag_shader_ptr->get_attribute_location("position");
		glVertexAttribPointer(location, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = rho_shader_ptr->get_attribute_location("position");
		glVertexAttribPointer(location, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = u_shader_ptr->get_attribute_location("position");
		glVertexAttribPointer(location, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = v_shader_ptr->get_attribute_location("position");
		glVertexAttribPointer(location, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = wall_shader_ptr->get_attribute_location("position");
		glVertexAttribPointer(location, 4, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		
		// index buffer object
		k=0;
		std::vector<unsigned int> indices(simulation_grid_ibo_size);
		for (int j=0; j<real_dim_y-1; ++j)
		{
			for (int i=0; i<real_dim_x-1; ++i)
			{
				indices[k++] = j*real_dim_x + i + 0;
				indices[k++] = j*real_dim_x + i + 1;
				indices[k++] = (j+1)*real_dim_x + i + 1;
				indices[k++] = (j+1)*real_dim_x + i + 0;
			}
		}
		// make buffer, bind it and fill with values
		simulation_grid_ibo = make_buffer(GL_ELEMENT_ARRAY_BUFFER, &indices[0], simulation_grid_ibo_size*sizeof(unsigned int), GL_STATIC_DRAW);
		
		// create field buffers
		// --------------------
		// make rho buffer, bind it and fill with values
		rho_vbo = make_buffer(GL_ARRAY_BUFFER, rho_data, field_vbo_size*sizeof(float), GL_DYNAMIC_DRAW);
		// make shader attribute
		location = rho_shader_ptr->get_attribute_location("s");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		
		// make u buffer, bind it and fill with values
		u_vbo = make_buffer(GL_ARRAY_BUFFER, u_data, field_vbo_size*sizeof(float), GL_DYNAMIC_DRAW);
		// make shader attribute
		location = vel_mag_shader_ptr->get_attribute_location("u");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = u_shader_ptr->get_attribute_location("s");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		
		// make v buffer, bind it and fill with values
		v_vbo = make_buffer(GL_ARRAY_BUFFER, v_data, field_vbo_size*sizeof(float), GL_DYNAMIC_DRAW);
		// make shader attribute
		location = vel_mag_shader_ptr->get_attribute_location("v");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		location = v_shader_ptr->get_attribute_location("s");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		
		// query walls
		for (unsigned int k=0; k<static_cast<unsigned int>(real_dim_x*real_dim_y); ++k)
		{
			wall[k] = (sim->l.properties.has_flag_property("wall",k) ? 1 : 0);
		}
		// make v buffer, bind it and fill with values
		wall_vbo = make_buffer(GL_ARRAY_BUFFER, wall_data, field_vbo_size*sizeof(float), GL_STATIC_DRAW);
		// make shader attribute
		location = wall_shader_ptr->get_attribute_location("wall");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		
		// set texture states for ibfv
		if (periodic_x) glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT); 
		else glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP); 
		if (periodic_y) glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT); 
		else glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP); 
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
		//glShadeModel(GL_FLAT);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		
		make_patterns();
		
		// setup 2d pixel plotting camera
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0f, width, 0.0f, height, 0.0f, 1.0f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		
		print_usage();
	}
	
public: // static functions

	static void initialize(simulation* sim, int argc, char *argv[])
	{
		if (!vis) vis = std::unique_ptr<visualization>( new visualization(sim,argc,argv) );
	}
	
	static visualization& get_instance()
	{
		if (vis) return *vis;
		else throw std::runtime_error("visualization is not initialized");
	}

	static void display_function() { get_instance().display(); glutPostRedisplay(); }
	
	static void idle_function() { get_instance().idle(); }
	
	static void reshape_function(int w, int h) { get_instance().reshape(w,h); }
	
	static void keybord_function(unsigned char key, int x, int y) 
	{ 
		get_instance().keybord(key,x,y); 
		if(key == 27) // ESC
		{
			vis.release();
			exit(0);
		}
	}
	
	static void motion_function(int x, int y)
	{
		get_instance().motion(x,y);
	}
	
	static void mouse_function(int button, int state, int x, int y)
	{
		get_instance().mouse(button,state,x,y);
	}
	
public: // opengl callback functions

	void run()
	{
		start_time = timer_type::now();
		last_display_time = start_time;
		glutMainLoop();
	}
	
	void reshape(int w, int h)
	{
		glViewport(0, 0, w, h);
		width = w;
		height = h;
		// setup 2d pixel plotting camera
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(0.0f, width, 0.0f, height, 0.0f, 1.0f);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
	}
	
	void display()
	{
		time_point current_time_point = timer_type::now();
		duration dt_d = current_time_point - last_display_time;
		last_display_time = current_time_point;
		current_time = std::chrono::duration_cast<milliseconds>(current_time_point - start_time).count();
		dt = std::chrono::duration_cast<milliseconds>(dt_d).count();
		fps_hist[fps_hist_index] = 1000.0/dt;
		fps_hist_index = (fps_hist_index + 1) % fps_hist.size();
		float fps(0);
		for (unsigned int i=0; i<fps_hist.size(); ++i) fps+=fps_hist[i];
		fps/=fps_hist.size();
		
		//static float dummy_pos = 0.0f;
		static int num_time_steps = 0;
		static const unsigned int num_steps = 15;//25;
		
		if (running)
		{
			for (unsigned int i=0; i<num_steps; ++i) sim->step();
			
			
			if (float_cast)
			{
				for (unsigned int k=0; k<static_cast<unsigned int>(real_dim_x*real_dim_y); ++k)
				{
					rho_data_cast[k] = static_cast<float>(sim->l.rho[k]);
					u_data_cast[k] = static_cast<float>(sim->l.u[k]);
					v_data_cast[k] = static_cast<float>(sim->l.v[k]);
				}
			}
			float rho_min = 999999; float rho_max=-999999;
			float u_min = 999999; float u_max=-999999;
			float v_min = 999999; float v_max=-999999;
			float vel2_min = 999999; float vel2_max=-999999;
			for (unsigned int j=0; j<static_cast<unsigned int>(dim_y); ++j)
			{
				for (unsigned int i=0; i<static_cast<unsigned int>(dim_x); ++i)
				{
					const float tmp = rho_data[(j+buffer_size)*real_dim_x + i + buffer_size];
					rho_min = std::min(rho_min,tmp);
					rho_max = std::max(rho_max,tmp);
				}
			}
			for (unsigned int j=0; j<static_cast<unsigned int>(dim_y); ++j)
			{
				for (unsigned int i=0; i<static_cast<unsigned int>(dim_x); ++i)
				{
					const float tmpu = u_data[(j+buffer_size)*real_dim_x + i + buffer_size];
					const float tmpv = v_data[(j+buffer_size)*real_dim_x + i + buffer_size];
					u_min = std::min(u_min,tmpu);
					u_max = std::max(u_max,tmpu);
					v_min = std::min(v_min,tmpv);
					v_max = std::max(v_max,tmpv);
					vel2_min = std::min(vel2_min,tmpu*tmpu + tmpv*tmpv);
					vel2_max = std::max(vel2_max,tmpu*tmpu + tmpv*tmpv);
				}
			}
			if (num_time_steps == 0)
			{
				rho_min_hist = std::vector<float>(30,rho_min);
				rho_max_hist = std::vector<float>(30,rho_max);
				u_min_hist = std::vector<float>(30,u_min);
				u_max_hist = std::vector<float>(30,u_max);
				v_min_hist = std::vector<float>(30,v_min);
				v_max_hist = std::vector<float>(30,v_max);
				vel_mag2_min_hist = std::vector<float>(30,vel2_min);
				vel_mag2_max_hist = std::vector<float>(30,vel2_max);
			}
			else
			{
				rho_min_hist[data_hist_index] = rho_min;
				rho_max_hist[data_hist_index] = rho_max;
				u_min_hist[data_hist_index] = u_min;
				u_max_hist[data_hist_index] = u_max;
				v_min_hist[data_hist_index] = v_min;
				v_max_hist[data_hist_index] = v_max;
				vel_mag2_min_hist[data_hist_index] = vel2_min;
				vel_mag2_max_hist[data_hist_index] = vel2_max;
			}
			max_rho = 0; min_rho = 0; max_u = 0; min_u = 0; max_v = 0; min_v = 0; max_vel2 = 0; min_vel2 = 0;
			for (unsigned int i=0; i<rho_max_hist.size(); ++i)
			{
				max_rho += rho_max_hist[i];
				min_rho += rho_min_hist[i];
				max_u += u_max_hist[i];
				min_u += u_min_hist[i];
				max_v += v_max_hist[i];
				min_v += v_min_hist[i];
				max_vel2 += vel_mag2_max_hist[i];
				min_vel2 += vel_mag2_min_hist[i];
			}
			max_rho /= rho_max_hist.size();
			min_rho /= rho_max_hist.size();
			max_u /= rho_max_hist.size();
			min_u /= rho_max_hist.size();
			max_v /= rho_max_hist.size();
			min_v /= rho_max_hist.size();
			max_vel2 /= rho_max_hist.size();
			min_vel2 /= rho_max_hist.size();
			data_hist_index = (data_hist_index+1) % rho_max_hist.size();
			
			num_time_steps += num_steps;
		}
		
		// clear screen
		glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		
		glViewport(0, 0, width, height);
		std::stringstream strstr;
		strstr << num_time_steps;
		screen_text(width-panel_size+20, height-230, 1.0, 1.0, 1.0, "step:");
		screen_text(width-panel_size+60, height-230, 1.0, 1.0, 1.0, strstr.str());
		strstr.str("");
		strstr << fps;
		screen_text(width-panel_size+20, height-245, 1.0, 1.0, 1.0, "fps:");
		screen_text(width-panel_size+60, height-245, 1.0, 1.0, 1.0, strstr.str());
		screen_text(width-panel_size+20, height-290, 1.0, 1.0, 1.0, "exit:");
		screen_text(width-panel_size+60, height-290, 1.0, 1.0, 1.0, "<esc>");
		screen_text(width-panel_size+20, height-305, 1.0, 1.0, 1.0, "start:");
		screen_text(width-panel_size+60, height-305, 1.0, 1.0, 1.0, "<space>");
		screen_text(width-panel_size+20, height-320, 1.0, 1.0, 1.0, "mode:");
		screen_text(width-panel_size+60, height-320, 1.0, 1.0, 1.0, "<m>");
		screen_text(width-panel_size+20, height-335, 1.0, 1.0, 1.0, "save:");
		screen_text(width-panel_size+60, height-335, 1.0, 1.0, 1.0, "<s>");
		if (mode >=0 && mode < 5)
		{
			screen_text(width-panel_size+20, height-365, 1.0, 1.0, 1.0, "use left mouse button to");
			screen_text(width-panel_size+20, height-380, 1.0, 1.0, 1.0, "create walls");
			screen_text(width-panel_size+20, height-395, 1.0, 1.0, 1.0, "use <r> to remove walls");
		}
		if (mode == 5)
		{
			screen_text(width-panel_size+20, height-365, 1.0, 1.0, 1.0, "use left mouse button to");
			screen_text(width-panel_size+20, height-380, 1.0, 1.0, 1.0, "create new particles");
		}
		if (mode == 7)
		{
			screen_text(width-panel_size+20, height-365, 1.0, 1.0, 1.0, "use <r> to restart");
		}
		
		glViewport(0, 0, width-panel_size, height);
		
		if (mode == 0) // velocity magnitude
		{
			glBindBuffer(GL_ARRAY_BUFFER, u_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), u_data, GL_DYNAMIC_DRAW);
			GLint location = vel_mag_shader_ptr->get_attribute_location("u");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), v_data, GL_DYNAMIC_DRAW);
			location = vel_mag_shader_ptr->get_attribute_location("v");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			vel_mag_shader_ptr->set_uniform("max_magnitude2",max_vel2);//max_u*max_u*0.2f);
			vel_mag_shader_ptr->set_uniform("min_magnitude2",min_vel2);
			vel_mag_shader_ptr->set_uniform("alpha",1.0f);
			vel_mag_shader_ptr->activate();
			GL_CHECK_AND_PRINT_ERROR
			glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
			GL_CHECK_AND_PRINT_ERROR
			vel_mag_shader_ptr->deactivate();
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"velocity field",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, std::sqrt(min_vel2), std::sqrt(max_vel2), "velocity magnitude");
			
		} else if (mode == 1) // velocity magnitude
		{
			if (running)
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				glTranslatef(-1.0, -1.0, 0.0); 
				glScalef(2.0, 2.0, 1.0);

				glEnable(GL_TEXTURE_2D);
				const float dx = 2.0/(dim_x-1) * ibfv_sample_rate_x;
				const float dy = 2.0/(dim_y-1) * ibfv_sample_rate_y;
				for (int j=-ibfv_additional_pts; j<ibfv_grid_dim_y-ibfv_additional_pts-1; ++j)
				{
					const float y0 = j*dy-0.5;
					const float y1 = (j+1)*dy-0.5;
					glBegin(GL_QUAD_STRIP);
					for (int i=-ibfv_additional_pts; i<ibfv_grid_dim_x-ibfv_additional_pts; ++i)
					{
						const float x = i*dx-0.5;
						unsigned int index0 = ((j*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						unsigned int index1 = (((j+1)*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						const float x0d = x + u_data[index0]*1.0/(dim_x)*num_steps;//*2.0/(dim_x-1)*ibfv_sample_rate_x;
						const float y0d = y0 + v_data[index0]*1.0/(dim_y)*num_steps;//(dim_y-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						const float x1d = x + u_data[index1]*1.0/(dim_x)*num_steps;//*ibfv_sample_rate_x;//*2.0/(dim_x-1);
						const float y1d = y1 + v_data[index1]*1.0/(dim_y)*num_steps;//-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						glTexCoord2f(x, y0); 
						glVertex2f(x0d, y0d);
						glTexCoord2f(x, y1); 
						glVertex2f(x1d, y1d);
					}
					glEnd();
				}
				++frame_index;
				glEnable(GL_BLEND); 
				glCallList(frame_index % (ibfv_num_pattern-1) + 2);
				glBegin(GL_QUAD_STRIP);
					glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
					glTexCoord2f(0.0,  height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(0.0, 1.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), 0.0);  glVertex2f(1.0, 0.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(1.0, 1.0);
				glEnd();
				//glDisable(GL_BLEND);
				//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
				glDisable(GL_TEXTURE_2D);
				glPopMatrix();
			}	
			glBindBuffer(GL_ARRAY_BUFFER, u_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), u_data, GL_DYNAMIC_DRAW);
			GLint location = vel_mag_shader_ptr->get_attribute_location("u");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), v_data, GL_DYNAMIC_DRAW);
			location = vel_mag_shader_ptr->get_attribute_location("v");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			vel_mag_shader_ptr->set_uniform("max_magnitude2",max_vel2);//max_u*max_u*0.2f);
			vel_mag_shader_ptr->set_uniform("min_magnitude2",min_vel2);
			vel_mag_shader_ptr->set_uniform("alpha",0.1f);
			vel_mag_shader_ptr->activate();
			GL_CHECK_AND_PRINT_ERROR
			glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
			GL_CHECK_AND_PRINT_ERROR
			vel_mag_shader_ptr->deactivate();
			
			if (running)
			{
				glDisable(GL_BLEND);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
			}
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"velocity field",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, std::sqrt(min_vel2), std::sqrt(max_vel2), "velocity magnitude");
		}
		else if (mode == 2) // x velocity component
		{
			glBindBuffer(GL_ARRAY_BUFFER, u_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), u_data, GL_DYNAMIC_DRAW);
			GLint location = u_shader_ptr->get_attribute_location("s");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			u_shader_ptr->set_uniform("max_s",max_u);
			u_shader_ptr->set_uniform("min_s",min_u);
			u_shader_ptr->set_uniform("alpha",1.0f);
			u_shader_ptr->activate();
			GL_CHECK_AND_PRINT_ERROR
			glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
			GL_CHECK_AND_PRINT_ERROR
			u_shader_ptr->deactivate();
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"x velocity",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, min_u, max_u, "velocity");
		}
		else if (mode == 3) // y velocity component
		{
			v_shader_ptr->activate();
			glBindBuffer(GL_ARRAY_BUFFER, v_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), v_data, GL_DYNAMIC_DRAW);
			GLint location = v_shader_ptr->get_attribute_location("s");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			v_shader_ptr->set_uniform("max_s",max_v);
			v_shader_ptr->set_uniform("min_s",min_v);
			v_shader_ptr->set_uniform("alpha",1.0f);
			GL_CHECK_AND_PRINT_ERROR
			glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
			GL_CHECK_AND_PRINT_ERROR
			v_shader_ptr->deactivate();
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"y velocity",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, min_v, max_v, "velocity");
		}
		else if (mode == 4) // density
		{
			glBindBuffer(GL_ARRAY_BUFFER, rho_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), rho_data, GL_DYNAMIC_DRAW);
			GLint location = rho_shader_ptr->get_attribute_location("s");
			glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
			glEnableVertexAttribArray(location);
			rho_shader_ptr->set_uniform("max_s",max_rho);
			rho_shader_ptr->set_uniform("min_s",min_rho);
			rho_shader_ptr->set_uniform("alpha",1.0f);
			rho_shader_ptr->activate();
			GL_CHECK_AND_PRINT_ERROR
			glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
			GL_CHECK_AND_PRINT_ERROR
			rho_shader_ptr->deactivate();
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"density",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, min_rho, max_rho, "density");
		}
		else if (mode == 5) // particles
		{
			if (start_particles)
			{
				particle_list.clear();
				start_particles = false;
				const float dx = static_cast<float>(dim_x)/10;
				const float dy = static_cast<float>(dim_y)/10;
				for (int j=0; j<10; ++j)
				{
					for (int i=0; i<10; ++i)
					{
						particle_list.push_back(particle(i*dx+dx/2,j*dy+dy/2));
					}
				}
			}
			//glViewport(0, 0, width, height);
			if (running)
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				glTranslatef(-1.0, -1.0, 0.0); 
				glScalef(2.0, 2.0, 1.0);

				glEnable(GL_TEXTURE_2D);
				const float dx = 2.0/(dim_x-1) * ibfv_sample_rate_x;
				const float dy = 2.0/(dim_y-1) * ibfv_sample_rate_y;
				for (int j=-ibfv_additional_pts; j<ibfv_grid_dim_y-ibfv_additional_pts-1; ++j)
				{
					const float y0 = j*dy-0.5;
					const float y1 = (j+1)*dy-0.5;
					glBegin(GL_QUAD_STRIP);
					for (int i=-ibfv_additional_pts; i<ibfv_grid_dim_x-ibfv_additional_pts; ++i)
					{
						const float x = i*dx-0.5;
						unsigned int index0 = ((j*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						unsigned int index1 = (((j+1)*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						const float x0d = x - u_data[index0]*1.0/(dim_x)*num_steps;//*2.0/(dim_x-1)*ibfv_sample_rate_x;
						const float y0d = y0 - v_data[index0]*1.0/(dim_y)*num_steps;//(dim_y-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						const float x1d = x - u_data[index1]*1.0/(dim_x)*num_steps;//*ibfv_sample_rate_x;//*2.0/(dim_x-1);
						const float y1d = y1 - v_data[index1]*1.0/(dim_y)*num_steps;//-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						glTexCoord2f(x, y0); 
						glVertex2f(x0d, y0d);
						glTexCoord2f(x, y1); 
						glVertex2f(x1d, y1d);
					}
					glEnd();
				}
				glEnable(GL_BLEND); 
				glCallList(1);
				glBegin(GL_QUAD_STRIP);
					glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
					glTexCoord2f(0.0,  height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(0.0, 1.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), 0.0);  glVertex2f(1.0, 0.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(1.0, 1.0);
				glEnd();
				//glDisable(GL_BLEND);
				//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
				glDisable(GL_TEXTURE_2D);
				glPopMatrix();
			}	
			for (auto iter=particle_list.begin(); iter!=particle_list.end(); ++iter)
			{
				const float vx = bilinear_interpolation(u_data,iter->x,iter->y);
				const float vy = bilinear_interpolation(v_data,iter->x,iter->y);
				if (running)
				{
					iter->x += vx*num_steps;
					iter->y += vy*num_steps;
					if (periodic_x && iter->x >= dim_x) iter->x -= dim_x;
					if (periodic_x && iter->x < 0) iter->x += dim_x;
					if (periodic_y && iter->y >= dim_y) iter->y -= dim_y;
					if (periodic_y && iter->y < 0) iter->y += dim_y;
				}
				float red, green, blue;
				get_jet_color(std::min(std::max((vx*vx+vy*vy - min_vel2)/(max_vel2-min_vel2),0.0f),1.0f), red, green, blue);
				//get_jet_color(std::min((vx*vx + vy*vy) / (max_u*max_u*0.2f),1.0f), red, green, blue);
				glColor4f(red,green,blue,0.5f);
				draw_circle((iter->x/dim_x)*(width), (iter->y/dim_y)*height, 3.0, 12);
			}
			if (running)
			{
				glDisable(GL_BLEND);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
			}
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"particles",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, std::sqrt(min_vel2), std::sqrt(max_vel2), "velocity magnitude");
		}
		else if (mode == 6) // "vectors"
		{
			//glViewport(0, 0, width, height);
			//if (use_ibfv && running)
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				glTranslatef(-1.0, -1.0, 0.0); 
				glScalef(2.0, 2.0, 1.0);

				glEnable(GL_TEXTURE_2D);
				const float dx = 2.0/(dim_x-1) * ibfv_sample_rate_x;
				const float dy = 2.0/(dim_y-1) * ibfv_sample_rate_y;
				for (int j=-ibfv_additional_pts; j<ibfv_grid_dim_y-ibfv_additional_pts-1; ++j)
				{
					const float y0 = j*dy-0.5;
					const float y1 = (j+1)*dy-0.5;
					glBegin(GL_QUAD_STRIP);
					for (int i=-ibfv_additional_pts; i<ibfv_grid_dim_x-ibfv_additional_pts; ++i)
					{
						const float x = i*dx-0.5;
						unsigned int index0 = ((j*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						unsigned int index1 = (((j+1)*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						const float x0d = x + u_data[index0]*1.0/(dim_x)*num_steps;//*2.0/(dim_x-1)*ibfv_sample_rate_x;
						const float y0d = y0 + v_data[index0]*1.0/(dim_y)*num_steps;//(dim_y-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						const float x1d = x + u_data[index1]*1.0/(dim_x)*num_steps;//*ibfv_sample_rate_x;//*2.0/(dim_x-1);
						const float y1d = y1 + v_data[index1]*1.0/(dim_y)*num_steps;//-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						glTexCoord2f(x, y0); 
						glVertex2f(x0d, y0d);
						glTexCoord2f(x, y1); 
						glVertex2f(x1d, y1d);
					}
					glEnd();
				}
				glEnable(GL_BLEND); 
				glCallList(1);
				glBegin(GL_QUAD_STRIP);
					glTexCoord2f(0.0,  0.0);  glVertex2f(0.0, 0.0);
					glTexCoord2f(0.0,  height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(0.0, 1.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), 0.0);  glVertex2f(1.0, 0.0);
					glTexCoord2f((width-panel_size)/(ibfv_scale*ibfv_texture_dim), height/(ibfv_scale*ibfv_texture_dim)); glVertex2f(1.0, 1.0);
				glEnd();
				//glDisable(GL_BLEND);
				//glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
				glDisable(GL_TEXTURE_2D);
				glPopMatrix();
			}	
			int nvx = (width-panel_size)/30;
			int nvy = (height)/30;
			const float dx = static_cast<float>(dim_x)/nvx;
			const float dy = static_cast<float>(dim_y)/nvy;
			for (int j=0; j<nvy; ++j)
			{
				for (int i=0; i<nvx; ++i)
				{
					const float x = (i*dx+dx/2);
					const float y = (j*dy+dy/2);
					const float xs = (x/dim_x)*(width);
					const float ys = (y/dim_y)*height;
					const float vx = bilinear_interpolation(u_data,x,y);
					const float vy = bilinear_interpolation(v_data,x,y);
					float red, green, blue;
					get_jet_color(std::min(std::max((vx*vx+vy*vy - min_vel2)/(max_vel2-min_vel2),0.0f),1.0f), red, green, blue);
					//get_jet_color(std::min((vx*vx + vy*vy) / (max_u*max_u*0.2f),1.0f), red, green, blue);
					glColor4f(red,green,blue,std::min(std::max((vx*vx+vy*vy - min_vel2)/(max_vel2-min_vel2),0.0f),1.0f));
					draw_circle(xs,ys, 3.0, 12);
				}
			}
			//if (use_ibfv && running)
			{
				glDisable(GL_BLEND);
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
			}
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"velocity vectors",GLUT_BITMAP_HELVETICA_18);
			draw_color_map_reference(width-panel_size+20, height-200, height-40, std::sqrt(min_vel2), std::sqrt(max_vel2), "velocity magnitude");
		}
		else if (mode == 7) // time line
		{
			if (running && start_timeline)
			{
				start_timeline = false;
				for (int i=40; i<width-40; i+=90)
				{
					float red0(1), green0(0), blue0(0);
					float red1(0), green1(0), blue1(1);
					glBegin(GL_QUAD_STRIP);
					for (int j=0; j<height+30; j+=30)
					{
						glColor4f(red0,green0,blue0,1.0f);
						glVertex2i(i-10,j);
						glVertex2i(i+10,j);
						std::swap(red0,red1);
						std::swap(green0,green1);
						std::swap(blue0,blue1);
					}
					glEnd();
				}
				for (int j=40; j<height-40; j+=90)
				{
					float red0(1), green0(0), blue0(0);
					float red1(0), green1(0), blue1(1);
					glBegin(GL_QUAD_STRIP);
					for (int i=0; i<width+30; i+=30)
					{
						glColor4f(red0,green0,blue0,1.0f);
						glVertex2i(i,j-10);    glVertex2i(i,j+10);
						std::swap(red0,red1);
						std::swap(green0,green1);
						std::swap(blue0,blue1);
					}
					glEnd();
				}
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
			}
			if (running)
			{
				glMatrixMode(GL_PROJECTION);
				glPushMatrix();
				glLoadIdentity();
				glTranslatef(-1.0, -1.0, 0.0); 
				glScalef(2.0, 2.0, 1.0);

				glEnable(GL_TEXTURE_2D);
				const float dx = 2.0/(dim_x-1) * ibfv_sample_rate_x;
				const float dy = 2.0/(dim_y-1) * ibfv_sample_rate_y;
				for (int j=-ibfv_additional_pts; j<ibfv_grid_dim_y-ibfv_additional_pts-1; ++j)
				{
					const float y0 = j*dy-0.5;
					const float y1 = (j+1)*dy-0.5;
					glBegin(GL_QUAD_STRIP);
					for (int i=-ibfv_additional_pts; i<ibfv_grid_dim_x-ibfv_additional_pts; ++i)
					{
						const float x = i*dx-0.5;
						unsigned int index0 = ((j*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						unsigned int index1 = (((j+1)*ibfv_sample_rate_y+real_dim_y) % real_dim_y)*real_dim_x + 
											  ((i*ibfv_sample_rate_x+real_dim_x) % real_dim_x);
						const float x0d = x + u_data[index0]*1.0/(dim_x)*num_steps;//*2.0/(dim_x-1)*ibfv_sample_rate_x;
						const float y0d = y0 + v_data[index0]*1.0/(dim_y)*num_steps;//(dim_y-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						const float x1d = x + u_data[index1]*1.0/(dim_x)*num_steps;//*ibfv_sample_rate_x;//*2.0/(dim_x-1);
						const float y1d = y1 + v_data[index1]*1.0/(dim_y)*num_steps;//-1)*ibfv_sample_rate_y;//*2.0/(dim_y-1);
						glTexCoord2f(x, y0); 
						glVertex2f(x0d, y0d);
						glTexCoord2f(x, y1); 
						glVertex2f(x1d, y1d);
					}
					glEnd();
				}
				glDisable(GL_TEXTURE_2D);
				glPopMatrix();
				glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, width-panel_size, height, 0);
			}
			
			draw_walls();
			
			glViewport(0, 0, width, height);
			screen_text(width-panel_size+20,height-20,1,1,1,"time lines",GLUT_BITMAP_HELVETICA_18);
		}
		
		
		glutSwapBuffers();
	}
	
	void draw_walls()
	{
		glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glBindBuffer(GL_ARRAY_BUFFER, wall_vbo);
		//glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), wall_data, GL_STATIC_DRAW);
		GLint location = wall_shader_ptr->get_attribute_location("wall");
		glVertexAttribPointer(location, 1, GL_FLOAT, GL_FALSE, 0, (GLvoid *) 0);
		glEnableVertexAttribArray(location);
		wall_shader_ptr->activate();
		GL_CHECK_AND_PRINT_ERROR
		glDrawElements(GL_QUADS, simulation_grid_ibo_size, GL_UNSIGNED_INT, 0);
		GL_CHECK_AND_PRINT_ERROR
		wall_shader_ptr->deactivate();
		glDisable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	
	void idle() { }
	
	void keybord(unsigned char key, int x, int y)
	{
		if(key == 32) // space
		{
			if (running) running = false;
			else running = true;
		} 
		else if (key == 109) // m
		{
			mode = (mode + 1) % num_modes;
			if (mode == 5) start_particles = true;
			if (mode == 7) start_timeline = true;
		}
		else if (key == 114) // r
		{
			if (mode >= 0 && mode < 5)
			{
				sim->l.delete_walls();
				// query walls
				for (unsigned int k=0; k<static_cast<unsigned int>(real_dim_x*real_dim_y); ++k)
				{
					wall[k] = (sim->l.properties.has_flag_property("wall",k) ? 1 : 0);
				}
				glBindBuffer(GL_ARRAY_BUFFER, wall_vbo);
				glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), wall_data, GL_STATIC_DRAW);
			}
			if (mode == 7) start_timeline = true;
		}
		else if (key == 115) // s
		{
			std::stringstream strstream;
			strstream << "LB2D_" << std::setw(4) << std::setfill('0') << picture_index++ << ".bmp";
			saveFrameBuffer(strstream.str());
		}
	}
	
	void motion(int x, int y)
	{
		if (wall_drawing_mode && x >=0 && x < width-panel_size && y >=0 && y < height)
		{
			const int i = (static_cast<float>(x)/(width-panel_size))*dim_x;
			const int j = (static_cast<float>(height-y)/height)*dim_y;
			sim->l.add_wall(coordinate<int>(i,j),coordinate<int>(i,j));
			// query walls
			for (unsigned int k=0; k<static_cast<unsigned int>(real_dim_x*real_dim_y); ++k)
			{
				wall[k] = (sim->l.properties.has_flag_property("wall",k) ? 1 : 0);
			}
			glBindBuffer(GL_ARRAY_BUFFER, wall_vbo);
			glBufferData(GL_ARRAY_BUFFER, field_vbo_size*sizeof(float), wall_data, GL_STATIC_DRAW);
		}
	}
	
	void mouse(int button, int state, int x, int y)
	{
		if ( mode >= 0 && mode < 5 && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && x < width-panel_size)
		{
			wall_drawing_mode = true;
		}
		if ( button == GLUT_LEFT_BUTTON && state == GLUT_UP && wall_drawing_mode)
		{
			wall_drawing_mode = false;
		}
		if (mode == 5 && button == GLUT_LEFT_BUTTON && state == GLUT_UP && x < width-panel_size)
		{
			particle_list.push_back(particle((static_cast<float>(x)/(width-panel_size))*dim_x, (static_cast<float>(height-y)/height)*dim_y ));
		}
	}
	
	void draw_color_map_reference(int x_min, int y_min, int y_max, float min_value, float max_value, std::string title) 
	{
		int x_max = x_min+30;
		glBegin(GL_QUAD_STRIP);
		for (int j=y_min+30; j<y_max+1-18; ++j)
		{
			const float tmp = (static_cast<float>(j)-y_min-30) / (y_max+1-18-y_min-30);
			float red,blue,green;
			get_jet_color(tmp,red,green,blue);
			glColor4f(red,green,blue,1.0);
			glVertex2i(x_min, j);
			glVertex2i(x_max, j);
		}
		glEnd();
		screen_text(x_min, y_min, 1.0, 1.0, 1.0, "min:");
		std::stringstream strstr;
		strstr << min_value;
		screen_text(x_max+10, y_min, 1.0, 1.0, 1.0, strstr.str());
		screen_text(x_min, y_min+15, 1.0, 1.0, 1.0, "max:");
		strstr.str("");
		strstr << max_value;
		screen_text(x_max+10, y_min+15, 1.0, 1.0, 1.0, strstr.str());
		screen_text(x_min, y_max-12, 1.0, 1.0, 1.0, title);
	}
	
	void get_jet_color(float value, float& red, float& green, float& blue) 
	{
		const float four_value = 4 * value;
		red   = std::max(std::min(std::min(four_value - 1.5f, -four_value + 4.5f),1.0f),0.0f);
		green = std::max(std::min(std::min(four_value - 0.5f, -four_value + 3.5f),1.0f),0.0f);
		blue  = std::max(std::min(std::min(four_value + 0.5f, -four_value + 2.5f),1.0f),0.0f);
	}
	
	void draw_circle(float cx, float cy, float r, int num_segments) 
	{ 
		glBegin(GL_POLYGON); 
		for(int ii = 0; ii < num_segments; ii++) 
		{ 
			float theta = 2.0f * 3.1415926f * float(ii) / float(num_segments);//get the current angle 

			float x = r * std::cos(theta);//calculate the x component 
			float y = r * std::sin(theta);//calculate the y component 

			glVertex2f(x + cx, y + cy);//output vertex 

		} 
		glEnd(); 
	}
	
	float bilinear_interpolation(float* field, float x, float y)
	{
		//std::cout << x << "," << y << " -> ";
		int xl(x); if (xl<0) --xl;
		int yl(y); if (yl<0) --yl;
		int xh(xl+1);
		int yh(yl+1);
		const float ax = (x-xl);
		const float ay = (y-yl);
		//std::cout << "(" << xl << "," << yl << " - " << xh << "," << yh << "), ax=" << ax << ", ay=" << ay << std::endl;
		if (periodic_x)
		{
			xl = (xl+dim_x) % dim_x;
			xh = (xh+dim_x) % dim_x;
		}
		else 
		{
			if (xl >= dim_x) xl = dim_x-1; if(xl < 0) xl=0;
			if (xh >= dim_x) xh = dim_x-1; if(xh < 0) xh=0;
		}
		if (periodic_y)
		{
			yl = (yl+dim_y) % dim_y;
			yh = (yh+dim_y) % dim_y;
		}
		else 
		{
			if (yl >= dim_y) yl = dim_y-1; if(yl < 0) yl=0;
			if (yh >= dim_y) yh = dim_y-1; if(yh < 0) yh=0;
		}
		const float v0 = (1.0-ax)*field[(yl+buffer_size)*real_dim_x + xl + buffer_size] + ax*field[(yl+buffer_size)*real_dim_x + xh + buffer_size];
		const float v1 = (1.0-ax)*field[(yh+buffer_size)*real_dim_x + xl + buffer_size] + ax*field[(yh+buffer_size)*real_dim_x + xh + buffer_size];
		return (1.0-ay)*v0 + ay*v1;
	}
	
	void screen_text(int x, int y, float r, float g, float b, std::string str, void* font = GLUT_BITMAP_HELVETICA_12)
	{
		/*GLUT_BITMAP_8_BY_13
		GLUT_BITMAP_9_BY_15
		GLUT_BITMAP_TIMES_ROMAN_10
		GLUT_BITMAP_TIMES_ROMAN_24
		GLUT_BITMAP_HELVETICA_10
		GLUT_BITMAP_HELVETICA_12
		GLUT_BITMAP_HELVETICA_18*/
		glColor3f( r, g, b );
		glRasterPos2f(x, y);
		for(unsigned int i=0; i<str.size(); ++i) 
		{
			glutBitmapCharacter(font, str[i]);
		}
	}
	
	void print_usage()
	{
		std::cout << "\n"
		<< "visualization usage\n"
		<< "-------------------\n"
		<< "exit:             <esc>\n"
		<< "start simulation: <space>\n"
		<< "switch mode:      <m>\n"
		<< "save screenshot:  <s>\n";
		std::cout << "\n- in velocty/density field mode use mouse/<r> to create/delete walls\n- in particle mode use mouse to create particles\n- in time line mode use <r> to restart\n";
	}
	
	void saveFrameBuffer(std::string outfile) 
	{
		std::vector<unsigned char> image(width*height*4);
		glPixelStorei(GL_PACK_ALIGNMENT, 1);
		glReadPixels(0,0,width,height, GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, &image[0]);
		BMP bmp;
		bmp.SetSize(width, height);
		bmp.SetBitDepth(32);
		for(int j=0; j<bmp.TellHeight(); ++j) 
		{
			bmp.Read32bitRow(&image[j*width*4], width*4, bmp.TellHeight()-j-1);
		}
		bmp.WriteToFile(outfile.c_str());
	}
	
private: // ibfv

	void make_patterns()
	{
		int lut[256];
		int phase[ibfv_texture_dim][ibfv_texture_dim];
		GLubyte pat[ibfv_texture_dim][ibfv_texture_dim][4];
		int i, j, k, t;

		for (i = 0; i < 256; i++) lut[i] = i < 127 ? 0 : 255;
		for (i = 0; i < ibfv_texture_dim; i++)
			for (j = 0; j < ibfv_texture_dim; j++) phase[i][j] = std::rand() % 256; 

		for (k = 0; k < ibfv_num_pattern; k++) 
		{
			t = k*256/ibfv_num_pattern;
			for (i = 0; i < ibfv_texture_dim; i++) 
				for (j = 0; j < ibfv_texture_dim; j++) 
				{
					if (k==0)
					{
						pat[i][j][0] =
						pat[i][j][1] =
						pat[i][j][2] = 0;
						pat[i][j][3] = ibfv_alpha;
					}
					else
					{
						pat[i][j][0] =
						pat[i][j][1] =
						pat[i][j][2] = lut[(t + phase[i][j]) % 255];
						pat[i][j][3] = ibfv_alpha;
					}
				}
			glNewList(k + 1, GL_COMPILE);
				glTexImage2D(GL_TEXTURE_2D, 0, 4, ibfv_texture_dim, ibfv_texture_dim, 0, GL_RGBA, GL_UNSIGNED_BYTE, pat);
			glEndList();
		}
	}
	
private: // impl

	#pragma pack(push, 1)
	// A vertex containing only a position value
	struct vertex
	{
		vertex() {}
		vertex(float x, float y, float z) : x(x), y(y), z(z) {}
		float x, y, z;
	};
	#pragma pack(pop)

	struct particle
	{
		particle(){}
		particle(float _x, float _y) : x(_x), y(_y) {}
		float x;
		float y;
	};

private: // members

	static std::unique_ptr<visualization> vis;
	// window properties
	const int panel_size;
	int width, height;
	// pointer to simulation
	simulation* sim;
	// simulation properties
	const int dim_x;
	const int dim_y;
	const int buffer_size;
	const int real_dim_x;
	const int real_dim_y;
	const bool periodic_x;
	const bool periodic_y;
	// if simulation does not use floats -> cast every time
	const bool float_cast;
	std::vector<float> rho_data_cast;
	std::vector<float> u_data_cast;
	std::vector<float> v_data_cast;
	// wall container
	std::vector<float> wall;
	// data pointers
	float* rho_data;
	float* u_data;
	float* v_data;
	float* wall_data;
	// opacity value for drawing
	float alpha0;
	// ibfv params
	const int ibfv_alpha;
	const int ibfv_texture_dim;
	const int ibfv_sample_rate_x;
	const int ibfv_sample_rate_y;
	const int ibfv_additional_pts;
	const int ibfv_grid_dim_x;
	const int ibfv_grid_dim_y;
	const float ibfv_scale;
	const int ibfv_num_pattern;
	int frame_index;
	std::vector<float> ibfv_u_data;
	std::vector<float> ibfv_v_data;
	std::vector<unsigned int> ibfv_x_index_lookup;
	std::vector<unsigned int> ibfv_y_index_lookup;
	// opengl buffers (location and sizes)
	GLuint simulation_grid_vbo;
	int simulation_grid_vbo_size;
	GLuint simulation_grid_ibo;
	int simulation_grid_ibo_size;
	GLuint rho_vbo;
	GLuint u_vbo;
	GLuint v_vbo;
	GLuint wall_vbo;
	int field_vbo_size;
	// holder for shader objects
	std::unique_ptr<shader> vel_mag_shader_ptr;
	std::unique_ptr<shader> rho_shader_ptr;
	std::unique_ptr<shader> u_shader_ptr;
	std::unique_ptr<shader> v_shader_ptr;
	std::unique_ptr<shader> wall_shader_ptr;
	// timings
	time_point start_time;
	time_point last_display_time;
	float current_time;
	float dt;
	std::vector<float> fps_hist;
	unsigned int fps_hist_index;
	// data extent
	std::vector<float> rho_max_hist;
	std::vector<float> rho_min_hist;
	std::vector<float> u_max_hist;
	std::vector<float> u_min_hist;
	std::vector<float> v_max_hist;
	std::vector<float> v_min_hist;
	std::vector<float> vel_mag2_max_hist;
	std::vector<float> vel_mag2_min_hist;
	unsigned int data_hist_index;
	// limits for simulation data
	float min_rho, max_rho, min_u, max_u, min_v, max_v, min_vel2, max_vel2;
	//particles
	std::list<particle> particle_list;
	// sate flags
	bool running;
	const int num_modes;
	int mode;
	bool start_particles;
	bool start_timeline;
	unsigned int picture_index;
	bool wall_drawing_mode;
};

std::unique_ptr<visualization> visualization::vis = nullptr;

} // lb


#endif // LB_VISUALIZATION_HPP_INCLUDED
