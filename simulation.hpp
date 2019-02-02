/**
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 *  this is the right folder
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

	#include "H_root.hpp"
	#include "lattice.hpp"
	#include <sstream>


	namespace lb {

		struct boundary_info
	    {
					int idx;
	        node n;
					int i;
					int j;
	        std::array<float_type,9> distance = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
					std::array<float_type,9> population_prev = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
					bool IS_BOUNDARY = 0;
	    };

		struct wall_info{
				int idx;
		    node n;
				int i;
				int j;
				bool IS_WALL = 0;
		};

		struct wall_info_basic{
				node n;
		};
		/**
		 *  @brief Simulation class implementing LB
		 *
		 *  This class holds a lattice as member (see @ref simulation::l) and
		 *  carries out the simulation steps on top of it. The main methods of
		 *  this class are @ref simulation::advect() and
		 *  @ref simulation::collide().
		 */

		class simulation
		{
		public: // ctor

		/**
		 *  @brief Construct from domain size and flow parameters
		 *  @param[in] nx    extent in x direction
		 *  @param[in] ny    extent in y direction
		 *  @param[in] _Re   Reynolds number
		 *  @param[in] _Vmax mean flow velocity
		 */
		simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax, unsigned int KBC_COLOR_, unsigned int BC_topbottom,bool WRITE_FILE_, float_type _Slope)
		: l(nx, ny),
		  shift(velocity_set().size),
		  Re(_Re),
			kbc_coluring_scheme(KBC_COLOR_),
			TopBotBC(BC_topbottom),
		  Vmax(_Vmax),
	    R(l.ny/10),
			x_c(l.nx/4),
			y_c(l.ny/2),
			Lx(x_c + l.nx/10),
	    visc(Vmax*R/Re),
	    beta(1/(((2*visc)/pow(1/std::sqrt(3.0),2))+1.0)),
			cs(1/std::sqrt(3.0)),
			slope(_Slope),
			Force_X(0),
			Force_Y(0),
		  time(0),
		  file_output(WRITE_FILE_), // set to true if you want to write files
		  output_freq(200),
		  output_index(0),
			wall_dummy(0)
		{
			// define amount to shift populations for advection
			for (unsigned int i=0; i<velocity_set().size; ++i)   //size is the number of velocities
			{
				shift[i] = velocity_set().c[0][i]+velocity_set().c[1][i]*l.real_nx;
			}
		}

		/**
		 *  @brief Initialize the flow field
		 *  Initialization includes defining initial density, velocity and
		 *  populations. You can use Taylor-Green vortex flow conditions.
		 */
		void initialize(){
			// **************************
			// * fill in your code here *
			// * the lines below are    *
			// * just examples          *
			// **************************
					// INITIAL VALUE OF RHO?
	        float_type rho=1.0;

			for (int j=0; j<static_cast<int>(l.ny); ++j){
				for (int i=0; i<static_cast<int>(l.nx); ++i){
	        l.get_node(i,j).u()   = Vmax;
					l.get_node(i,j).v()   = 0;
          l.get_node(i,j).rho() =	rho;
					l.get_node(i,j).gamma() = 2. ;
          velocity_set().equilibrate(l.get_node(i,j));
				}
			}

					//set_SimpleNoseconeBoundary();
					//set_BoattailNoseconeBoundary();
					//set_FinNoseconeBoundary();
	        std::cout << "found " << boundary_nodes.size() << " boundary nodes" << std::endl;
		return ;
		}

			void set_SimpleNoseconeBoundary(){
				float_type tol = 0;
		    float_type a, b, c, dx, D;
		    for (unsigned int y=0; y<l.ny+1; ++y){            //direction y , y is y_0
		        for (unsigned int x=0; x<l.nx+1; ++x){        //direction x , x is x_0
								if ( x >= (x_c - 10) && x <= (Lx+10)  &&  (y >= y_c - R - 10)  && y <= (y_c + R + 10) ){
										wall_info wi;
										if ((float_type)x <= Lx) { 										// DO THE REGULAR THING
												if ((y-y_c)*(y-y_c) - (R*R)*(x-x_c)/(Lx - x_c) < tol  ){
													  wi.IS_WALL = true;
														wi.n = l.get_node(x,y);
														wi.i = x;
														wi.j = y;
														wall_nodes.push_back(wi);
												} else {
													wi.IS_WALL = false;
													wi.n = l.get_node(x,y);
													wi.i = x;
													wi.j = y;
													wall_nodes.push_back(wi);
												}
										} else {
											wi.IS_WALL = false;
											wi.n = l.get_node(x,y);
											wi.i = x;
											wi.j = y;
											wall_nodes.push_back(wi);
										}
								}
						}
				}

				for (auto wi : wall_nodes){
						int x = wi.i; int y = wi.j; //node n = wi.n;
						if (!wi.IS_WALL){
								boundary_info bi;
								bool found_intersection = false;
								if (x <= Lx){
										for (int i=1; i<9; ++i){        //9 directions
												a = velocity_set().c[1][i]*velocity_set().c[1][i];
												b = 2*velocity_set().c[1][i]*(y - y_c) - R*R*velocity_set().c[0][i]/(Lx - x_c);
												c =  (y-y_c)*(y-y_c) - R*R*(x-x_c)/(Lx - x_c);
												D = b*b-4*a*c;

												if (D>=0){
													const auto dx1 =(-b + std::sqrt(b*b-4*a*c))/(2*a);
													const auto dx2 =(-b - std::sqrt(b*b-4*a*c))/(2*a);
													dx = 0;
													if (dx1 > 0) dx = dx1;
													if (dx2 > 0 && dx2 < dx1) dx = dx2;

													if (dx>0 && dx<=1){
														bi.distance[i] = dx;
														found_intersection = true;
													}
										}
									}
								} else {
										if ((y <= y_c + R) && (y >= y_c - R) ){
												for (int i=1; i<9; ++i){
														dx = (Lx-x)/(velocity_set().c[0][i]);
														if (dx>0 && dx<=1){
																bi.distance[i] = dx;
																found_intersection = true;
														}
												}
										}
								}
								if (found_intersection){
										bi.n = l.get_node(x,y);
										bi.i = x;
										bi.j = y;
										bi.IS_BOUNDARY = true;
										boundary_nodes.push_back(bi);
								}
						}
				}

		}

		void set_BoattailNoseconeBoundary(){
			float_type tol = 0;
			float_type a, b, c, dx, D;
			float_type m = slope; // slope of the Boattail Section
			float_type H = R*(1-m);
			for (unsigned int y=0; y<l.ny+1; ++y){            //direction y , y is y_0
					for (unsigned int x=0; x<l.nx+1; ++x){        //direction x , x is x_0
							if ( x >= (x_c - 10) && x <= (Lx+R+5)  &&  (y >= y_c - R - 10)  && y <= (y_c + R + 10) ){
									wall_info wi;
									if ((float_type)x <= Lx) { 										// DO THE REGULAR THING
											if ((y-y_c)*(y-y_c) - (R*R)*(x-x_c)/(Lx - x_c) < tol  ){
													wi.IS_WALL = true;
													wi.n = l.get_node(x,y);
													wi.i = x;
													wi.j = y;
													wall_nodes.push_back(wi);
											} else {
													wi.IS_WALL = false;
													wi.n = l.get_node(x,y);
													wi.i = x;
													wi.j = y;
													wall_nodes.push_back(wi);
											}
									} else if ((float_type)x > Lx && (float_type)x <= Lx + R){
											if ((y < y_c - m*(x-Lx) + R) && (y > y_c + m*(x-Lx) - R) ){
													wi.IS_WALL = true;
													wi.n = l.get_node(x,y);
													wi.i = x;
													wi.j = y;
													wall_nodes.push_back(wi);
											} else {
													wi.IS_WALL = false;
													wi.n = l.get_node(x,y);
													wi.i = x;
													wi.j = y;
													wall_nodes.push_back(wi);
											}
									} else if ( (float_type)x > Lx + R ){
											wi.IS_WALL = false;
											wi.n = l.get_node(x,y);
											wi.i = x;
											wi.j = y;
											wall_nodes.push_back(wi);
									}
							}
					}
			}
			for (auto wi : wall_nodes){
					int x = wi.i; int y = wi.j; //node n = wi.n;
					if (!wi.IS_WALL){
							boundary_info bi;
							bool found_intersection = false;
							if ((float_type)x <= Lx){
									for (int i=1; i<9; ++i){        //9 directions
											a = velocity_set().c[1][i]*velocity_set().c[1][i];
											b = 2*velocity_set().c[1][i]*(y - y_c) - R*R*velocity_set().c[0][i]/(Lx - x_c);
											c =  (y-y_c)*(y-y_c) - R*R*(x-x_c)/(Lx - x_c);
											D = b*b-4*a*c;
											if (D>=0){
													const auto dx1 =(-b + std::sqrt(b*b-4*a*c))/(2*a);
													const auto dx2 =(-b - std::sqrt(b*b-4*a*c))/(2*a);
													dx = 0;
													if (dx1 > 0) dx = dx1;
													if (dx2 > 0 && dx2 < dx1) dx = dx2;

													if (dx>0 && dx<=1){
														bi.distance[i] = dx;
														found_intersection = true;
													}
											}
						   		}
							} else if ((float_type)x > Lx && (float_type)x <= Lx + R){
									for (int i=1; i<9; ++i){
											float_type dx = 0;
											if ((y >= y_c - m*(x-Lx) + R)) {
													dx = (-m*(x-Lx)+R-(y-y_c))/(velocity_set().c[1][i] + m*velocity_set().c[0][i]);
											} else if ((y <= y_c + m*(x-Lx) - R) ){
													dx = ( m*(x-Lx)-R-(y-y_c))/(velocity_set().c[1][i] - m*velocity_set().c[0][i]);
											} else {std::cerr<<"ERROR CALCULATING DISTANCE TO WALL IN BOATTAL SECTION"<<'\n';}
											if (dx>0 && dx<=1){
													bi.distance[i] = dx;
													found_intersection = true;
											}
									}
							} else if ( (float_type)x >= Lx + R ){
								// after Boattail | wall | change y limits -> (yc+H,yc-H)
								if ((y <= y_c + H ) && (y >= y_c-H) ){
										for (int i=1; i<9; ++i){
												dx = (Lx+R-x)/(velocity_set().c[0][i]);
												if (dx>0 && dx<=1){
														bi.distance[i] = dx;
														found_intersection = true;
												}
										}
								}
							}
							if (found_intersection){
									bi.n = l.get_node(x,y);
									bi.i = x;
									bi.j = y;
									bi.IS_BOUNDARY = true;
									boundary_nodes.push_back(bi);
							}
					}
			}
		}

		/*
		void set_FinNoseconeBoundary(){

			// coefficients of a.dx^2 + b.dx + c = 0
			float_type a, b, c, dx, D;
			float_type H = R/5;
			float_type L_fin = 2*R;
			for (unsigned int y=0; y<l.ny+1; ++y){            //direction y , y is y_0
					for (unsigned int x=0; x<l.nx+1; ++x){        //direction x , x is x_0
							if ( x >= (x_c - 5) && x <= (Lx+L_fin+5)  &&  (y >= y_c - R - 5)  && y <= (y_c + R + 5) ){
									boundary_info bi;
									bool found_intersection = false;

									if ((float_type)x <= Lx) {
										// DO THE REGULAR THING
										if ((y-y_c)*(y-y_c) - (R*R)*(x-x_c)/(Lx - x_c) > 0  ){ //check if node is outside nosecone, search for nearest boundary node
											for (int i=1; i<9; ++i){        //9 directions
												a = velocity_set().c[1][i]*velocity_set().c[1][i];
												b = 2*velocity_set().c[1][i]*(y - y_c) - R*R*velocity_set().c[0][i]/(Lx - x_c);
												c =  (y-y_c)*(y-y_c) - R*R*(x-x_c)/(Lx - x_c);
												D = b*b-4*a*c;

												if (D>=0){
													const auto dx1 =(-b + std::sqrt(b*b-4*a*c))/(2*a);
													const auto dx2 =(-b - std::sqrt(b*b-4*a*c))/(2*a);
													dx = 0;
													if (dx1 > 0) dx = dx1;
													if (dx2 > 0 && dx2 < dx1) dx = dx2;

													if (dx>0 && dx<=1){
														bi.distance[i] = dx;
														found_intersection = true;
													}
												}
											}
										}		else	{
											wall_info_basic wi{l.get_node(x,y)};
											wall_nodes.push_back(wi);
										}	 		// CHANGES HERE ONWARDS
									}	else if ( (float_type)x >= Lx && (float_type)x <= Lx + L_fin) {
											// Boattail section y = +-m*x + R
											if ((y > y_c + H)  ){
													for (int i=1; i<9; ++i){
														float_type dx=0;
														float_type dx_H=0;
														float_type dx_V=0;
															dx_H = (H - y + y_c)/(velocity_set().c[1][i]);
															if (y <= y_c + R) dx_V = (Lx - x)/(velocity_set().c[0][i]);

															if (dx_H > 0) dx = dx_H;
															if (dx_V > 0 && dx_V < dx_H) dx = dx_V;

															if (dx>0 && dx<=1){
																	bi.distance[i] = dx;
																	found_intersection = true;
															}
													}
												} else if ((y < y_c - H) ){
														for (int i=1; i<9; ++i){
															float_type dx=0;
															float_type dx_H=0;
															float_type dx_V=0;
															dx_H = (-H - y + y_c)/(velocity_set().c[1][i]);
															if (y >= y_c - R) dx_V = (Lx - x)/(velocity_set().c[0][i]);

															if (dx_H > 0) dx = dx_H;
															if (dx_V > 0 && dx_V < dx_H) dx = dx_V;

															if (dx>0 && dx<=1){
																	bi.distance[i] = dx;
																	found_intersection = true;
															}
														}
												}	else	{
													wall_info_basic wi{l.get_node(x,y)};
													wall_nodes.push_back(wi);
												}
									} else if ( (float_type)x >= Lx + L_fin ){
											// after Boattail | wall | change y limits -> (yc+H,yc-H)
											if ((y < y_c + H ) && (y > y_c - H) ){
													for (int i=1; i<9; ++i){
															dx = (Lx+L_fin-x)/(velocity_set().c[0][i]);
															if (dx>0 && dx<=1){
																	bi.distance[i] = dx;
																	found_intersection = true;
															}
													}
											}
									}
									if (found_intersection){
											bi.n = l.get_node(x,y);
											bi.i = x;
											bi.j = y;
											boundary_nodes.push_back(bi);
									}
							}
					}
			}

		}
		*/
		/**
		 *  @brief advect the populations
		 *
		 *  Include periodic boundary conditions here also
		 */
		void advect(){
	        //first initialize buffers

				// horizontal faces
			for (int ix=0;ix< static_cast<int>(l.nx);++ix){
				for (unsigned int m=0; m<velocity_set().size; ++m){
					l.f[m][l.index(ix,-1)] = l.f[m][l.index(ix,0)];
					l.f[m][l.index(ix,l.ny)] = l.f[m][l.index(ix,l.ny-1)];
				}
			}
					// vertical faces faces
			for (int iy=0;iy < static_cast<int>(l.ny);++iy){
				for (unsigned int m=0; m < velocity_set().size; ++m){
					l.f[m][l.index(-1,iy)] = l.f[m][l.index(l.nx-1,iy)];
					l.f[m][l.index(l.nx,iy)] = l.f[m][l.index(0,iy)];
				}
			}
				/*	// corners
			for (unsigned int m=0; m < velocity_set().size; ++m){
				l.f[m][l.index(-1,-1)] = l.f[m][l.index(l.nx-1,l.ny-1)];
				l.f[m][l.index(l.nx,l.ny)] = l.f[m][l.index(0,0)];
				l.f[m][l.index(l.nx,-1)] = l.f[m][l.index(0,l.ny-1)];
				l.f[m][l.index(-1,l.ny)] = l.f[m][l.index(l.nx-1,0)];
			}*/
	        //shift populations in the opposite directions to the velocities so that not the same population is copied across the whole lattice
			for (unsigned int m=0; m<velocity_set().size; ++m){
				if (shift[m] > 0){
					for (unsigned int n = l.index(l.nx-1,l.ny-1); n>=l.index(0,0);--n){
						l.f[m][n] = l.f[m][n-shift[m]];
					}
				}
				else if (shift[m] < 0){
					for (unsigned int n = l.index(0,0);n<=l.index(l.nx-1,l.ny-1);++n){
						l.f[m][n] = l.f[m][n-shift[m]];
					}
				}
			}

		return ;
		}

		/**  @brief apply wall boundary conditions */
		void wall_bc(){
					// bounce back bc by looping over all boundary cells
					/*
		  for (auto& bi : boundary_nodes){
		      for (int i=1;i<9;++i){
			        if (bi.distance[i] > 0)
			          bi.n.f(velocity_set().rflct_latticeVelocity[i]) = bi.n.f(i);
		      }
		  } 	*/

					// GRAD BOUNDARY CONDITIONS
			for (auto& bi : boundary_nodes){
					int count = 0;
					float_type rho = 0.0;
					float_type u = 0.0;
					float_type v = 0.0;
					const int i = bi.i;
					const int j = bi.j;

					for (int m=1;m<9;++m){
							if (bi.distance[m] > 0){
								bi.n.f(velocity_set().rflct_latticeVelocity[m]) = bi.n.f(m);
								int i_interp,j_interp;
								velocity_set().interpolation_node(i,j,m,i_interp,j_interp);
								auto n_interp = l.get_node(i_interp,j_interp);
										// get velocity via interpolation and get density after bounce back
								u +=  bi.distance[m]*n_interp.u()/( 1 + bi.distance[m] );
								v +=  bi.distance[m]*n_interp.v()/( 1 + bi.distance[m] );
								count++;
							}
							rho += bi.n.f(m);
					}
					bi.n.u() = u/count;
					bi.n.v() = v/count;
					/*
					float_type f_tgt[9];
					velocity_set().f_eq(f_tgt,rho,u,v);
					 rho = 0.0; u = 0.0; v = 0.0;
					for (int m=1;m<9;++m){
							if (bi.distance[m] > 0){
								bi.n.f(velocity_set().rflct_latticeVelocity[m]) = f_tgt[velocity_set().rflct_latticeVelocity[m]];
							}
							rho += bi.n.f(m);
							u += bi.n.f(m)*velocity_set().c[0][m];
							v += bi.n.f(m)*velocity_set().c[1][m];
					}
					u = u/(rho);
					v = v/(rho);

					float_type f_loc[9];
					velocity_set().f_eq(f_loc,rho,u,v);
					rho = 0.0; u = 0.0; v = 0.0;
					for (int m=1;m<9;++m){
							if (bi.distance[m] > 0){
									bi.n.f(velocity_set().rflct_latticeVelocity[m]) += f_tgt[velocity_set().rflct_latticeVelocity[m]]  - f_loc[velocity_set().rflct_latticeVelocity[m]] ;
							}
					}
				} */
					float_type dudx,dudy,dvdx,dvdy;

					if ( bi.distance[5]>0 ){
							dudx = (- l.get_node(i-1,j).u() )/(1+bi.distance[5]);
							dvdx = (- l.get_node(i-1,j).v() )/(1+bi.distance[5]);
							if ( l.get_node(i-1,j).u() == wall_dummy || l.get_node(i-1,j).v() == wall_dummy)
									printf("error proper dist 5\n" );
							dudy = (- l.get_node(i,j-1).u() )/(1+bi.distance[5]);
							dvdy = (- l.get_node(i,j-1).v() )/(1+bi.distance[5]);
							if (  l.get_node(i,j-1).u() == wall_dummy || l.get_node(i,j-1).v() == wall_dummy)
									printf("error proper dist 5\n" );
					} else if (bi.distance[7]>0) {
							dudx = ( l.get_node(i,j).u() )/(1+bi.distance[7]);
							dvdx = ( l.get_node(i,j).v() )/(1+bi.distance[7]);
							if (  l.get_node(i,j).u() == wall_dummy || l.get_node(i,j).v() == wall_dummy)
									printf("error proper dist 7a\n" );
							dudy = (  l.get_node(i,j+1).u() )/(1+bi.distance[7]);
							dvdy = (  l.get_node(i,j+1).v() )/(1+bi.distance[7]);
							if (  l.get_node(i,j+1).u()  == wall_dummy || l.get_node(i,j+1).v() == wall_dummy)
									printf("error proper dist 7b\n" );
					}

					if ( bi.distance[6]>0 ) {
							dudx = ( l.get_node(i,j).u() )/(1+bi.distance[6]);
							dvdx = ( l.get_node(i,j).v() )/(1+bi.distance[6]);
							if (  l.get_node(i,j).u() == wall_dummy || l.get_node(i,j).v() == wall_dummy)
								printf("error proper dist 6a\n" );
							dudy = (- l.get_node(i,j-1).u() )/(1+bi.distance[6]);
							dvdy = (- l.get_node(i,j-1).v() )/(1+bi.distance[6]);
							if (  l.get_node(i,j-1).u() == wall_dummy || l.get_node(i,j-1).v() == wall_dummy)
								printf("error proper dist 6b\n" );
					} else if (bi.distance[8]>0) {
							dudx = (- l.get_node(i-1,j).u() )/(1+bi.distance[8]);
							dvdx = (- l.get_node(i-1,j).v() )/(1+bi.distance[8]);
							if ( l.get_node(i-1,j).u() == wall_dummy || l.get_node(i-1,j).v() == wall_dummy)
								printf("error proper dist 8\n" );
							dudy = (  l.get_node(i,j+1).u() )/(1+bi.distance[8]);
							dvdy = (  l.get_node(i,j+1).v() )/(1+bi.distance[8]);
							if (  l.get_node(i,j+1).u()  == wall_dummy || l.get_node(i,j+1).v() == wall_dummy)
									printf("error proper dist 8\n" );
					}

					if (bi.distance[1] > 0){
						dudx = (- l.get_node(i-1,j).u() )/(1+bi.distance[1]);
						dvdx = (- l.get_node(i-1,j).v() )/(1+bi.distance[1]);
						if ( l.get_node(i-1,j).u() == wall_dummy || l.get_node(i-1,j).v() == wall_dummy)
							printf("error proper dist 1\n" );
					} else if (bi.distance[3] > 0){
						dudx = ( l.get_node(i+1,j).u() )/(1+bi.distance[3]);
						dvdx = ( l.get_node(i+1,j).v() )/(1+bi.distance[3]);
						if (  l.get_node(i+1,j).u() == wall_dummy || l.get_node(i+1,j).v() == wall_dummy)
							printf("error proper dist 3\n" );
					}
					if (bi.distance[2] > 0){
						dudy = (- l.get_node(i,j-1).u() )/(1+bi.distance[2]);
						dvdy = (- l.get_node(i,j-1).v() )/(1+bi.distance[2]);
						if (  l.get_node(i,j-1).u() == wall_dummy || l.get_node(i,j-1).v() == wall_dummy)
							printf("error proper dist 2\n" );
					}  else if (bi.distance[4] > 0){
						dudy = (  l.get_node(i,j+1).u() )/(1+bi.distance[4]);
						dvdy = (  l.get_node(i,j+1).v() )/(1+bi.distance[4]);
						if (  l.get_node(i,j+1).u()  == wall_dummy || l.get_node(i,j+1).v() == wall_dummy)
							printf("error proper dist 4\n" );
					}
					/*
					int ip,ip2,jp,jp2,im,im2,jm,jm2;

						ip = i + 1 ; ip2 = i + 2;
						im = i - 1 ; im2 = i - 2;
						jp = j + 1; jp2 = j + 2;
						jm = j - 1; jm2 = j - 2;
						// One Sided Second Order Approx of dUdX //

					if (bi.distance[1] > 0){
						//backwards
						dudx = ( 3*l.get_node(i,j).u() - 4*l.get_node(im,j).u() + l.get_node(im2,j).u() )/2;
						dvdx = ( 3*l.get_node(i,j).v() - 4*l.get_node(im,j).v() + l.get_node(im2,j).v() )/2;
					} else if (bi.distance[3] > 0){
						//forwards
						dudx = (-3*l.get_node(i,j).u() + 4*l.get_node(ip,j).u() - l.get_node(ip2,j).u() )/2;
						dvdx = (-3*l.get_node(i,j).v() + 4*l.get_node(ip,j).v() - l.get_node(ip2,j).v() )/2;
					} else {
						dudx = ( 3*l.get_node(i,j).u() - 4*l.get_node(im,j).u() + l.get_node(im2,j).u() )/2;
						dvdx = ( 3*l.get_node(i,j).v() - 4*l.get_node(im,j).v() + l.get_node(im2,j).v() )/2;
					}

					if (bi.distance[2] > 0){
						//backwards
						dudy = (+3*l.get_node(i,j).u() - 4*l.get_node(i,jm).u() + l.get_node(i,jm2).u() )/2;
						dvdy = (+3*l.get_node(i,j).v() - 4*l.get_node(i,jm).v() + l.get_node(i,jm2).v() )/2;
					} else if (bi.distance[4] > 0) {
						//forward
						dudy = (-3*l.get_node(i,j).u() + 4*l.get_node(i,jp).u() - l.get_node(i,jp2).u() )/2;
						dvdy = (-3*l.get_node(i,j).v() + 4*l.get_node(i,jp).v() - l.get_node(i,jp2).v() )/2;
					} else {
						dudy = (+3*l.get_node(i,j).u() - 4*l.get_node(i,jm).u() + l.get_node(i,jm2).u() )/2;
						dvdy = (+3*l.get_node(i,j).v() - 4*l.get_node(i,jm).v() + l.get_node(i,jm2).v() )/2;
					}
					*/
					float_type Pxx,Pyy,Pxy;
					Pxx = bi.n.rho()*( cs*cs + bi.n.u()*bi.n.u() - (cs*cs*(2*dudx)/(2*beta)) );
					Pyy = bi.n.rho()*( cs*cs + bi.n.v()*bi.n.v() - (cs*cs*(2*dvdy)/(2*beta)) );
					Pxy = bi.n.rho()*( bi.n.u()*bi.n.v() - (cs*cs*(dudy + dvdx)/(2*beta)) );

					for (int m=1;m<9;++m){
							// only for nodes with missing data i.e if bi.distance[m]>0 ==>
							// pop. to be replaced is bi.n.f(velocity_set().rflct_latticeVelocity[m])
							if (bi.distance[m]>0){
								bi.n.f(velocity_set().rflct_latticeVelocity[m]) = velocity_set().W[m]*(
								bi.n.rho()*( 1+(bi.n.u()*velocity_set().c[0][m]+bi.n.v()*velocity_set().c[1][m])/(cs*cs) ) +
								( ( (Pxx-bi.n.rho()*cs*cs)*(velocity_set().c[0][m]*velocity_set().c[0][m]-cs*cs) ) +
								  ( (Pyy-bi.n.rho()*cs*cs)*(velocity_set().c[1][m]*velocity_set().c[1][m]-cs*cs) ) +
								  (2*(Pxy*velocity_set().c[0][m]*velocity_set().c[1][m]) ) )/(2*std::pow(cs,4))	);
							}
					}
			}

			switch (TopBotBC) { 	// user input
				case 0: // no slip
					for (unsigned int i=1 ; i < l.nx-1 ; ++i){
							auto n_topwall = l.get_node(i,l.ny-1);
							n_topwall.u()   = 0;
							n_topwall.v()   = 0;
							n_topwall.rho() =	1;
							velocity_set().equilibrate(n_topwall);

							auto n_botwall = l.get_node(i,0);
							n_botwall.u()   = 0;
							n_botwall.v()   = 0;
							n_botwall.rho() =	1;
							velocity_set().equilibrate(n_botwall);
					}					// ADD BC HERE WHEN DONE
					break;

				case 1: // slip // FIX THIS IF NEEDED ACC TO FABIAN
					for (unsigned int i=1 ; i < l.nx-1 ; ++i){
						l.get_node(l.index(i,l.ny-1)).f(8)=l.get_node(l.index(i-1,l.ny-2)).f(5);
						l.get_node(l.index(i,l.ny-1)).f(7)=l.get_node(l.index(i+1,l.ny-2)).f(6);
						l.get_node(l.index(i,l.ny-1)).f(4)=l.get_node(l.index(i,l.ny-2)).f(2);

						l.get_node(l.index(i,l.ny-1)).f(1)=l.get_node(l.index(i-1,l.ny-1)).f(1);
						l.get_node(l.index(i,l.ny-1)).f(3)=l.get_node(l.index(i+1,l.ny-1)).f(3);

						l.get_node(l.index(i,l.ny-1)).f(5)=l.get_node(l.index(i-1,l.ny)).f(8);
			      l.get_node(l.index(i,l.ny-1)).f(6)=l.get_node(l.index(i+1,l.ny)).f(7);
			      l.get_node(l.index(i,l.ny-1)).f(2)=l.get_node(l.index(i,l.ny)).f(4);

			      //bottom wall
						l.get_node(l.index(i,0)).f(5)=l.get_node(l.index(i-1,1)).f(8);
						l.get_node(l.index(i,0)).f(6)=l.get_node(l.index(i+1,1)).f(7);
						l.get_node(l.index(i,0)).f(2)=l.get_node(l.index(i,1)).f(4);

						l.get_node(l.index(i,0)).f(1)=l.get_node(l.index(i-1,0)).f(1);
						l.get_node(l.index(i,0)).f(3)=l.get_node(l.index(i+1,0)).f(3);

						l.get_node(l.index(i,0)).f(8)=l.get_node(l.index(i-1,-1)).f(5);
						l.get_node(l.index(i,0)).f(7)=l.get_node(l.index(i+1,-1)).f(6);
						l.get_node(l.index(i,0)).f(4)=l.get_node(l.index(i,-1)).f(2);

						/*	for (int m=0; m < velocity_set().size; ++m){
									l.get_node(i,l.ny-1).f(m)=l.get_node(i,l.ny-2).f(m);
									l.get_node(i,0).f(m)=l.get_node(i,1).f(m);
							}*/
					}
					break;
			}

			for (unsigned int i=0; i<l.ny; ++i){
					// INLET BOUNDARY - EQUILIBRIATE AT INITIAL/BOUNDARY VALUE
				auto n_inlet = l.get_node(0,i);
				n_inlet.u()   = Vmax;
				n_inlet.v()   = 0;
				n_inlet.rho() =	1;
				velocity_set().equilibrate(n_inlet);

					// OUTLET BOUNDARY - JUST TAKE VALUE OF PREVIOUS CELL
				auto n_outlet    = l.get_node(l.nx-1,i);
				auto n_outlet_m1 = l.get_node(l.nx-2,i);
				n_outlet.u()   = n_outlet_m1.u();
				n_outlet.v()   = n_outlet_m1.v();
				n_outlet.rho() = n_outlet_m1.rho();
				velocity_set().equilibrate(n_outlet);
			}

		return ;
		}


	  void collide(){

				#pragma omp parallel for
		    for (int y = 0; y<(int)l.ny; ++y){
		        for (int x = 0; x<(int)l.nx; ++x){
		            float_type rho = 0.0;
		            float_type u = 0.0;
		            float_type v = 0.0;

								auto n = l.get_node(x,y);
		            for (int i=0; i<9; ++i){
		                rho += n.f(i);
		                u   += n.f(i)*velocity_set().c[0][i];
		                v   += n.f(i)*velocity_set().c[1][i];
		            }
		            u /= rho;
		            v /= rho;

		            n.rho() = rho;
		            n.u() = u;
		            n.v() = v;
		            float_type feq[9];
		            velocity_set().f_eq(feq, rho, u, v);  //this finds the equilibrium at all physical nodes for all populations

		            //Finding the moments KBC
		            float_type dM_xy=0.0;  // the difference in moments M-M_eq
		            float_type dM_yy=0.0;
		            float_type dM_xx=0.0;
		            float_type dM_xyy=0.0;
		            float_type dM_xxy=0.0;
		            float_type dM_xxyy=0.0;

		            for (int i=0; i<9; ++i)
		            {
		                dM_xy +=velocity_set().c[0][i]*velocity_set().c[1][i]*(n.f(i)-feq[i]);
		                dM_xx +=velocity_set().c[0][i]*velocity_set().c[0][i]*(n.f(i)-feq[i]);
		                dM_yy +=velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-feq[i]);
		                dM_xyy +=velocity_set().c[0][i]*velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-feq[i]);
		                dM_xxy +=velocity_set().c[0][i]*velocity_set().c[0][i]*velocity_set().c[1][i]*(n.f(i)-feq[i]);
		                dM_xxyy +=velocity_set().c[0][i]*velocity_set().c[0][i]*velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-feq[i]);
		            }

		            dM_xy=dM_xy/rho;
		            dM_yy=dM_yy/rho;
		            dM_xx=dM_xx/rho;
		            dM_xxy=dM_xxy/rho;
		            dM_xyy=dM_xyy/rho;
		            dM_xxyy=dM_xxyy/rho;

		            float_type delS[9];
		            float_type delH[9];

								//unsigned int kbc_coluring_scheme=0;
								//kbc_coluring_scheme=1;
								switch (kbc_coluring_scheme) {
									case 0:
											//KBC minimalistic grouping // lecture 7 p.8)
											delS[0] = 0.0;
											delS[1] = 0.5*rho*0.5*(dM_xx-dM_yy);
											delS[2] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
											delS[3] = 0.5*rho*0.5*(dM_xx-dM_yy);
											delS[4] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
											delS[5] = 0.25*rho*velocity_set().c[0][5]*velocity_set().c[1][5]*dM_xy;
											delS[6] = 0.25*rho*velocity_set().c[0][6]*velocity_set().c[1][6]*dM_xy;
											delS[7] = 0.25*rho*velocity_set().c[0][7]*velocity_set().c[1][7]*dM_xy;
											delS[8] = 0.25*rho*velocity_set().c[0][8]*velocity_set().c[1][8]*dM_xy;
											break;

									case 1:
											// classical
											delS[0] = -1.0*rho*(dM_xx+dM_yy);
											delS[1] = 0.5*rho*0.5*(2*dM_xx);
											delS[2] = 0.5*rho*0.5*(2*dM_yy);
											delS[3] = 0.5*rho*0.5*(2*dM_xx);
											delS[4] = 0.5*rho*0.5*(2*dM_yy);
											delS[5] = 0.25*rho*velocity_set().c[0][5]*velocity_set().c[1][5]*dM_xy;
											delS[6] = 0.25*rho*velocity_set().c[0][6]*velocity_set().c[1][6]*dM_xy;
											delS[7] = 0.25*rho*velocity_set().c[0][7]*velocity_set().c[1][7]*dM_xy;
											delS[8] = 0.25*rho*velocity_set().c[0][8]*velocity_set().c[1][8]*dM_xy;
											break;

									default :
											std::cerr << "INVALID CHOICE OF COLORING SCHEME, CHOSING MINIMALISTIC" << '\n';
											//KBC minimalistic grouping // lecture 7 p.8)
											delS[0] = 0.0;
											delS[1] = 0.5*rho*0.5*(dM_xx-dM_yy);
											delS[2] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
											delS[3] = 0.5*rho*0.5*(dM_xx-dM_yy);
											delS[4] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
											delS[5] = 0.25*rho*velocity_set().c[0][5]*velocity_set().c[1][5]*dM_xy;
											delS[6] = 0.25*rho*velocity_set().c[0][6]*velocity_set().c[1][6]*dM_xy;
											delS[7] = 0.25*rho*velocity_set().c[0][7]*velocity_set().c[1][7]*dM_xy;
											delS[8] = 0.25*rho*velocity_set().c[0][8]*velocity_set().c[1][8]*dM_xy;
											break;

								}
								// fi = ki + si + hi
								//calculating delH = fi - fi_eq - delS
		            for (int i=0; i<9; ++i){
		                delH[i] = n.f(i)-feq[i]-delS[i];
		            }

		            float_type entScalarProd_dSdH=0;
		            float_type entScalarProd_dHdH=0;

		            // Calculating gamma
		            for (int i=0; i<9; ++i){
		                    entScalarProd_dSdH += delS[i]*delH[i]/feq[i];
		                    entScalarProd_dHdH += delH[i]*delH[i]/feq[i];
		            }

		            float_type gamma = (1.0/beta) - (2.0 - 1.0/beta)*(entScalarProd_dSdH/entScalarProd_dHdH);
		            if (entScalarProd_dHdH == 0)  gamma = 2; // gamma = 1/beta; ->	regularised 	// gamma = 2 -> LBGK
								n.gamma() = gamma;

		            for (int i=0; i<9; ++i)
		            n.f(i)=feq[i]+(1-2.0*beta)*delS[i]+(1.0-gamma*beta)*delH[i];

		        }
		    }

		    for (auto& wi : wall_nodes){
						if (wi.IS_WALL){
							wi.n.rho() = 10;
							wi.n.u() = wall_dummy;
							wi.n.v() = wall_dummy;
							for (int i=0; i<9; ++i) wi.n.f(i) = 0;
						}
		    }

		return ;
    }


		/** @brief LB step */
		void step(){

						// store populations at previous time step for force calculations //
					for (auto& bi : boundary_nodes){
							for (unsigned int m=1;m<velocity_set().size;++m){
									if (bi.distance[m]>0)
											bi.population_prev[m] = bi.n.f(m);
							}
					}
					advect();
					wall_bc();
					Force_X = 0 ;
					Force_Y = 0;
						// calculate forces using momentum exchange method
					for (auto& bi : boundary_nodes){
							for (unsigned int m=1;m<velocity_set().size;++m){
									if (bi.distance[m]>0){
											int rflct = velocity_set().rflct_latticeVelocity[m];
											Force_X += velocity_set().c[0][rflct]*( bi.n.f(rflct) + bi.population_prev[m]);
											Force_Y += velocity_set().c[1][rflct]*( bi.n.f(rflct) + bi.population_prev[m]);
									}
							}
					}

					collide();

				// write instantaneous data like Body Forces, Velocity at Probe Locations etc.
				l.write_TimeMonitors(time,Force_X,Force_Y);

				// write entire field data
				if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) ){
						std::cout << "T = "<<time<<"\tWriting File ... " << '\n';
						write_fields();
						++output_index;
				}
				++time;
		return ;
		}

		public: // write to file

			/** write macroscopic variables to ascii file */
			void write_fields(){
				std::stringstream fns;
				fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".tec";
				l.write_fields(fns.str());
			}

		public: // print

			/** print to output stream */
			friend std::ostream& operator<<(std::ostream& os, const simulation& sim){
				os << "simulation parameters\n"
				   << "---------------------\n";
				os << "domain: " << sim.l.nx << " x " << sim.l.ny << "\n";
				os << "Re:     " << sim.Re << "\n";
				os << "Vmax:   " << sim.Vmax << "\n";
				os << "visc:   " << sim.visc << "\n";
				os << "beta:   " << sim.beta << "\n";
				return os;
			}

		public: // members

			lattice l;                 ///< lattice
			std::vector<int> shift;    ///< amount of nodes to shift each population in data structure during advection
			const float_type Re;       ///< Reynolds number
			unsigned int kbc_coluring_scheme;
			unsigned int TopBotBC;
			const float_type Vmax;     ///< mean flow velocity
		  const float_type R; 			 /// Radius of Nosecone/Geometry
			const float_type x_c; 		 /// X coord of centre of Nosecone/Geometry
			const float_type y_c;			 /// Y coord of centre of Nosecone/Geometry
			const float_type Lx;			 /// Length of Nosecone
			// DECLARE XC YC P1 P2 HERE;
			float_type wall_dummy;
			const float_type visc;     ///< viscosity
			const float_type beta;     ///< LB parameter beta
			const float_type cs; 			 /// speed of sound
			const float_type slope;		 /// slope of boattail section
			float_type Force_X;				 /// Drag Force
			float_type Force_Y;				 /// Lift Force
			unsigned int time;         ///< simulation time
			bool file_output;          ///< flag whether to write files
			unsigned int output_freq;  ///< file output frequency
			unsigned int output_index; ///< index for file naming
		  std::vector<boundary_info> boundary_nodes;
		  std::vector<wall_info> wall_nodes;
		};

	} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
