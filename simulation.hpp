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
	        node n;
	        std::array<float_type,9> distance = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
	    };

		struct wall_info{
		    node n;
		};

	  boundary_info info_n_boundary;

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
		simulation(unsigned int nx, unsigned int ny, float_type _Re, float_type _Vmax)
		: l(nx, ny),
		  shift(velocity_set().size),
		  Re(_Re),

		  Vmax(_Vmax),
	      R(l.nx/10),
	      visc(Vmax*R/Re),
	      //visc( /*fill in your code here*/ 0.001),
		  beta(1/(((2*visc)/pow(1/std::sqrt(3.0),2))+1.0)),
	      //beta( /*fill in your code here*/ 0.9),
		  time(0),
		  file_output(true), // set to true if you want to write files
		  output_freq(250),
		  output_index(0)
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
          l.get_node(i,j).rho() 	=	rho;
          velocity_set().equilibrate(l.get_node(i,j));
				}
			}


	        //we make our way through every node in the lattice from left to right and bottom to top. For every node we iterate over each direction and calculate the distance from that node to the circle in the specific direction. We solve for dx after substituting x=x0 +ci*dx into the equation for the circle. If dx is between 0 and 1 then we know that the node is close to the circle. We then figure out whether it is inside the circle or on the outside. We want the nodes that are on the outside of the circle. Once we have these we can want arrays containing the location of these nodes and their distance from the circle for each direction.

	        //define variables for position and radius of cylinder
	        //const float_type R(l.nx/8);
	        const float_type x_c(l.nx/2);
	        const float_type y_c(l.ny/2);
					const float_type Lx(x_c + l.nx/10);
	        node n_boundary;
						// coefficients of a.dx^2 + b.dx + c = 0
	        float_type a = 0;
	        float_type b = 0;
	        float_type c = 0;
	        float_type dx= 0;
	        float_type D=0;

	        	//info_n_boundary.distance = {};                      //will be of type array

	        for (int y=0; y<l.ny+1; ++y){            //direction y , y is y_0
	            for (int x=0; x<l.nx+1; ++x){        //direction x , x is x_0
									if ( x >= (x_c - 5) && x <= (Lx+5)  &&  (y >= y_c - R - 5)  && y <= (y_c + R + 5) ){
											boundary_info bi;
											bool found_intersection = false;

											if ((float_type)x <= Lx) {
												// DO THE REGULAR THING
												if ((y-y_c)*(y-y_c) - (R*R)*(x-x_c)/(Lx - x_c) > 0  ){ //check if node is inside nosecone, if it is it cant be a boundary node
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
																// ON LARGER STENCIL, WILL CHANGE AND INCLUDE dx>1 ALSO : depends on the lattice velocity i.e consider dx>1 when c[1,2][i] > 1

																//need to find a way to add the location of this node and its distance value to the array holding the boundary nodes (for each direction)
																//info_n_boundary.distance.push_back(dx);
																//info_n_boundary.n = (l,x,y);              //x refers to x0 and y refers to y0
																//info_n_boundary.distance = {dx};

																bi.distance[i] = dx;
																found_intersection = true;
															}
														}
													}
												}		else	{
													wall_info wi{l.get_node(x,y)};
													wall_nodes.push_back(wi);
												}
											}	else if ((float_type)x >= Lx) { 									// LONGER THAN NOSECONE X >= LX
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
													boundary_nodes.push_back(bi);
											}
									}
	            }
	        }
	        std::cout << "found " << boundary_nodes.size() << " boundary nodes" << std::endl;
		return ;
		}

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
					l.f[m][l.index(ix,-1)] = l.f[m][l.index(ix,l.ny-1)];
					l.f[m][l.index(ix,l.ny)] = l.f[m][l.index(ix,0)];
				}
			}
					// vertical faces faces
			for (int iy=0;iy< static_cast<int>(l.ny);++iy){
				for (unsigned int m=0; m < velocity_set().size; ++m){
					l.f[m][l.index(-1,iy)] = l.f[m][l.index(l.nx-1,iy)];
					l.f[m][l.index(l.nx,iy)] = l.f[m][l.index(0,iy)];
				}
			}
					// corners
			for (unsigned int m=0; m < velocity_set().size; ++m){
				l.f[m][l.index(-1,-1)] = l.f[m][l.index(l.nx-1,l.ny-1)];
				l.f[m][l.index(l.nx,l.ny)] = l.f[m][l.index(0,0)];
				l.f[m][l.index(l.nx,-1)] = l.f[m][l.index(0,l.ny-1)];
				l.f[m][l.index(-1,l.ny)] = l.f[m][l.index(l.nx-1,0)];
			}
	        //shift populations in the opposite directions to the velocities so that not the same population is copied across the whole lattice
			for (unsigned int m=0; m<velocity_set().size; ++m){
				if (shift[m] > 0){
					for (int n = l.index(l.nx-1,l.ny-1); n>=l.index(0,0);--n){
						l.f[m][n] = l.f[m][n-shift[m]];
					}
				}
				else if (shift[m] < 0){
					for (int n = l.index(0,0);n<=l.index(l.nx-1,l.ny-1);++n){
						l.f[m][n] = l.f[m][n-shift[m]];
					}
				}
			}

		return ;
		}

		/**  @brief apply wall boundary conditions */
		void wall_bc(){
					// bounce back bc by looping over all boundary cells
			#pragma omp parallel for
		  for (auto& bi : boundary_nodes){
		      for (int i=1;i<9;++i){
			        if (bi.distance[i] > 0)
			          bi.n.f(velocity_set().rflct_latticeVelocity[i]) = bi.n.f(i);
		      }
		  }
			    //channel flow
	        //you can turn this into a free slip condition pretty easily, just switch the populations
	        //no slip boundary condition on top and bottom wall
			for (int i=0; i<l.nx; ++i){
					  //l.get_node(i,l.ny-1).set_flag_property("wall");
					  //l.get_node(i,0).set_flag_property("wall");
					  //l.add_wall((global().c0,l.ny),(l.nx-1,l.ny));
					  //for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
		      	//top wall
					l.get_node(l.index(i,l.ny-1)).f(1)=l.get_node(l.index(i+1,l.ny-1)).f(3);
		      l.get_node(l.index(i,l.ny-1)).f(5)=l.get_node(l.index(i+1,l.ny)).f(7);
		      l.get_node(l.index(i,l.ny-1)).f(8)=l.get_node(l.index(i+1,l.ny-2)).f(6);

		      l.get_node(l.index(i,l.ny-1)).f(3)=l.get_node(l.index(i-1,l.ny-1)).f(1);
		      l.get_node(l.index(i,l.ny-1)).f(6)=l.get_node(l.index(i-1,l.ny)).f(8);
		      l.get_node(l.index(i,l.ny-1)).f(7)=l.get_node(l.index(i-1,l.ny-2)).f(5);

		      l.get_node(l.index(i,l.ny-1)).f(2)=l.get_node(l.index(i,l.ny)).f(4);
		      l.get_node(l.index(i,l.ny-1)).f(4)=l.get_node(l.index(i,l.ny-2)).f(2);

		      	//bottom wall
		      l.get_node(l.index(i,0)).f(1)=l.get_node(l.index(i+1,0)).f(3);
		      l.get_node(l.index(i,0)).f(5)=l.get_node(l.index(i+1,1)).f(7);
		      l.get_node(l.index(i,0)).f(8)=l.get_node(l.index(i+1,-1)).f(6);

		      l.get_node(l.index(i,0)).f(3)=l.get_node(l.index(i-1,0)).f(1);
		      l.get_node(l.index(i,0)).f(6)=l.get_node(l.index(i-1,1)).f(8);
		      l.get_node(l.index(i,0)).f(7)=l.get_node(l.index(i-1,-1)).f(5);

		      l.get_node(l.index(i,0)).f(2)=l.get_node(l.index(i,1)).f(4);
		      l.get_node(l.index(i,0)).f(4)=l.get_node(l.index(i,-1)).f(2);
			}

			// OVERWRITE BOUNDARY CONDITION AT EXIT - TAKE VALUE OF PREVIOUS CELL
			// OVERWRITE BOUNDARY CONDITION AT ENTRY -
		return ;
		}

		/** @brief collide the populations */

	    /*
		void collide()
		{
			// **************************
			// * fill in your code here *
			// **************************
	#pragma omp parallel for
	        for (int y = 0; y<(int)l.ny; ++y)
	        {
	            for (int x = 0; x<(int)l.nx; ++x)
	            {
	                auto n = l.get_node(x,y);


	                float_type rho = 0.0;
	                float_type u = 0.0;
	                float_type v = 0.0;
	                for (int i=0; i<9; ++i)
	                {
	                    rho += n.f(i);
	                    u   += n.f(i)*velocity_set().c[0][i];
	                    v   += n.f(i)*velocity_set().c[1][i]; //
	                }
	                u /= rho;
	                v /= rho;

	                n.rho() = rho;
	                n.u() = u;
	                n.v() = v;
	                float_type eq[9];
	                velocity_set().f_eq(eq, rho, u, v);

	                for (int i=0; i<9; ++i)
	                {
	                //insert the equation with (1-beta) from lecture 1 that updates the population     after the collision
	                n.f(i) = (1-beta)*n.f(i) + beta*(2*eq[i]-n.f(i));
	                }


	            }

	        }
	        //creates a dark space inside the cylinder, just to be able to visualize it better (also safety precaution)
	        for (auto& wi : wall_nodes)
	        {
	            wi.n.rho() = 1;
	            wi.n.u() = 0;
	            wi.n.v() = 0;
	            for (int i=0; i<9; ++i) wi.n.f(i) = 0;
	        }

		}
	     */



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
		            float_type eq[9];
		            velocity_set().f_eq(eq, rho, u, v);  //this finds the equilibrium at all physical nodes for all populations

		            //Finding the moments KBC
		            float_type dM_xy=0.0;  // the difference in moments M-M_eq
		            float_type dM_yy=0.0;
		            float_type dM_xx=0.0;
		            float_type dM_xyy=0.0;
		            float_type dM_xxy=0.0;
		            float_type dM_xxyy=0.0;

		            //float_type dM_x0=u;    Lower order moments not used
		            //float_type dM_0y=v;
		            //int dM_00=1;

		            for (int i=0; i<9; ++i)
		            {
		                dM_xy +=velocity_set().c[0][i]*velocity_set().c[1][i]*(n.f(i)-eq[i]);
		                dM_xx +=velocity_set().c[0][i]*velocity_set().c[0][i]*(n.f(i)-eq[i]);
		                dM_yy +=velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-eq[i]);
		                dM_xyy +=velocity_set().c[0][i]*velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-eq[i]);
		                dM_xxy +=velocity_set().c[0][i]*velocity_set().c[0][i]*velocity_set().c[1][i]*(n.f(i)-eq[i]);
		                dM_xxyy +=velocity_set().c[0][i]*velocity_set().c[0][i]*velocity_set().c[1][i]*velocity_set().c[1][i]*(n.f(i)-eq[i]);    //all these give M_...*rho not the actual moments
		            }

		            dM_xy=dM_xy/rho;
		            dM_yy=dM_yy/rho;
		            dM_xx=dM_xx/rho;
		            dM_xxy=dM_xxy/rho;
		            dM_xyy=dM_xyy/rho;
		            dM_xxyy=dM_xxyy/rho;

		            float_type delta_s[9];
		            float_type delta_h[9];

		            delta_s[0] = 0.0;
		            delta_s[1] = 0.5*rho*0.5*(dM_xx-dM_yy);
		            delta_s[2] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
		            delta_s[3] = 0.5*rho*0.5*(dM_xx-dM_yy);
		            delta_s[4] = 0.5*rho*0.5*-1.0*(dM_xx-dM_yy);
		            delta_s[5] = 0.25*rho*velocity_set().c[0][5]*velocity_set().c[1][5]*dM_xy;     //sigma= velocity_set.c[0][i] lambda=velocity_set.c[1][i];
		            delta_s[6] = 0.25*rho*velocity_set().c[0][6]*velocity_set().c[1][6]*dM_xy;
		            delta_s[7] = 0.25*rho*velocity_set().c[0][7]*velocity_set().c[1][7]*dM_xy;
		            delta_s[8] = 0.25*rho*velocity_set().c[0][8]*velocity_set().c[1][8]*dM_xy;

		            for (int i=0; i<9; ++i){ //KBC colision (minimalistic lecture 7 p.8)
		                //calculating deltah
		                delta_h[i] = n.f(i)-eq[i]-delta_s[i];
		            }

		            float_type delta_s_delta_h=0;
		            float_type delta_h_delta_h=0;

		            // Calculating gamma
		            for (int l=0; l<9; ++l){
		                    delta_s_delta_h=delta_s_delta_h+delta_s[l]*delta_h[l]/eq[l];
		                    delta_h_delta_h=delta_h_delta_h+delta_h[l]*delta_h[l]/eq[l];
		            }

		            float_type gamma= 1.0/beta -(2.0-1.0/beta)*(delta_s_delta_h/delta_h_delta_h);
		            if (delta_h_delta_h == 0) gamma = 2;		// LBGK

		            // Equilibrating
		            //auto f_mirr=k+(2*Delta_s)+(1-gamma*beta)*Delta_h;
		            //n.f(i)=(1-beta)*n.f(i)+beta*(f_mirr);

		            for (int i=0; i<9; ++i) //KBC colision (minimalistic lecture 7 p.8)
		            n.f(i)=eq[i]+(1-2.0*beta)*delta_s[i]+(1.0-gamma*beta)*delta_h[i];

		        }
		    }

		    for (auto& wi : wall_nodes){
		        wi.n.rho() = 1;
		        wi.n.u() = 0;
		        wi.n.v() = 0;
		        for (int i=0; i<9; ++i) wi.n.f(i) = 0;
		    }

		return ;
    }


		/** @brief LB step */
		void step(){
				advect();
				wall_bc();
				collide();

				// file io
				if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) ){
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
				fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
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
			const float_type Vmax;     ///< mean flow velocity
		    const float_type R;
			const float_type visc;     ///< viscosity
			const float_type beta;     ///< LB parameter beta
			unsigned int time;         ///< simulation time
			bool file_output;          ///< flag whether to write files
			unsigned int output_freq;  ///< file output frequency
			unsigned int output_index; ///< index for file naming
		    std::vector<boundary_info> boundary_nodes;
		    std::vector<wall_info> wall_nodes;
		};

	} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
