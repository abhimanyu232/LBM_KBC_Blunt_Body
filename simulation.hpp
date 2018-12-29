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
	      R(l.ny/8.0),
	      visc(Vmax*R/Re),
	      //visc( /*fill in your code here*/ 0.001),
		  beta(1/(((2*visc)/pow(1/std::sqrt(3.0),2))+1.0)),
	      //beta( /*fill in your code here*/ 0.9),
		  time(0),
		  file_output(false), // set to true if you want to write files
		  output_freq(10),
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

	        float_type rho=1.0;

			for (int j=0; j<static_cast<int>(l.ny); ++j){
				for (int i=0; i<static_cast<int>(l.nx); ++i){
	        l.get_node(i,j).u()   = Vmax;
					l.get_node(i,j).v()   = 0;
					/* 	taylor green vortex
					//l.get_node(i,j).rho() = 1.0;
	        if (j>static_cast<int>(l.ny/2)){
						l.get_node(i,j).u()   = 0.005;
	          l.get_node(i,j).u()   = Vmax*std::tanh(kap*(0.75-(1.0*j)/l.nx));
	          //l.get_node(i,j).v() = delta*Vmax*std::sin(2*pi*((1.0*i)/l.ny+0.25));
	          l.get_node(i,j).v()		= 0;
					} else	{
            l.get_node(i,j).u()   = Vmax*std::tanh(kap*((1.0*j)/l.nx-0.25));
            l.get_node(i,j).v() 	= delta*Vmax*std::sin(2*pi*((1.0*i)/l.ny+0.25));
            l.get_node(i,j).v()		= 0;
          }*/
          l.get_node(i,j).rho() 	=	rho;
          velocity_set().equilibrate(l.get_node(i,j));
				}
			}


	        //we make our way through every node in the lattice from left to right and bottom to top. For every node we iterate over each direction and calculate the distance from that node to the circle in the specific direction. We solve for t after substituting x=x0 +ci*t into the equation for the circle. If t is between 0 and 1 then we know that the node is close to the circle. We then figure out whether it is inside the circle or on the outside. We want the nodes that are on the outside of the circle. Once we have these we can want arrays containing the location of these nodes and their distance from the circle for each direction.

	        //define variables for position and radius of cylinder
	        //const float_type R(l.nx/8);
	        const float_type x_c(l.nx/2);
	        const float_type y_c(l.ny/2);
	        node n_boundary;

	        float_type A=0;
	        float_type B=0;
	        float_type C=0;
	        float_type t=0;
	        float_type discr=0;

	        //info_n_boundary.distance = {};                      //will be of type array

	        for (int y=0; y<l.ny+1; ++y){            //direction y , y is y_0
	            for (int x=0; x<l.nx+1; ++x){        //direction x , x is x_0

	                boundary_info bi;
	                bool found_intersection = false;

	                if ((x-x_c)*(x-x_c) + (y-y_c)*(y-y_c)>R*R){ //check if node is inside cylinder, if it is it cant be a boundary node
	                    for (int i=1; i<9; ++i){        //9 directions
	                        A = velocity_set().c[0][i]*velocity_set().c[0][i] + velocity_set().c[1][i]*velocity_set().c[1][i];
	                        B = 2*velocity_set().c[0][i]*(x - x_c) + 2*velocity_set().c[1][i]*(y - y_c);
	                        C = x*x + y*y + x_c*x_c + y_c*y_c -2*x*x_c -2*y*y_c - R*R;
	                        discr = B*B-4*A*C;

	                        if (discr>=0){
	                            const auto t1 =(-B + std::sqrt(B*B-4*A*C))/(2*A);
	                            const auto t2 =(-B - std::sqrt(B*B-4*A*C))/(2*A);
	                            t = 0;
	                            if (t1 > 0) t = t1;
	                            if (t2 > 0 && t2 < t1) t = t2;

	                            if (t>0 && t<=1){
																	// ON LARGER STENCIL, WILL CHANGE AND INCLUDE t>1 ALSO : depends on the lattice velocity i.e consider t>1 when c[1,2][i] > 1

	                                //need to find a way to add the location of this node and its distance value to the array holding the boundary nodes (for each direction)
	                                //info_n_boundary.distance.push_back(t);
	                                //info_n_boundary.n = (l,x,y);              //x refers to x0 and y refers to y0
	                                //info_n_boundary.distance = {t};

                                bi.distance[i] = t;
                                found_intersection = true;
	                            }
	                        }
	                    }
	                }	else	{
	                    wall_info wi{l.get_node(x,y)};
	                    wall_nodes.push_back(wi);
	                }
	                if (found_intersection){
	                    bi.n = l.get_node(x,y);
	                    boundary_nodes.push_back(bi);
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
		void advect()
		{
			// **************************
			// * fill in your code here *
			// **************************

	        //first initialize buffers

	        //bottom buffer
	        for (int b=0; b<l.nx; ++b)
	        {
	            for (int k=0; k<9; ++k)
	            {
	                l.get_node(l.index(b,-1)).f(k)=l.get_node(l.index(b,l.ny-1)).f(k);
	            }
	        }

	        //top buffer
	        for (int b=0; b<l.nx; ++b)
	        {
	            for (int k=0; k<9; ++k)
	            {
	                l.get_node(l.index(b,l.ny)).f(k)=l.get_node(l.index(b,0)).f(k);
	            }
	        }

	        //left buffer
	        for (int b=0; b<l.ny; ++b)
	        {
	            for (int k=0; k<9; ++k)
	            {
	                l.get_node(l.index(-1,b)).f(k)=l.get_node(l.index(l.nx-1,b)).f(k);
	            }
	        }

	        //right buffer
	        for (int b=0; b<l.ny; ++b)
	        {
	            for (int k=0; k<9; ++k)
	            {
	                l.get_node(l.index(l.nx,b)).f(k)=l.get_node(l.index(0,b)).f(k);
	            }
	        }

	        //bottom left corner buffer
	        for (int k=0; k<9; ++k)
	        {
	            l.get_node(l.index(-1,-1)).f(k)=l.get_node(l.index(l.nx-1,l.ny-1)).f(k);
	        }

	        //bottom right corner buffer
	        for (int k=0; k<9; ++k)
	        {
	            l.get_node(l.index(l.nx,-1)).f(k)=l.get_node(l.index(0,l.ny-1)).f(k);
	        }

	        //top right corner buffer
	        for (int k=0; k<9; ++k)
	        {
	            l.get_node(l.index(l.nx,l.ny)).f(k)=l.get_node(l.index(0,0)).f(k);
	        }

	        //top left corner buffer
	        for (int k=0; k<9; ++k)
	        {
	            l.get_node(l.index(-1,l.ny)).f(k)=l.get_node(l.index(l.nx-1,0)).f(k);
	        }


	        //
	        //shift populations in the opposite directions to the velocities so that not the same population is copied across the whole lattice
	        for (int j=l.ny-1; j>-1; --j)
	        {
	            for (int i=0; i<l.nx; ++i)
	            {
	                //f2 and f6
	                l.get_node(l.index(i,j)).f(2)=l.get_node(l.index(i,j)-shift[2]).f(2);
	                l.get_node(l.index(i,j)).f(6)=l.get_node(l.index(i,j)-shift[6]).f(6);
	            }
	        }

	        for (int j=l.ny-1; j>-1; --j)
	        {
	            for (int i=l.nx-1; i>-1; --i)
	            {
	                //f1 and f5
	                l.get_node(l.index(i,j)).f(1)=l.get_node(l.index(i,j)-shift[1]).f(1);
	                l.get_node(l.index(i,j)).f(5)=l.get_node(l.index(i,j)-shift[5]).f(5);
	            }
	        }

	        for (int j=0; j<l.ny; ++j)
	        {
	            for (int i=l.nx-1; i>-1; --i)
	            {
	                //f4 and f8
	                l.get_node(l.index(i,j)).f(4)=l.get_node(l.index(i,j)-shift[4]).f(4);
	                l.get_node(l.index(i,j)).f(8)=l.get_node(l.index(i,j)-shift[8]).f(8);
	            }
	        }

	        for (int j=0; j<l.ny; ++j)
	        {
	            for (int i=0; i<l.nx; ++i)
	            {
	                //f3 and f7
	                l.get_node(l.index(i,j)).f(3)=l.get_node(l.index(i,j)-shift[3]).f(3);
	                l.get_node(l.index(i,j)).f(7)=l.get_node(l.index(i,j)-shift[7]).f(7);
	            }
	        }

		return ;
		}

		/**  @brief apply wall boundary conditions */
		void wall_bc(){
			  for (auto& bi : boundary_nodes){
			      for (int i=1;i<9;++i){
				        if (bi.distance[i] > 0)
				          bi.n.f(velocity_set().opposite_c_index[i]) = bi.n.f(i);
									//bi.n.f(velocity_set().opposite_c_index(i)) = bi.n.f(i);
			      }
			  }
			//#pragma omp parallel for
			//for (unsigned int i=0; i<l.wall_nodes.size(); ++i){
		        /*
		        //flat plate obstacle in middle of flow
		        int q = l.nx/2;
		        int nobst = l.ny/3;
		        int jbot = l.ny/2 - nobst/2;
		        int jtop = l.ny/2 + nobst/2 +1;

		        for (int j=jbot; j<jtop+1; ++j){
		            l.get_node(l.index(q,j)).f(1)=l.get_node(l.index(q+1,j)).f(3);
		            l.get_node(l.index(q,j)).f(5)=l.get_node(l.index(q+1,j+1)).f(7);
		            l.get_node(l.index(q,j)).f(8)=l.get_node(l.index(q+1,j-1)).f(6);
		            l.get_node(l.index(q,j)).f(3)=l.get_node(l.index(q-1,j)).f(1);
		            l.get_node(l.index(q,j)).f(7)=l.get_node(l.index(q-1,j)).f(5);
		            l.get_node(l.index(q,j)).f(6)=l.get_node(l.index(q-1,j)).f(8);

		            l.get_node(l.index(q,jtop)).f(2)=l.get_node(l.index(q,jtop+1)).f(4);
		            l.get_node(l.index(q,jtop)).f(6)=l.get_node(l.index(q-1,jtop+1)).f(8);

		            l.get_node(l.index(q,jbot)).f(4)=l.get_node(l.index(q,jbot-1)).f(2);
		            l.get_node(l.index(q,jbot)).f(7)=l.get_node(l.index(q-1,jbot-1)).f(5);
		        }
				}*/

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
	        // ************************
	        // * fill in your code here *
	        // ************************
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
	                if (delta_h_delta_h == 0) gamma = 2;

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
