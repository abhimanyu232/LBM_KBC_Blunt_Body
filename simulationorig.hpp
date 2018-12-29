/**
 *  @file
 *  @author Fabian Bösch
 *  @brief simulation
 */

#ifndef LB_SIMULATION_HPP_INCLUDED
#define LB_SIMULATION_HPP_INCLUDED

#include "H_root.hpp"
#include "lattice.hpp"
#include <sstream>

namespace lb {

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
	  visc( _Vmax*nx/_Re /*0.0000001*/),
	  beta( 1/(2*visc/std::pow(velocity_set().cs,2)+1)),
	  time(0),
	  file_output(true), // set to true if you want to write files
	  output_freq(100),
	  output_index(0)
	{
		// define amount to shift populations for advection
		for (unsigned int i=0; i<velocity_set().size; ++i)
		{
			// **************************
			// * fill in your code here *
			// **************************
			shift[i] = velocity_set().c[0][i] + l.real_nx*velocity_set().c[1][i] ;
		}
	}

	/**
	 *  @brief Initialize the flow field
	 *
	 *  Initialization includes defining initial density, velocity and
	 *  populations. You can use Taylor-Green vortex flow conditions.
	 */
	void initialize()
	{
		// **************************
		// * fill in your code here *
		// * the lines below are    *
		// * just examples          *
		// **************************

		const float_type kappa(80);
		const float_type pi(std::acos(-1.0));
		const float_type delta(0.05);
		for (int j=0; j< static_cast<int>(l.ny); ++j){
			for (int i=0; i< static_cast<int>(l.nx); ++i){
				if (j < static_cast<int>(l.ny/2.) ){
					l.get_node(i,j).u() = Vmax*std::tanh(kappa*(((float_type)j)/l.ny - 0.25));
				}
				else {
					l.get_node(i,j).u() = Vmax*std::tanh(kappa*(0.75 - ((float_type)j)/l.ny));
				}

				l.get_node(i,j).v() = delta*Vmax*std::sin(2*pi*(0.25 + ((float_type)i)/l.nx));
				l.get_node(i,j).rho() = 1.0;
				velocity_set().equilibrate(l.get_node(i,j));
			}
		}
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
		// populate buffers

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


	}

	/**  @brief apply wall boundary conditions */
	void wall_bc()
	{
		#pragma omp parallel for
		for (unsigned int i=0; i<l.wall_nodes.size(); ++i)
		{
			// **************************
			// * fill in your code here *
			// **************************
		}
	}

	/** @brief collide the populations */
	void collide()
	{

		for (int iy=0;iy< static_cast<int>(l.ny);++iy){
			for (int ix=0;ix< static_cast<int>(l.nx);++ix){
				auto n = l.get_node(ix,iy);
				n.rho() = 0.;
				float_type mom_x = 0;
				float_type mom_y = 0;

				for (unsigned int m=0; m<velocity_set().size; m++){
					n.rho() += n.f(m);
					mom_x += n.f(m)*velocity_set().c[0][m];
					mom_y += n.f(m)*velocity_set().c[1][m];
				}

				n.u() = mom_x/n.rho();
				n.v() = mom_y/n.rho();

				float_type feq[9];
				velocity_set().f_eq(feq,n.rho(),n.u(),n.v());

				for (unsigned int m=0; m < velocity_set().size; m++){
					float_type fmirr = 2*feq[m] - n.f(m);
					n.f(m) = (1-beta)*n.f(m) + beta*fmirr;
				}
			}
		}
		// **************************
		// * fill in your code here *
		// **************************
		// calculate rho u v initial, then loop over all velocities
		// equilibrate
		//repeat

	}

	/** @brief LB step */
	void step()
	{
		advect();
		wall_bc();
		collide();

		// file io
		if ( file_output && ( ((time+1) % output_freq) == 0 || time == 0 ) )
		{
			write_fields();
			++output_index;
		}

		++time;
	}

public: // write to file

	/** write macroscopic variables to ascii file */
	void write_fields()
	{
		std::stringstream fns;
		fns << "output/data_" << std::setfill('0') << std::setw(4) << output_index << ".txt";
		l.write_fields(fns.str());
	}

public: // print

	/** print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const simulation& sim)
	{
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
	const float_type visc;     ///< viscosity
	const float_type beta;     ///< LB parameter beta
	unsigned int time;         ///< simulation time
	bool file_output;          ///< flag whether to write files
	unsigned int output_freq;  ///< file output frequency
	unsigned int output_index; ///< index for file naming
};

} // lb

#endif // LB_SIMULATION_HPP_INCLUDED
