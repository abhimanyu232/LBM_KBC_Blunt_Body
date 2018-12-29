
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>



int main(int argc, char *argv[])
{
	omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));

	int Lx_,Ly_;
	float Re_,Vmax_;
	if (argc == 5 ){
	Lx_ = std::stoi(argv[1]);
	Ly_ = std::stoi(argv[2]);
	Re_ = std::stof(argv[3]);
	Vmax_ = std::stof(argv[4]);
	}
	else {
		std::cerr << "PLEASE PROVIDE LX LY RE AND VMAX AS COMMAND LINE ARGUEMENTS" << '\n';
		return 0;
	}
	std::cout << "Length: " << Lx_ <<"x"<<Ly_<< "\nRe:"<<Re_<< "\nVmax:" << Vmax_ << '\n';
	lb::simulation* sim = new lb::simulation(Lx_,Ly_,Re_,Vmax_); // 100,100,20000,0.04
	sim->initialize();
	std::cout << *sim << std::endl;


	#ifdef USE_OPENGL_VISUALIZATION

		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();
		for (unsigned int i=0; i<sim->l.nx/0.04; ++i)
		{
			sim->step();
		}


	#else

		// Here are some hints for getting aquainted with the lattice class
		// ================================================================

		// how to print the lattice:
		// -------------------------

		std::cout << sim->l << std::endl;

		// how to access the lattice:
		// --------------------------

		// 1) access via node proxy
		sim->l.get_node(1,0).f(0) = 2;

		// 2) access data directly (make sure you know what you're doing)
		sim->l.f[0][sim->l.index(2,0)] = 3;

		// 3) using iterators to nodes
		(sim->l.begin() + sim->l.index(0,0))->f(0) = 1;



		std::cout << sim->l << std::endl;



		// use a loop like this to run the simulation

		for (unsigned int i=0; i<sim->l.nx/0.04; ++i)
		{
			sim->step();
		}

	#endif

	return 0;
}
