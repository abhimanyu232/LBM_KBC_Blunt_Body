
#include "simulation.hpp"
#ifdef USE_OPENGL_VISUALIZATION
#include "visualization.hpp"
#endif
#include <omp.h>



int main(int argc, char *argv[])
{
	omp_set_num_threads(std::max(omp_get_max_threads(),omp_get_num_procs()));

	int Lx_,Ly_;
	unsigned int KBC_COLOR_;
	unsigned int BC_topbottom;
	float Re_,Vmax_;
	bool WRITE_FILE_;
	int iter = 5000;
	if (argc == 8 ){
	Lx_ = std::stoi(argv[1]);
	Ly_ = std::stoi(argv[2]);
	Re_ = std::stof(argv[3]);
	Vmax_ = std::stof(argv[4]);
	KBC_COLOR_ = std::stoi(argv[5]);	 // 0 : MINIMALISTIC  // 1 : CLASSICAL
	BC_topbottom = std::stoi(argv[6]); // 0 : No Slip 			// 1 : Slip
	WRITE_FILE_ = std::stoi(argv[7]);	 // 0 : NO WRITE 			// 1 : WRITE
	}
	else {
		std::cerr << "PLEASE PROVIDE LX LY RE VMAX KBC_COLOR_SCHEME(0:MINIMALISTIC ; 1:CLASSICAL) TOP_BOTTOM_BC(0:NO SLIP ; 1:SLIP) AND WHETHER TO WRITE DATA TO FILE(0:NO WRITE ; 1:WRITE) AS COMMAND LINE ARGUEMENTS" << '\n';
		return 0;
	}

	if (KBC_COLOR_ < 0 || KBC_COLOR_ > 1 ){
			std::cout << "Invalid Choice for Coloring Scheme \t Please from Below \n"
			"0: MINIMALISTIC GROUPING \n1: CLASSICAL GROUPING\n KBC_COLOR_= " ;
			std::cin >> KBC_COLOR_;
	}

	// Initialize Time Series Output Needed to Calculate Vortex Shedding Frequency
	std::ofstream ofile("output/monitors.txt", std::ios::out );
	if (!ofile.is_open()){ throw std::runtime_error("could not write to file");	}

	lb::simulation* sim = new lb::simulation(Lx_,Ly_,Re_,Vmax_,KBC_COLOR_,BC_topbottom,WRITE_FILE_); // 100,100,20000,0.04
	sim->initialize();
	std::cout << *sim << std::endl;


	#ifdef USE_OPENGL_VISUALIZATION

		lb::visualization::initialize(sim,argc,argv);
		lb::visualization::get_instance().run();

	#else

		// Here are some hints for getting aquainted with the lattice class
		// ================================================================

		// how to print the lattice:
		// -------------------------
		/*
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
		*/


		// use a loop like this to run the simulation

		for (unsigned int i=0; i<sim->l.nx/0.04; ++i)
		{
			sim->step();
		}

	#endif

	return 0;
}
