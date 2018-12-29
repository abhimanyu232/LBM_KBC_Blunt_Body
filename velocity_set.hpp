/**
 *  @file
 *  @author Fabian Bösch
 *  @brief velocity set
 */

#ifndef LB_VELOCITY_SET_HPP_INCLUDED
#define LB_VELOCITY_SET_HPP_INCLUDED

#include "global.hpp"
#include <array>
#include <cmath>

namespace lb {

struct v9;                // forward declaration
const v9& velocity_set(); // forward declaration

/**
 *  @brief Lattice parameters for 9 velocity model.
 *
 *  This class models a the singleton design pattern. That means there
 *  exists only one single instance throughout the lifetime of the
 *  program. To instantiate and access this object use the free function
 *  @ref velocity_set.
 *
 *  This class holds parameters like lattice weights, molecular
 *  velocities and speed of sound. It also exposes member functions to
 *  compute the equilibrium populations.
 */
struct v9 // singleton
{
private:

	/** @brief Default constructor */
	v9(){};
	/** @brief Function for instantiating the singleton is a friend */
	friend const v9& lb::velocity_set();

public:

	v9(const v9&) = delete;
	v9& operator=(const v9&) = delete;


	//                                                     0,       1,       2,       3,       4,       5,       6,       7,       8
	const std::array<float_type, 9>         W =   {{ 16.0/36,  4.0/36,  4.0/36,  4.0/36,  4.0/36,  1.0/36,  1.0/36,  1.0/36,  1.0/36}};   ///< Lattice weights

	const std::array<std::array<int, 9>, 2> c = {{{{       0,       1,       0,      -1,       0,       1,      -1,      -1,       1}},
	                                              {{       0,       0,       1,       0,      -1,       1,       1,      -1,      -1}}}}; ///< Molecular velocities
	const std::array<int, 9> opposite_c_index= {{0, 3, 4, 1, 2, 7, 8, 5, 6}};

	const float_type cs = 1.0/std::sqrt(3.0);   ///< Speed of sound

	const unsigned int size = 9;                ///< Number of velocities

	/**
	 *  @brief Compute equilibrium.
	 *
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  (see also @ref v9::equilibrate).
	 *  @param[in,out] f_eq Pointer to an array of size 9 to store the computed values
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	inline void f_eq(float_type* f_eq, float_type rho, float_type u, float_type v) const
	{
		// **************************
		// * fill in your code here *
		// * the lines below are    *
		// * just examples          *
		// **************************
		f_eq[0] = rho*W[0]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v));
		f_eq[1] = rho*W[1]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*(2*u+sqrt(1+3*u*u))/(1-u);
		f_eq[2] = rho*W[2]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*(2*v+sqrt(1+3*v*v))/(1-v);
		f_eq[3] = rho*W[3]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*pow((2*u+sqrt(1+3*u*u))/(1-u),-1);
		f_eq[4] = rho*W[4]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*pow((2*v+sqrt(1+3*v*v))/(1-v),-1);
		f_eq[5] = rho*W[5]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*
											((2*u+sqrt(1+3*u*u))/(1-u))*((2*v+sqrt(1+3*v*v))/(1-v));
		f_eq[6] = rho*W[6]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*
											pow((2*u+sqrt(1+3*u*u))/(1-u),-1)*(2*v+sqrt(1+3*v*v))/(1-v);
		f_eq[7] = rho*W[7]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*
											pow((2*u+sqrt(1+3*u*u))/(1-u),-1)*pow((2*v+sqrt(1+3*v*v))/(1-v),-1);
		f_eq[8] = rho*W[8]*(2-sqrt(1+3*u*u))*(2-sqrt(1+3*v*v))*
											(2*u+sqrt(1+3*u*u))/(1-u)*pow((2*v+sqrt(1+3*v*v))/(1-v),-1);
	}

	/**
	 *  @brief Equilibrate a node.
	 *
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also
	 *  @ref v9::f_eq).
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 *  @param[in]     rho  Local density
	 *  @param[in]     u    Local flow velocity in x-direction
	 *  @param[in]     v    Local flow velocity in y-direction
	 */
	template <typename Node>
	inline void equilibrate(Node& n, float_type rho, float_type u, float_type v) const
	{
		float_type feq[9];
		f_eq(feq, rho, u, v);
		n.f(0) = feq[0];
		n.f(1) = feq[1];
		n.f(2) = feq[2];
		n.f(3) = feq[3];
		n.f(4) = feq[4];
		n.f(5) = feq[5];
		n.f(6) = feq[6];
		n.f(7) = feq[7];
		n.f(8) = feq[8];
	}

	/**
	 *  @brief Equilibrate a node.
	 *
	 *  Compute f_eq from the locally conserved quantities rho, u and v
	 *  and set the node's population to that equilibrium ( see also
	 *  @ref v9::f_eq and v9::equilibrate). The locally conserved
	 *  quantities are taken form the node object itself.
	 *  @tparam        Node A node type
	 *  @param[in,out] n    Reference to a Node object
	 */
	template <typename Node>
	inline void equilibrate(Node& n) const
	{
		return equilibrate(n, n.rho(), n.u(), n.v());
	}
};

/**
 *  @brief Get a reference single instance of the velocity set.
 *  @return 9-velocity set
 */
inline const v9& velocity_set()
{
	static v9 v_set;
	return v_set;
}

} // lb

#endif // LB_VELOCITY_SET_HPP_INCLUDED
