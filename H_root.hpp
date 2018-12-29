/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief compute root of H function
 */

#ifndef LB_H_ROOT_HPP_INCLUDED
#define LB_H_ROOT_HPP_INCLUDED

#include "velocity_set.hpp"
#include <algorithm>

namespace lb {

/**
 *  @brief Find over-relaxation parameter alpha
 * 
 *  This function looks for a state with the same entropy/H-function 
 *  value as the current state in the hyperplane connecting the current
 *  f with the f_eq.
 * 
 *  This procedure involves finding the root of the entropy function:
 *  H(f + alpha*(f_eq - f)) - H(f) = 0
 *  
 *  You should also account for preventing the populations from becoming
 *  negative and for performance reasons not computing alpha for states
 *  close to the equilibrium.
 *  
 *  For the root finding of the nonlinear equation use Newton-Raphson.
 * 
 *  @tparam Node node type
 *  @param[in] n node object
 *  @return alpha
 */
template <typename Node>
inline float_type H_root(const Node& n)
{
	// **************************
	// * fill in your code here *
	// **************************
	
	return 2.0;
}

} // lb

#endif // LB_H_ROOT_HPP_INCLUDED
