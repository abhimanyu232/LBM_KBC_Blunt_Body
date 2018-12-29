/** 
 *  @file
 *  @author Fabian Bösch
 *   @brief global typdefs etc
 */

#ifndef LB_GLOBAL_HPP_INCLUDED
#define LB_GLOBAL_HPP_INCLUDED

#include <iostream>
#include <iomanip>
#include <chrono>


/** Top level namespace */
namespace lb {

/** Floating point type (single or double precision) */
typedef float float_type;

/**
 *  @brief Coordinate in 2D
 *  @tparam Element type
 */
template<typename T>
struct coordinate
{
	/** @brief Default constructor */
	coordinate(){}
	
	/** 
	 *  @brief Construct from x and y coordinates.
	 *  @param[in] _i x coordinate
	 *  @param[in] _j y coordinate
	 */
	coordinate(T _i, T _j): i(_i), j(_j) {}
	
	/** @brief Print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const coordinate& c)
	{
		os << "{" << c.i << ", " << c.j << "}";
		return os;
	}
	
	T i; ///< x coordinate
	T j; ///< y coordinate
};

/** High resolution clock type */
typedef std::chrono::high_resolution_clock timer_type;
/** High resolution clock units type */
typedef std::chrono::duration<float, std::milli> milliseconds;
/** High resolution clock time point type */
typedef typename timer_type::time_point time_point;
/** High resolution clock time duration type */
typedef typename timer_type::duration duration;

} // lb

#endif // LB_GLOBAL_HPP_INCLUDED
