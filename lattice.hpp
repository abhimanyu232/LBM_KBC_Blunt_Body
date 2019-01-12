/**
 *  @file
 *  @author Fabian Bösch
 *  @brief lattice and node
 */

#ifndef LB_LATTICE_HPP_INCLUDED
#define LB_LATTICE_HPP_INCLUDED

#include "velocity_set.hpp"
//#include "distance.hpp"
#include "property_array.hpp"
#include <vector>
#include <fstream>

namespace lb {

class lattice; // forward declaration

/**
 *  @brief Node representing one lattice site.
 *
 *  Easy access to bundled quantities and properties (works as proxy to
 *  the lattice class).
 */
struct node
{
public: // ctors

	/** @brief Default constructor */
	node() {};

	/**
	 *  @brief Construct from lattice and position
	 *  @param[in] lat Pointer to the lattice
	 *  @param[in] i   x coordinate
	 *  @param[in] j   y coordinate
	 *  @pre coordinates are in domain
	 */
	node(lattice* lat, int i, int j);

	node(const node&) = default;

public: // init

	/**
	 *  @brief Set lattice and position.
	 *  @param[in] lat Pointer to lattice
	 *  @param[in] i   x coordinate
	 *  @param[in] j   y coordinate
	 *  @pre coordinates are in domain
	 */
	void set(lattice* lat, int i, int j);

public: // access populations and macroscopic quantities

	/**
	 *  @brief Get population.
	 *  @param i Population index
	 *  @return Value of distribution function
	 *  @pre population index exists
	 */
	inline float_type f(unsigned int i) const;

	/**
	 *  @brief Get/set population.
	 *  @param i Population index
	 *  @return Reference to value of distribution function
	 *  @pre population index exists
	 */
	inline float_type& f(unsigned int i);

	/**
	 *  @brief Get density.
	 *  @return Local density
	 */
	inline float_type rho() const;

	/**
	 *  @brief Get/set density.
	 *  @return Reference to local density
	 */
	inline float_type& rho();

	/**
	 *  @brief Get x-velocity.
	 *  @return Local flow velocity in x direction
	 */
	inline float_type u() const;

	/**
	 *  @brief Get/set x-velocity.
	 *  @return Reference to local flow velocity in x direction
	 */
	inline float_type& u();

	/**
	 *  @brief Get y-velocity.
	 *  @return Local flow velocity in y direction
	 */
	inline float_type v() const;

	/**
	 *  @brief Get/set y-velocity.
	 *  @return Reference to local flow velocity in y direction
	 */
	inline float_type& v();

	/**
	 *  @brief Get overRelaxation.
	 *  @return local overRelaxation
	 */
	inline float_type gamma() const;

	/**
	 *  @brief Get/set overRelaxation.
	 *  @return Reference to local overRelaxation
	 */
	inline float_type& gamma();


public: // query and access properties

	/**
	 *  @brief Query for flag property.
	 *  Query whether a flag is set for the node.
	 *  @param[in] name Flag name
	 *  @return True if flag is set, otherwise false
	 */
	inline bool has_flag_property(std::string name) const;

	/**
	 *  @brief Set a flag
	 *  Set the flag "name" to true
	 *  @param[in] name Flag name
	 *  @return True if flag exists, otherwise false
	 */
	inline bool set_flag_property(std::string name);

	/**
	 *  @brief Unset a flag
	 *  Set the flag "name" to false
	 *  @param[in] name Flag name
	 *  @return True if flag exists, otherwise false
	 */
	inline bool unset_flag_property(std::string name);

	/**
	 *  @brief Query for data property.
	 *  Query whether data property (object) is stored for the node.
	 *  @param[in] name Data property name
	 *  @return True if thre is such a data property, otherwise false
	 */
	inline bool has_data_property(std::string name) const;

	/**
	 *  @brief Store a data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @param[in] property Data property object
	 *  @return True if data property exists, otherwise false
	 */
	template<typename T>
	inline bool set_data_property(std::string name, const T& property);

	/**
	 *  @brief Delete a data property
	 *  @param[in] name Data property name
	 *  @return True if data property exists, otherwise false
	 */
	bool unset_data_property(std::string name);

	/**
	 *  @brief Get data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @return Reference to data property object
	 */
	template <typename T>
	T& get_data_property(std::string name);

	/**
	 *  @brief Get data property
	 *  @tparam T Type of the data property
	 *  @param[in] name Data property name
	 *  @return Reference to data property object
	 */
	template <typename T>
	const T& get_data_property(std::string name) const;

public: // members

	lattice* l;            ///< Pointer to a lattice object
	unsigned int index;    ///< Index for looking up data in the lattice
	coordinate<int> coord; ///< Coordinate of node's position
};

/**
 *  @brief Lattice containing the populations.
 *
 *  The lattice is constructed using the function @ref velocity_set()
 *  which returns a velocity set object. Hence, the number of
 *  populations is defined through that function. Data structures are
 *  set up accordingly.
 *
 *  The basic data structure for the population and the macroscopic
 *  qunatities are one dimensional arrays (vectors) interpreted as two
 *  dimensional planes. The x (i) dimension varies first and the y (j)
 *  dimension last.
 *
 *  This class does provide access to the data through node iterators or
 *  through direct access of the public members. The node iterators
 *  return a @ref node object that provides easy access to all local
 *  quantities according to the 2d lattice coordinate.
 *
 *  There are buffer regions (extent is one in all directions) around
 *  the data to make the advection procedure easier.
 *
 *  The data is indexed in the range [0, nx-1][0, ny-1]; including
 *  buffers indices span the range [-1, nx][-1, ny], repectively.
 */
class lattice
{
public: // typedefs

	/** @brief Iterator type */
	typedef typename std::vector<node>::iterator node_iterator;
	/** @brief Const iterator type */
	typedef typename std::vector<node>::const_iterator const_node_iterator;
	/** @brief Reverse iterator type */
	typedef typename std::vector<node>::reverse_iterator reverse_node_iterator;
	/** @brief Const reverse iterator type */
	typedef typename std::vector<node>::const_reverse_iterator const_reverse_node_iterator;

public: // ctor

	/**
	 *  @brief Construct the lattice with given extent
	 *  @param[in] _nx Number of nodes in x direction
	 *  @param[in] _ny Number of nodes in y direction
	 */
	lattice(unsigned int _nx, unsigned int _ny);

public: // coordinates to index conversion

	/**
	 *  @brief Convert a coordinate to a unique index
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return unique index
	 *  @pre Coordinates are in the domain
	 */
	inline unsigned int index(int i, int j) const;

public: // node access

	/** @brief Iterator pointing to the beginning  @return iterator */
	node_iterator begin();
	/** @brief Const iterator pointing to the beginning @return const iterator */
	const_node_iterator begin() const;
	/** @brief Iterator pointing to the end @return iterator */
	node_iterator end();
	/** @brief Const iterator pointing to the end @return const iterator */
	const_node_iterator end() const;
	/** @brief Reverse iterator pointing to the end @return reverse iterator */
	reverse_node_iterator rbegin();
	/** @brief Const reverse iterator pointing to the end @return const reverse iterator */
	const_reverse_node_iterator rbegin() const;
	/** @brief Reverse iterator pointing to the beginning @return reverse iterator */
	reverse_node_iterator rend();
	/** @brief Const reverse iterator pointing to the beginning @return const reverse iterator */
	const_reverse_node_iterator rend() const;

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return reference to node at coordinate (i,j)
	 *  @pre coordinates are in domain
	 */
	inline node& get_node(int i, int j);

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] i x coordinate
	 *  @param[in] j y coordinate
	 *  @return const reference to node at coordinate (i,j)
	 *  @pre coordinates are in domain
	 */
	inline const node& get_node(int i, int j) const;

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] idx unique node index
	 *  @return reference to node at coordinate (i,j)
	 *  @pre idx is between [0, @ref lattice::real_size )
	 */
	inline node& get_node(unsigned int idx);

	/**
	 *  @brief Get node at coordinate (i,j)
	 *  @param[in] idx unique node index
	 *  @return const reference to node at coordinate (i,j)
	 *  @pre idx is between [0, @ref lattice::real_size )
	 */
	inline const node& get_node(unsigned int idx) const;

public: // walls

	/**
	 *  @brief Add a solid wall
	 *
	 *  Creates wall flags in the coordinate rectangle defined by
	 *  min_coord and max_coord. The corresponding nodes get the flag
	 *  "wall" and they are also stored in the vector
	 *  @ref lattice:wall_nodes for convienience.
	 *
	 *  @param[in] min_coord minimum bounding rectangle corner
	 *  @param[in] max_coord maximum bounding rectangle corner
	 *  @pre (min_coord, max_coord) define a rectangle
	 *  @pre Both min_coord and max_coord are in the domain
	 */
	void add_wall(coordinate<int> min_coord, coordinate<int> max_coord);

	/** @brief Delete all existing walls */
	void delete_walls();

public: // file dump

	/**
	 *  @brief Write fields to file
	 *
	 *  Write macroscopic variables to simple ascii file.
	 *
	 *  @param[in] file_name file name
	 */
	void write_fields(std::string file_name);

	/**
	 *  @brief Write Ux at specified probes to file
	 *
	 *  Write macroscopic variables to simple ascii file.
	 *
	 */
	void write_TimeMonitors(unsigned int time,float_type Fx, float_type Fy);

public: // print

	/** @brief print to output stream, useful for debugging only */
	friend std::ostream& operator<<(std::ostream& os, const lattice& l);

public: // members

	const unsigned int nx;                    ///< extent in x direction (excluding buffers)
	const unsigned int ny;                    ///< extent in y direction (excluding buffers)
	const unsigned int size;                  ///< total number of nodes (excluding buffers)
	const unsigned int buffer_size;           ///< buffer width (equal to one)
	const unsigned int real_nx;               ///< extent in x direction including buffers
	const unsigned int real_ny;               ///< extent in y direction including buffers
	const unsigned int real_size;             ///< total number of nodes including buffers
	const unsigned int n_populations;         ///< number of populations
	std::vector<std::vector<float_type> > f;  ///< population data
	std::vector<float_type> rho;              ///< density data
	std::vector<float_type> u;                ///< flow x-velocity data
	std::vector<float_type> v;                ///< flow y-velocity data
	std::vector<float_type> gamma;            ///< flow y-velocity data
	std::vector<node> nodes;                  ///< array holding all node objects
	std::vector<node> wall_nodes;             ///< array holding node objects belonging to a solid wall
	property_array properties;                ///< properties datastructure (can hold many different properties per node)
	const bool periodic_x;                    ///< flag whether to use periodicity in x direction
	const bool periodic_y;                    ///< flag whether to use periodicity in y direction
};

// implementation
// --------------

// node

node::node(lattice* lat, int i, int j)
: l(lat), index(l->real_nx*(j+l->buffer_size) + i + l->buffer_size), coord(i,j) { }

void node::set(lattice* lat, int i, int j)
{
	l = lat;
	index = l->real_nx*(j+l->buffer_size) + i + l->buffer_size;
	coord.i = i;
	coord.j = j;
}

inline float_type node::f(unsigned int i) const { return l->f[i][index]; }
inline float_type& node::f(unsigned int i) { return l->f[i][index]; }
inline float_type node::rho() const { return l->rho[index]; }
inline float_type& node::rho() { return l->rho[index]; }
inline float_type node::u() const { return l->u[index]; }
inline float_type& node::u() { return l->u[index]; }
inline float_type node::v() const { return l->v[index]; }
inline float_type& node::v() { return l->v[index]; }
inline float_type node::gamma() const { return l->gamma[index]; }
inline float_type& node::gamma() { return l->gamma[index]; }

inline bool node::has_flag_property(std::string name) const { return l->properties.has_flag_property(name, index); }
inline bool node::set_flag_property(std::string name) { return l->properties.set_flag_property(name, index); }
inline bool node::unset_flag_property(std::string name) { return l->properties.unset_flag_property(name, index); }
inline bool node::has_data_property(std::string name) const { return l->properties.has_data_property(name,index); }
template<typename T>
inline bool node::set_data_property(std::string name, const T& property) { return l->properties.set_data_property(name, index, property); }
bool node::unset_data_property(std::string name) { return l->properties.unset_data_property(name, index); }
template <typename T>
T& node::get_data_property(std::string name) { return l->properties.get_data_property<T>(name, index); }
template <typename T>
const T& node::get_data_property(std::string name) const { return l->properties.get_data_property<T>(name, index); }



// lattice

lattice::lattice(unsigned int _nx, unsigned int _ny)
: nx(_nx), ny(_ny), size(nx*ny), buffer_size(1), real_nx(nx+2*buffer_size), real_ny(ny+2*buffer_size),
  real_size(real_nx*real_ny), n_populations(velocity_set().size),
  f( n_populations, std::vector<float_type>(real_size, 0) ),
  rho(real_size, 0), u(real_size, 0), v(real_size, 0), gamma(real_size,0), nodes(real_size),
  properties(real_size), periodic_x(true), periodic_y(true)
{
	// register some properties
	properties.register_flag_property("fluid");
	properties.register_flag_property("buffer");
	properties.register_flag_property("wall");

	// set up nodes and properties
	unsigned int k(0);
	for (unsigned int j=0; j<real_ny; ++j)
	{
		for (unsigned int i=0; i<real_nx; ++i)
		{
			nodes[k].set(this, static_cast<int>(i)-buffer_size, static_cast<int>(j)-buffer_size);
			if (i<buffer_size || i>=real_nx-buffer_size || j<buffer_size || j>=real_ny-buffer_size)
				properties.set_flag_property("buffer",nodes[k].index);
			else properties.set_flag_property("fluid",nodes[k].index);
			++k;
		}
	}
}

lattice::node_iterator lattice::begin() { return nodes.begin(); }
lattice::const_node_iterator lattice::begin() const { return nodes.begin(); }
lattice::node_iterator lattice::end() { return nodes.end(); }
lattice::const_node_iterator lattice::end() const { return nodes.end(); }
lattice::reverse_node_iterator lattice::rbegin() { return nodes.rbegin(); }
lattice::const_reverse_node_iterator lattice::rbegin() const { return nodes.rbegin(); }
lattice::reverse_node_iterator lattice::rend() { return nodes.rend(); }
lattice::const_reverse_node_iterator lattice::rend() const { return nodes.rend(); }

inline unsigned int lattice::index(int i, int j) const { return real_nx*(j+buffer_size) + i + buffer_size; }

inline node& lattice::get_node(int i, int j) { return nodes[real_nx*(j+buffer_size) + i + buffer_size]; }

inline const node& lattice::get_node(int i, int j) const { return nodes[real_nx*(j+buffer_size) + i + buffer_size]; }

inline node& lattice::get_node(unsigned int idx) { return nodes[idx]; }

inline const node& lattice::get_node(unsigned int idx) const { return nodes[idx]; }


std::ostream& operator<<(std::ostream& os, const lattice& l)
{
	for (unsigned int p=0; p<l.n_populations; ++p)
	{
		os << " f" << std::setw(2) << p << ":";
		os << "\n" << "   y  " << std::setw((l.nx+2*l.buffer_size)*12+1) << std::setfill('-') << "" << "\n" << std::setfill(' ');
		for (int j=static_cast<int>(l.ny+l.buffer_size-1); j>-static_cast<int>(l.buffer_size+1); --j)
		{
			os << std::setw(4) << j << " |";
			for (int i=-static_cast<int>(l.buffer_size); i<static_cast<int>(l.nx+l.buffer_size); ++i)
			{
				const unsigned int index = (j+l.buffer_size)*l.real_nx + i + l.buffer_size;
				if (i>=0 && i<static_cast<int>(l.nx) && j>=0 && j<static_cast<int>(l.ny))
				os << std::setw(12) << std::setprecision(5) /*<< std::scientific*/ << l.f[p][index];
				else
				os << std::setw(12) << "*";
			}
			os << " |" <<  "\n";
		}
		os << std::setw(6) << "" << std::setw((l.nx+2*l.buffer_size)*12+1) << std::setfill('-') << "" << "\n" << std::setfill(' ') << std::setw(6) << "";
		for (int i=-static_cast<int>(l.buffer_size); i<static_cast<int>(l.nx+l.buffer_size); ++i) os << std::setw(12) << i;
		os << " x\n";
	}
	os << l.properties;
	return os;
}

void lattice::add_wall(coordinate<int> min_coord, coordinate<int> max_coord)
{
	for (int j=min_coord.j; j<=max_coord.j; ++j)
	{
		for (int i=min_coord.i; i<=max_coord.i; ++i)
		{
			// check if node not yet labelled as wall
			if (!get_node(i,j).has_flag_property("wall"))
			{
				// set wall property
				get_node(i,j).set_flag_property("wall");
				wall_nodes.push_back( get_node(i,j) );
			}
		}
	}
}

void lattice::delete_walls()
{
	for (node n : wall_nodes)
	{
		n.unset_flag_property("wall");
	}
	wall_nodes.clear();
}

void lattice::write_fields(std::string file_name)
{
	std::ofstream ofs(file_name.c_str());
	if (ofs.is_open())
	{
		// write tecplot header (comment that part if necessary)
		ofs << "TITLE=\"LB2D\"\n";
    ofs << "VARIABLES = \"X\", \"Y\", \"RHO\", \"U\", \"V\",\"GAMMA\" \n";
    ofs << "ZONE T = \"BIG ZONE\", I=" << nx << ", J=" << ny << ", F=POINT" << std::endl;
		// write body
		for (unsigned int j=0; j<ny; ++j){
			for (unsigned int i=0; i<nx; ++i){
				ofs << i << " " << j << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].rho() << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].u() << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].v() << " "
				    << std::scientific << nodes[(j+buffer_size)*real_nx + i + buffer_size].gamma() << "\n";
			}
		}
	}
	else throw std::runtime_error("could not write to file");
}

/** WRITE INSTANTANEOUS TIME DATA FOR ANY VARIABLE AT PROBE LOCATION, WRITE FORCE DATA etc**/
void lattice::write_TimeMonitors(unsigned int time,float_type Fx, float_type Fy){
		std::ofstream ofile("output/monitors.txt", std::ios::app);
		float_type X,Y1,Y2,Y3;
		// PROBE LOCATIONS, CAN EDIT LATER
		X  = nx/4 + 2*nx/10;
		Y1 = ny/2. + ny/10;			// P+R
		Y2 = ny/2. - ny/10;			// P-R
		Y3 = ny/2.; 						// P0

		if (ofile.is_open()){
				if (time == 0 ) {
						ofile << "U_P+R" << '\t' << "U_P-R" << '\t' << "U_P0" << '\t' << "X_Force" << '\t' << "Y_Force" << "\n" ;
				}

				ofile << std::scientific << get_node(X,Y1).u() <<'\t'<< get_node(X,Y2).u() <<'\t'<< get_node(X,Y3).u() << '\t' << Fx << '\t' << Fy << "\n";
		}
		else throw std::runtime_error("could not write to file");
		return;
}

} // lb

#endif //LB_LATTICE_HPP_INCLUDED
