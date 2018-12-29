/** 
 *  @file
 *  @author Fabian Bösch
 *  @brief property array
 */
 
#ifndef LB_PROPERTY_ARRAY_HPP_INCLUDED
#define LB_PROPERTY_ARRAY_HPP_INCLUDED

#include <string>
#include <vector>
#include <typeinfo>
#include <typeindex>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <stdexcept>

namespace lb {

namespace detail {
	
template<typename SizeType, typename UIntType>
struct flag_property_description
{
	typedef SizeType size_type;
	typedef UIntType uint_type;
	
	flag_property_description(std::string name_, size_type index_) 
	: name(name_), index(index_), mask(static_cast<uint_type>(1) << (index)) {}
	flag_property_description(const flag_property_description&) = default;
	flag_property_description& operator=(const flag_property_description&) = default;
	flag_property_description() = default;
	
	std::string name;
	size_type index;
	uint_type mask;
};

template<typename SizeType, typename UIntType>
struct data_property_description
{
	typedef SizeType size_type;
	typedef UIntType uint_type;
	
	data_property_description(std::string name_, std::type_index t_index_, size_type index_) 
	: name(name_), t_index(t_index_), index(index_) {}
	data_property_description(const data_property_description&) = default;
	data_property_description& operator=(const data_property_description&) = default;
	virtual ~data_property_description() = default;
	
	virtual void release_property(void*) const {}
	
	virtual void* clone_property(void*) const { return nullptr; }
	
	virtual data_property_description* clone() const { return nullptr; }
	
	std::string name;
	std::type_index t_index;
	size_type index;
};

template<typename SizeType, typename UIntType, typename T>
struct data_property_description_impl : public data_property_description<SizeType, UIntType>
{
	typedef SizeType size_type;
	typedef UIntType uint_type;
	
	data_property_description_impl(std::string name, size_type index) 
	: data_property_description<SizeType, UIntType>(name, std::type_index(typeid(T)), index) {}
	
	data_property_description_impl(const data_property_description_impl& other) : data_property_description<SizeType, UIntType>(other) {}
	
	data_property_description_impl& operator=(const data_property_description_impl& other)
	{
		data_property_description<SizeType, UIntType>::operator=(other);
		return *this;
	}
	
	virtual ~data_property_description_impl() = default;
	
	virtual void release_property(void* ptr) const { delete reinterpret_cast<T*>(ptr); }
	
	virtual void* clone_property(void* ptr) const { return new T(*(reinterpret_cast<T*>(ptr))); }
	
	virtual data_property_description_impl* clone() const { return new data_property_description_impl(*this); }
};
	
struct find_flag_property_by_name
{
	find_flag_property_by_name(const std::string& name_) : name(name_) {}
	template <typename FlagProperty>
	inline bool operator()(const FlagProperty& p) const { return (p.name == name); }
	const std::string& name;
};

struct find_data_property_by_name
{
	find_data_property_by_name(const std::string& name_) : name(name_) {}
	template <typename DataProperty>
	inline bool operator()( const DataProperty * p ) const { return (p->name == name); }
	const std::string& name;
};
	
} // detail

/**
 *  @brief This class allows you to store properties in an array.
 * 
 *  For every index in the property array you can register several 
 *  flag properties (boolean values) and several custom objects of any
 *  type. As long as you do not store a data property for a given index
 *  no space will be used.
 * 
 *  Example usage for flag property:
 *  @code
 *  	property_array pa(5);                                            // property array of size 5
 *  	pa.register_register_flag_property("test_flag",false);           // register a flag with name "test_flag", default value false
 * 
 *  	pa.set_flag_property("test_flag",2);                             // set the flag at position 2 to true
 *  	// alternatively you can use the property index
 *  	property_array::size_type idx;
 *  	pa.flag_property_index("test_flag", idx);
 *  	pa.set_flag_property(idx,2);                                     // set the flag at position 2 to true
 * 
 *  	bool value = pa.has_flag_property("test_flag",2);                // get the value at position 2 for your flag
 *  	// using the index
 *  	value = pa.has_flag_property(idx,2);                             // get the value at position 2 for your flag
 * 
 *  	pa.unset_flag_property("test_flag",2);                           // set flag at position 2 to false
 *  	// using the index
 *  	pa.unset_flag_property(idx,2);                                   // set flag at position 2 to false
 *  @endcode
 * 
 *  Example usage for data property:
 *  @code
 *  	property_array pa(5);                                                           // property array of size 5
 *  	// assume you want to store a std::vector<float>
 *  	pa.register_register_data_property<std::vector<float> >("some_data",false);     // register a data property with name "some_data", default value: null
 * 
 *  	std::vector<float> test_property(9,1.0);
 *  	pa.set_data_property("some_data", 3, test_property);                            // store data property "some_data" at position 3 to be equal to variable test_property
 *  	// again you could use an index instead of a string to access the property
 *  	property_array::size_type idx;
 *  	pa.data_property_index("some_data", idx);
 *  	pa.set_data_property("some_data", 3, test_property);
 *  	
 *  	
 *  	if (pa.has_data_property("some_data", 3))                                       // check whether there is a data property stored at index 3
 *  	{
 *  		std::vector<float> return_value;
 *  		return_value=pa.get_data_property<std::vector<float> >("some_data",3);  // retrieve the data property
 *  	}
 *  	
 *  	pa.unset_data_property("some_data", 3);                                         // delete the data property at index 3
 *  @endcode
 */
class property_array
{
public: // typedefs

	/** @brief enumeration and size type */
	typedef unsigned long int size_type;
	
private: // typedefs

	typedef unsigned long int uint_type;
	typedef detail::flag_property_description<size_type,uint_type> flag_property_description_type;
	typedef detail::data_property_description<size_type,uint_type> data_property_description_type;
	template <typename T> using data_property_description_impl_type = detail::data_property_description_impl<size_type,uint_type, T>;
	
	typedef std::vector<flag_property_description_type> flag_description_vector_type;
	typedef std::vector<data_property_description_type*> data_description_ptr_vector_type;
	typedef std::vector<uint_type> set_vector_type;
	typedef std::vector<std::vector<void*> > data_property_storage_vector_type;

private: // members	

	size_type size_;
	set_vector_type flag_set_vector;
	flag_description_vector_type flag_description_vector;
	data_description_ptr_vector_type data_description_ptr_vector;
	data_property_storage_vector_type data_property_storage_vector;
	
public: // ctors

	/**
	 *  @brief Construct with array size
	 *  @param[in] size array size
	 */
	property_array(size_type size) : size_(size), flag_set_vector(size_,0) {}
	
	/** @brief Copy construct */
	property_array(const property_array& other)
	: size_(other.size_), flag_set_vector(other.flag_set_vector), 
	  flag_description_vector(other.flag_description_vector),
	  data_description_ptr_vector(other.data_description_ptr_vector), 
	  data_property_storage_vector(other.data_property_storage_vector.size(), std::vector<void*>(size_,0))
	{
		for (size_type k=0; k<num_data_properties(); ++k)
		{
			data_description_ptr_vector[k] = data_description_ptr_vector[k]->clone();
			for (size_type i=0; i<size_; ++i)
			{
				if (other.data_property_storage_vector[k][i]) 
					data_property_storage_vector[k][i] = data_description_ptr_vector[k]->clone_property(other.data_property_storage_vector[k][i]);
			}
		}
	}
	
	/** @brief Swap internal state */
	void swap(property_array& other)
	{
		std::swap(size_, other.size_);
		std::swap(flag_set_vector, other.flag_set_vector);
		std::swap(flag_description_vector, other.flag_description_vector);
		std::swap(data_description_ptr_vector, other.data_description_ptr_vector);
		std::swap(data_property_storage_vector, other.data_property_storage_vector);
	}
	
	/** @brief Assignement */
	property_array& operator=(property_array other)
	{
		this->swap(other);
		return *this;
	}
	
	property_array(property_array&&) = default;
	
	/** Destruct */
	~property_array()
	{
		for (size_type k=0; k<num_data_properties(); ++k)
		{
			for (size_type i=0; i<size_; ++i)
			{
				if (data_property_storage_vector[k][i]) data_description_ptr_vector[k]->release_property(data_property_storage_vector[k][i]);
			}
			delete data_description_ptr_vector[k];
		}
	}

public: // register property types

	/** 
	 *  @brief register a new flag property
	 *  @param[in] name flag property name
	 *  @param[in] set  default state
	 *  @return true if property does not exist yet
	 */
	bool register_flag_property(std::string name, bool set = false)
	{
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) return false;
		flag_description_vector.push_back(flag_property_description_type(name, flag_description_vector.size()));
		if (set) for (uint_type& i : flag_set_vector) i |= flag_description_vector.back().mask;
		return true;
	}
	
	/** 
	 *  @brief register a new data property
	 *  @tparam T type for data property
	 *  @param[in] name  data property name
	 *  @param[in] set   default behaviour: initialized or null
	 *  @param[in] value value for default initialization
	 *  @return true if property does not exist yet
	 */
	template <typename T>
	bool register_data_property(std::string name, bool set = false, const T& value = T())
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name,d_property_desc_ptr) || (typeid(void) == typeid(T))) return false;
		data_description_ptr_vector.push_back(new data_property_description_impl_type<T>(name, data_description_ptr_vector.size()) );
		data_property_storage_vector.push_back( std::vector<void*>(size_,0) );
		if (set) for (size_type i=0; i<size_; ++i) data_property_storage_vector[data_property_storage_vector.size()-1][i] = new T(value);
		return true;
	}
	
public: // set and unset properties

	/** 
	 *  @brief set flag property
	 *  @param[in] name        flag property name
	 *  @param[in] array_index position in array
	 *  @return true if flag property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool set_flag_property(std::string name, size_type array_index)
	{
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) { flag_set_vector[array_index] |= f_property_desc_ptr->mask; return true; }
		else return false;
	}
	
	/** 
	 *  @brief set data property
	 *  @tparam T type for data property
	 *  @param[in] name        data property name
	 *  @param[in] array_index position in array
	 *  @param[in] property    data to store
	 *  @return true if data property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template<typename T>
	bool set_data_property(std::string name, size_type array_index, const T& property)
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name,d_property_desc_ptr)) 
		{ 
			if (data_property_storage_vector[d_property_desc_ptr->index][array_index] != 0)
			{
				d_property_desc_ptr->release_property(data_property_storage_vector[d_property_desc_ptr->index][array_index]);
			}
			data_property_storage_vector[d_property_desc_ptr->index][array_index] = new T(property);
			return true;
		}
		else return false;
	}
	
	/** 
	 *  @brief set flag property
	 *  @param[in] property_index  flag property index
	 *  @param[in] array_index     position in array
	 *  @return true if flag property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool set_flag_property(size_type property_index, size_type array_index)
	{
		if (property_index < num_flag_properties()) { flag_set_vector[array_index] |= flag_description_vector[property_index].mask; return true; }
		else return false;
	}
	
	/** 
	 *  @brief set data property
	 *  @tparam T type for data property
	 *  @param[in] property_index data property index
	 *  @param[in] array_index    position in array
	 *  @param[in] property       data to store
	 *  @return true if data property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template<typename T>
	bool set_data_property(size_type property_index, size_type array_index, const T& property)
	{
		if (property_index < num_data_properties()) 
		{ 
			if (data_property_storage_vector[property_index][array_index] != 0)
			{
				data_description_ptr_vector[property_index]->release_property(data_property_storage_vector[property_index][array_index]);
			}
			data_property_storage_vector[property_index][array_index] = new T(property);
			return true;
		}
		else return false;
	}
	
	/** 
	 *  @brief unset flag property
	 *  @param[in] name        flag property name
	 *  @param[in] array_index position in array
	 *  @return true if flag property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool unset_flag_property(std::string name, size_type array_index)
	{
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) { flag_set_vector[array_index] &= ~(f_property_desc_ptr->mask); return true; }
		else return false;
	}
	
	/** 
	 *  @brief unset data property
	 *  @param[in] name        data property name
	 *  @param[in] array_index position in array
	 *  @return true if data property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool unset_data_property(std::string name, size_type array_index)
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name,d_property_desc_ptr)) 
		{ 
			if (data_property_storage_vector[d_property_desc_ptr->index][array_index] != 0)
			{
				d_property_desc_ptr->release_property(data_property_storage_vector[d_property_desc_ptr->index][array_index]);
				data_property_storage_vector[d_property_desc_ptr->index][array_index] = 0;
			}
			return true;
		}
		else return false;
	}
	
	/** 
	 *  @brief unset flag property
	 *  @param[in] property_index  flag property index
	 *  @param[in] array_index     position in array
	 *  @return true if flag property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool unset_flag_property(size_type property_index, size_type array_index)
	{
		if (property_index < num_flag_properties()) { flag_set_vector[array_index] &= ~(flag_description_vector[property_index].mask); return true; }
		else return false;
	}
	
	/** 
	 *  @brief unset data property
	 *  @param[in] property_index  data property index
	 *  @param[in] array_index     position in array
	 *  @return true if data property exists
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool unset_data_property(size_type property_index, size_type array_index)
	{
		if (property_index < num_data_properties()) 
		{ 
			if (data_property_storage_vector[property_index][array_index] != 0)
			{
				data_description_ptr_vector[property_index]->release_property(data_property_storage_vector[property_index][array_index]);
				data_property_storage_vector[property_index][array_index] = 0;
			}
			return true;
		}
		else return false;
	}

public: // element access for data properties

	/**
	 *  @brief access data property object
	 *  @tparam T data property type
	 *  @param[in] name        data property name
	 *  @param[in] array_index position in array
	 *  @return reference to data property object
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template <typename T>
	T& get_data_property(std::string name, size_type array_index)
	{
		data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name, d_property_desc_ptr) && 
		    (d_property_desc_ptr->t_index == std::type_index(typeid(T))) &&
		    (data_property_storage_vector[d_property_desc_ptr->index][array_index] != 0))
		{
			return *(reinterpret_cast<T*>(data_property_storage_vector[d_property_desc_ptr->index][array_index]));
		}
		else throw std::range_error("invalid acces to data property!");
	}
	
	/**
	 *  @brief access data property object
	 *  @tparam T data property type
	 *  @param[in] name        data property name
	 *  @param[in] array_index position in array
	 *  @return const reference to data property object
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template <typename T>
	const T& get_data_property(std::string name, size_type array_index) const
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name, d_property_desc_ptr) && 
		    (d_property_desc_ptr->t_index == std::type_index(typeid(T))) &&
		    (data_property_storage_vector[d_property_desc_ptr->index][array_index] != 0))
		{
			return *(reinterpret_cast<const T*>(data_property_storage_vector[d_property_desc_ptr->index][array_index]));
		}
		else throw std::range_error("invalid acces to data property!");
	}
	
	/**
	 *  @brief access data property object
	 *  @tparam T data property type
	 *  @param[in] property_index data property index
	 *  @param[in] array_index    position in array
	 *  @return reference to data property object
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template <typename T>
	T& get_data_property(size_type property_index, size_type array_index)
	{
		data_property_description_type * d_property_desc_ptr;
		if ((property_index < num_data_properties()) && 
		    (data_description_ptr_vector[property_index]->t_index == std::type_index(typeid(T))) &&
		    (data_property_storage_vector[property_index][array_index] != 0))
		{
			return *(reinterpret_cast<T*>(data_property_storage_vector[property_index][array_index]));
		}
		else throw std::range_error("invalid acces to data property!");
	}
	
	/**
	 *  @brief access data property object
	 *  @tparam T data property type
	 *  @param[in] property_index data property index
	 *  @param[in] array_index    position in array
	 *  @return const reference to data property object
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	template <typename T>
	const T& get_data_property(size_type property_index, size_type array_index) const
	{
		const data_property_description_type * d_property_desc_ptr;
		if ((property_index < num_data_properties()) && 
		    (data_description_ptr_vector[property_index]->t_index == std::type_index(typeid(T))) &&
		    (data_property_storage_vector[property_index][array_index] != 0))
		{
			return *(reinterpret_cast<const T*>(data_property_storage_vector[property_index][array_index]));
		}
		else throw std::range_error("invalid acces to data property!");
	}

public: // name - property_index mapping
	
	/**
	 *  @brief get flag property index from flag property name
	 *  @param[in]  name           flag property name
	 *  @param[out] property_index flag property index
	 *  @return true if property exists
	 */
	bool flag_property_index(std::string name, size_type& property_index) const
	{
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) { property_index = f_property_desc_ptr->index; return true; }
		else return false;
	}
	
	/**
	 *  @brief get data property index from data property name
	 *  @param[in]  name           data property name
	 *  @param[out] property_index data property index
	 *  @return true if property exists
	 */
	bool data_property_index(std::string name, size_type& property_index) const
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name,d_property_desc_ptr)) { property_index = d_property_desc_ptr->index; return true; }
		else return false;
	}
	
	/**
	 *  @brief get flag property name from flag property index
	 *  @param[in]  property_index flag property index
	 *  @param[out] name           flag property name
	 *  @return true if property exists
	 */
	bool flag_property_name(size_type property_index, std::string& name) const
	{
		if (property_index < num_flag_properties()) { name = flag_description_vector[property_index].name; return true; }
		else return false;
	}
	
	/**
	 *  @brief get data property name from data property index
	 *  @param[in]  property_index data property index
	 *  @param[out] name           data property name
	 *  @return true if property exists
	 */
	bool data_property_name(size_type property_index, std::string& name) const
	{
		if (property_index < num_data_properties()) { name = data_description_ptr_vector[property_index]->name; return true; }
		else return false;
	}
	
public: // queries

	/** @brief property array size */
	size_type size() const { return size_; }

	/** @brief number of registered flag properties */
	size_type num_flag_properties() const { return flag_description_vector.size(); }
	
	/** @brief number of registered data properties */
	size_type num_data_properties() const { return data_description_ptr_vector.size(); }
	
	/** @brief number of registered flag and data properties */
	size_type num_properties() const { return num_flag_properties() + num_data_properties(); }
	
	/**
	 *  @brief Check whether flag property is set
	 *  @param[in] name        flag property name
	 *  @param[in] array_index position in array
	 *  @result true if property is set
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool has_flag_property(std::string name, size_type array_index) const 
	{ 
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) return (flag_set_vector[array_index] & f_property_desc_ptr->mask);
		else return false;
	}
	
	/**
	 *  @brief Check whether data property is set
	 *  @param[in] name        data property name
	 *  @param[in] array_index position in array
	 *  @result true if property is set
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool has_data_property(std::string name, size_type array_index) const
	{ 
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name, d_property_desc_ptr)) return ( data_property_storage_vector[d_property_desc_ptr->index][array_index] != 0 );
		else return false;
	}
	
	/**
	 *  @brief Check whether flag property is set
	 *  @param[in] property_index flag property index
	 *  @param[in] array_index    position in array
	 *  @result true if property is set
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool has_flag_property(size_type property_index, size_type array_index) const 
	{ 
		if (property_index < num_flag_properties()) return ( flag_set_vector[array_index] & flag_description_vector[property_index].mask );
		else return false;
	}
	
	/**
	 *  @brief Check whether data property is set
	 *  @param[in] property_index data property index
	 *  @param[in] array_index    position in array
	 *  @result true if property is set
	 *  @pre array_index is in range [ 0, @ref property_array::size() )
	 */
	bool has_data_property(size_type property_index, size_type array_index) const
	{ 
		if (property_index < num_data_properties()) return ( data_property_storage_vector[property_index][array_index] != 0 );
		else return false;
	}
	
	/**
	 *  @brief count flag properties which are set
	 *  @param[in]  name           flag property index
	 *  @param[out] count          number of flag properties set
	 *  @return true if flag property exists
	 */
	bool count_set_flag_properties(std::string name, size_type& count) const 
	{
		const flag_property_description_type* f_property_desc_ptr;
		if (find_flag_property(name,f_property_desc_ptr)) 
		{
			count = 0;
			for (size_type i=0; i<size_; ++i) if (flag_set_vector[i] & f_property_desc_ptr->mask) ++count;
			return true;
		}
		else return false;
	}
	
	/**
	 *  @brief count data properties which are set
	 *  @param[in]  name           data property name
	 *  @param[out] count          number of data properties set
	 *  @return true if data property exists
	 */
	bool count_set_data_properties(std::string name, size_type& count) const 
	{
		const data_property_description_type * d_property_desc_ptr;
		if (find_data_property(name, d_property_desc_ptr))
		{ 
			count = 0;
			for (size_type i=0; i<size_; ++i) if (data_property_storage_vector[d_property_desc_ptr->index][i] != 0) ++count;
			return true;
		}
		else return false;
	}
	
	/**
	 *  @brief count flag properties which are set
	 *  @param[in]  property_index flag property index
	 *  @param[out] count          number of flag properties set
	 *  @return true if flag property exists
	 */
	bool count_set_flag_properties(size_type property_index, size_type& count) const 
	{
		if (property_index < num_flag_properties()) 
		{
			count = 0;
			for (size_type i=0; i<size_; ++i) if (flag_set_vector[i] & flag_description_vector[property_index].mask) ++count;
			return true;
		}
		else return false;
	}
	
	/**
	 *  @brief count data properties which are set
	 *  @param[in]  property_index data property index
	 *  @param[out] count          number of data properties set
	 *  @return true if data property exists
	 */
	bool count_set_data_properties(size_type property_index, size_type& count) const 
	{
		if (property_index < num_data_properties())
		{ 
			count = 0;
			for (size_type i=0; i<size_; ++i) if (data_property_storage_vector[property_index][i] != 0) ++count;
			return true;
		}
		else return false;
	}
	
	/**
	 *  @brief check whether data property exists
	 *  @param[in] name data property name
	 *  @return true if exists
	 */
	bool exist_data_property(std::string name) const 
	{
		const data_property_description_type * d_property_desc_ptr;
		return (find_data_property(name, d_property_desc_ptr));
	}
	
	/**
	 *  @brief check whether flag property exists
	 *  @param[in] name flag property name
	 *  @return true if exists
	 */
	bool exist_flag_property(std::string name) const 
	{
		const flag_property_description_type * f_property_desc_ptr;
		return (find_flag_property(name, f_property_desc_ptr));
	}
	
public: // print

	/** @brief print to output stream */
	friend std::ostream& operator<<(std::ostream& os, const property_array& pa)
	{
		os << "property array of size " << pa.size() << " with " << pa.num_properties() << " properties (" 
		   << pa.num_flag_properties() << " flag and " << pa.num_data_properties() << " data properties)\n";
		os << "  flag properties\n";
		typename property_array::size_type i=0;
		for (const auto& f_desc : pa.flag_description_vector)
		{
			size_type count;
			pa.count_set_flag_properties(i,count);
			os << "    property " << std::setw(3) << std::setfill(' ') << i << ": " << f_desc.name << " (" << count << " are set)\n";
			++i;
		}
		os << "  data properties\n";
		i=0;
		for (const auto d_desc_ptr : pa.data_description_ptr_vector)
		{
			size_type count;
			pa.count_set_data_properties(i, count);
			os << "    property " << std::setw(3) << std::setfill(' ') << i << ": " << d_desc_ptr->name << " with type " << d_desc_ptr->t_index.name() << " (" << count << " are set)\n";
			++i;
		}
		return os;
	}
	
private: // impl

	inline bool find_flag_property(std::string name, const flag_property_description_type * & f_property_desc_ptr) const 
	{
		flag_description_vector_type::const_iterator iter = 
			std::find_if(flag_description_vector.begin(), flag_description_vector.end(), detail::find_flag_property_by_name(name));
		if (iter != flag_description_vector.end()) { f_property_desc_ptr = &(*iter); return true; }
		else { f_property_desc_ptr = 0; return false; }
	}
	
	inline bool find_data_property(std::string name, const data_property_description_type * & d_property_desc_ptr) const 
	{
		data_description_ptr_vector_type::const_iterator iter = 
			std::find_if(data_description_ptr_vector.begin(), data_description_ptr_vector.end(), detail::find_data_property_by_name(name));
		if (iter != data_description_ptr_vector.end()) { d_property_desc_ptr = (*iter); return true; }
		else { d_property_desc_ptr = 0; return false; }
	}
	
	inline bool find_flag_property(std::string name, flag_property_description_type * & f_property_desc_ptr) 
	{
		flag_description_vector_type::iterator iter = 
			std::find_if(flag_description_vector.begin(), flag_description_vector.end(), detail::find_flag_property_by_name(name));
		if (iter != flag_description_vector.end()) { f_property_desc_ptr = &(*iter); return true; }
		else { f_property_desc_ptr = 0; return false; }
	}
	
	inline bool find_data_property(std::string name, data_property_description_type * & d_property_desc_ptr) 
	{
		data_description_ptr_vector_type::iterator iter = 
			std::find_if(data_description_ptr_vector.begin(), data_description_ptr_vector.end(), detail::find_data_property_by_name(name));
		if (iter != data_description_ptr_vector.end()) { d_property_desc_ptr = (*iter); return true; }
		else { d_property_desc_ptr = 0; return false; }
	}

};

} // lb


#endif // LB_PROPERTY_ARRAY_HPP_INCLUDED
