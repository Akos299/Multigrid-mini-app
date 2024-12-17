#ifndef __ND_ARRAY__
#define __ND_ARRAY__

//================================================
// @file         ndArray.h
// @author       Jonathan Hadida
// @contact      Jonathan.hadida [at] dtc.ox.ac.uk
//================================================

#include <iostream>
#include <cstdlib>
#include <cstdarg>
#include <cstdint>
#include <cstdio>

#include <exception>
#include <stdexcept>

#include <memory>
#include <algorithm>
#include <type_traits>
#include <initializer_list>

// Comment for (slightly) faster access
#define ND_ARRAY_SAFE_ACCESS
#define ND_ARRAY_USING_MATLAB

// Protect 1D access
#ifdef ND_ARRAY_SAFE_ACCESS
	#define ND_ARRAY_PROTECT(k,n) (k % n)
#else
	#define ND_ARRAY_PROTECT(k,n) k
#endif



        /********************     **********     ********************/
        /********************     **********     ********************/



// // Include Matlab Mex library
// #ifdef ND_ARRAY_USING_MATLAB

// // #include "mex.h"

// /**
//  * Convert numeric types to mex types.
//  */
// template <typename T>
// struct mx_type
// {
// 	static const char *name;
// 	static const mxClassID id;
// };

// // ------------------------------------------------------------------------

// template <> const char* mx_type<char>::name = "int8";
// template <> const char* mx_type<unsigned char>::name = "uint8";
// template <> const char* mx_type<short>::name = "int16";
// template <> const char* mx_type<unsigned short>::name = "uint16";
// template <> const char* mx_type<int>::name = "int32";
// template <> const char* mx_type<unsigned int>::name = "uint32";
// template <> const char* mx_type<float>::name = "single";
// template <> const char* mx_type<double>::name = "double";

// template <> const mxClassID mx_type<char>::id = mxINT8_CLASS;
// template <> const mxClassID mx_type<unsigned char>::id = mxUINT8_CLASS;
// template <> const mxClassID mx_type<short>::id = mxINT16_CLASS;
// template <> const mxClassID mx_type<unsigned short>::id = mxUINT16_CLASS;
// template <> const mxClassID mx_type<int>::id = mxINT32_CLASS;
// template <> const mxClassID mx_type<unsigned int>::id = mxUINT32_CLASS;
// template <> const mxClassID mx_type<float>::id = mxSINGLE_CLASS;
// template <> const mxClassID mx_type<double>::id = mxDOUBLE_CLASS;

// // ------------------------------------------------------------------------

// template <> const char* mx_type<const char>::name = "int8";
// template <> const char* mx_type<const unsigned char>::name = "uint8";
// template <> const char* mx_type<const short>::name = "int16";
// template <> const char* mx_type<const unsigned short>::name = "uint16";
// template <> const char* mx_type<const int>::name = "int32";
// template <> const char* mx_type<const unsigned int>::name = "uint32";
// template <> const char* mx_type<const float>::name = "single";
// template <> const char* mx_type<const double>::name = "double";

// template <> const mxClassID mx_type<const char>::id = mxINT8_CLASS;
// template <> const mxClassID mx_type<const unsigned char>::id = mxUINT8_CLASS;
// template <> const mxClassID mx_type<const short>::id = mxINT16_CLASS;
// template <> const mxClassID mx_type<const unsigned short>::id = mxUINT16_CLASS;
// template <> const mxClassID mx_type<const int>::id = mxINT32_CLASS;
// template <> const mxClassID mx_type<const unsigned int>::id = mxUINT32_CLASS;
// template <> const mxClassID mx_type<const float>::id = mxSINGLE_CLASS;
// template <> const mxClassID mx_type<const double>::id = mxDOUBLE_CLASS;

// #endif



//         /********************     **********     ********************/
//         /********************     **********     ********************/



namespace nd {

// #ifdef ND_ARRAY_USING_MATLAB
// 	using index_t   = mwIndex;
// #else
// 	using index_t   = uint64_t;
// #endif

    using index_t   = uint64_t;
	using dimen_t   = uint8_t;
	using index_ptr = const index_t*;

	/**
	 * Convert nd coordinates to 1d index.
	 * Two main variants are provided:
	 * - Taking an ARRAY as coordinates (size input by template)
	 * - Taking a VA_LIST as a list of coordinate inputs (cf
	 * operator() below).
	 */
	template <dimen_t N>
	index_t sub2ind( index_ptr subs, index_ptr size, index_ptr strides )
	{
		index_t ind = 0;

		for (dimen_t i = 0; i < N; ++i)
			ind += ND_ARRAY_PROTECT(subs[i],size[i]) * strides[i];

		return ind;
	}

	template <dimen_t N>
	index_t sub2ind( va_list& vl, index_ptr size, index_ptr strides )
	{
		index_t ind = 0;

		for (dimen_t i = 1; i < N; ++i)
			ind += ND_ARRAY_PROTECT(va_arg(vl,index_t),size[i]) * strides[i];

		va_end(vl); return ind;
	}

	template <> inline index_t
	sub2ind<0>( index_ptr, index_ptr, index_ptr )
		{ return 0; }
	template <> inline index_t
	sub2ind<1>( index_ptr subs, index_ptr size, index_ptr strides )
		{ return
			ND_ARRAY_PROTECT(subs[0],size[0])*strides[0]; }
	template <> inline index_t
	sub2ind<2>( index_ptr subs, index_ptr size, index_ptr strides )
		{ return
			ND_ARRAY_PROTECT(subs[0],size[0])*strides[0] +
			ND_ARRAY_PROTECT(subs[1],size[1])*strides[1]; }
	template <> inline index_t
	sub2ind<3>( index_ptr subs, index_ptr size, index_ptr strides )
		{ return
			ND_ARRAY_PROTECT(subs[0],size[0])*strides[0] +
			ND_ARRAY_PROTECT(subs[1],size[1])*strides[1] +
			ND_ARRAY_PROTECT(subs[2],size[2])*strides[2]; }

	// ------------------------------------------------------------------------

	/**
	 * Simple singleton.
	 */
	template <typename T> struct singleton { static T instance; };
	template <typename T> T singleton<T>::instance = T();

	// ------------------------------------------------------------------------

	/**
	 * Dummy deleter functor.
	 * This litterally does nothing to the input pointer; it can be used
	 * safely with shared pointers for either statically allocated memory
	 * (eg fixed-size arrays) or externally managed memory (eg Matlab in/out).
	 */
	template <typename T>
	struct no_delete { inline void operator() ( T* ptr ) const {} };



	        /********************     **********     ********************/
	        /********************     **********     ********************/



	/**
	 * n-dimensional array.
	 * NOTE: T can be CONST (underlying elements non-assignable:
	 * suitable for Matlab inputs for instance), or NON-CONST
	 * (underlying elements assignable, suitable for Matlab
	 * outputs or owned memory allocations for instance).
	 */
	template <typename T, dimen_t N>
	class ndArray
	{
	public:

		typedef T value_type;
		typedef T* pointer;
		typedef T& reference;

		typedef typename std::add_const<T>::type const_value;
		typedef const_value* const_pointer;
		typedef const_value& const_reference;

		typedef std::shared_ptr<value_type> shared;
		typedef ndArray<T,N> self;



		// Constructors
		ndArray() { reset(); }  
		ndArray( pointer ptr, index_ptr size, bool manage ) { assign(ptr,size,manage); } /*TESTED*/
        


		// Copy constructor
		ndArray( const self& other ) { operator=(other); } /*TESTED*/
		self& operator= ( const self& other );   

		// Check pointer validity
		inline bool empty() const { return !((bool) m_data); } /*TESTED*/
		inline operator bool() const { return m_data; }

		// Print array dimensions
		void info() const;    /*TESTED*/

        //Arithmeticals operators
       self& operator+=(const self& other); /*TESTED*/
       self& operator+=(const T x); /*TESTED*/
       self& operator*=(const T x); /*TESTED*/
	   self& operator= (const T x);
      


	// #ifdef ND_ARRAY_USING_MATLAB

	// 	// Build from Matlab's mxArray
	// 	ndArray( const mxArray *A ) { assign(A); }
	// 	void assign( const mxArray *A );

	// #endif


		// ------------------------------------------------------------------------


		// Clear contents
		void clear();  /*TESDTED*/
		void reset();  /*TESTED*/


		// Assign either an mxArray or a pointer
		void assign( pointer ptr, index_ptr size, bool manage );


		// Swap contents with another array
		void swap( self& other );  /*TESTED*/


		// Copy from another array
		template <typename U>
		void copy( const ndArray<U,N>& other );  /*TESTED*/


		// ------------------------------------------------------------------------


		// 1D access
		inline reference operator[] ( index_t n ) const
			{ return data()[ ND_ARRAY_PROTECT(n,m_numel) ]; }

		// ND access
		reference operator() ( index_ptr subs ) const
			{ return data()[ sub2ind<N>(subs, m_size, m_strides) ]; }

		reference operator() ( std::initializer_list<index_t> subs ) const
			{
	#ifdef ND_ARRAY_SAFE_ACCESS
				if ( subs.size() != N )
					throw std::length_error("Invalid coordinates length.");
	#endif
				return data()[ sub2ind<N>(subs.begin(), m_size, m_strides) ];
			}

		// Coordinates access
		reference operator() ( index_t i, ... ) const
			{
				va_list vl; va_start(vl,i);
				return data()[ (i*m_strides[0]) + sub2ind<N>(vl, m_size, m_strides) ];
			}


		// Access data directly
		inline const_pointer cdata() const { return m_data.get(); }
		inline pointer data() const { return m_data.get(); }

		// Iterators
		inline const_pointer cbegin() const { return data(); }
		inline const_pointer cend() const { return data() + m_numel; }

		inline pointer begin() const { return data(); }
		inline pointer end() const { return data() + m_numel; }


		// ------------------------------------------------------------------------


		// Dimensions
		inline index_ptr size() const { return m_size; }
		inline index_t size( dimen_t n ) const { return m_size[ ND_ARRAY_PROTECT(n,N) ]; }
		inline index_ptr strides() const { return m_strides; }
		inline index_t stride( dimen_t n ) const { return m_strides[ ND_ARRAY_PROTECT(n,N) ]; }
		inline index_t numel() const { return m_numel; }
		inline index_t ndims() const { return N; }


        void set_zero()      /*TESTED*/
        {
            for(index_t i = 0; i < numel(); i++)
            (*this)({i}) = 0;
        }



	protected:

		void assign_shared( pointer ptr, bool manage );

		index_t m_numel;
		index_t m_size[N];
		index_t m_strides[N];

		shared m_data;
	};

	// ------------------------------------------------------------------------
	// 
	
	/**
 * Assignment operator.
 * Performs a shallow copy of the other instance.
 * For deep-copies, see copy() below.
 */
template <typename T, dimen_t N>
ndArray<T, N> &ndArray<T, N>::operator=(const self &other)
{
	if (other.m_data != m_data)
	{
		// Clear current instance first
		clear();

		// Copy data from other
		m_data = other.m_data;
		m_numel = other.m_numel;
		std::copy_n(other.m_size, N, m_size);
		std::copy_n(other.m_strides, N, m_strides);
	}

	return *this;
}

/**
 * Binary addition operator
 * between two ndArray
 *
 */
template <typename T, dimen_t N>
ndArray<T, N> &ndArray<T, N>::operator+=(const self &other)
{
	for (index_t i = 0; i < numel(); i++)
		(*this)({i}) += other({i});
	return (*this);
}

/**
 * @brief Add a constant to the current instance of ndArray
 */
template <typename T, dimen_t N>
ndArray<T, N> &ndArray<T, N>::operator+=(const T x)
{
	for (index_t i = 0; i < numel(); i++)
		(*this)({i}) += x;
	return (*this);
}

/**
 * @brief Multiply the current instance of ndArray by a constant 
 */
template <typename T, dimen_t N>
ndArray<T, N> &ndArray<T, N>::operator*=(const T x)
{
	for (index_t i = 0; i < numel(); i++)
		(*this)({i}) *= x;
	return (*this);
}

template<typename T, dimen_t N>
ndArray<T,N> &ndArray<T,N>::operator=(const T x)
{
	for(index_t i=0; i < numel(); i++)
		(*this)({i}) = x;
	return (*this);
}
// ------------------------------------------------------------------------

/**
 * Performs a deep-copy of another instance with possibly
 * different value-type. Deep copies are allowed only if
 * the pointer type is non-const (otherwise no assignment
 * is possible).
 *
 * To perform the copy, a new memory allocation is requested
 * to store as many values as other.m_numel; the current
 * instance takes ownership of this new memory.
 *
 * Note that subsequent shallow copies (see assignment operator)
 * will simply share this ownership (reference counting).
 * Note also that this code might generate warnings because of
 * the value-cast performed on the values of other.
 */
template <typename T, dimen_t N>
template <typename U>
void ndArray<T, N>::copy(const ndArray<U, N> &other)
{
	if (!std::is_const<T>::value)
	{
		// Create new allocation only if necessary
		if (other.numel() == m_numel)
		{
			// Otherwise simply copy dimensions
			std::copy_n(other.size(), N, m_size);
			std::copy_n(other.strides(), N, m_strides);
		}
		else
			assign(new T[other.numel()], other.size(), true);

		// Copy data
		auto dst = begin();
		auto src = other.cbegin();
		for (; src != other.cend(); ++src, ++dst)
			*dst = (T)*src;
	}
	else
		throw std::logic_error("Const values cannot be assigned!");
}

// ------------------------------------------------------------------------

/**
 * Reset shared pointer.
 * This will trigger the deletion of the underlying memory if
 * m_data is unique (m_data.use_count() == 1). Note that if the
 * data was assigned with 'manage' set to false (see below),
 * the deleter(no_delete functor) will NOT release the memory.
 */
template <typename T, dimen_t N>
void ndArray<T, N>::clear()
{
	m_data.reset();
}

// ------------------------------------------------------------------------

/**
 * More thorough cleanup. Calls clear() (see above), and sets
 * all the rest to 0.
 */
template <typename T, dimen_t N>
void ndArray<T, N>::reset()
{
	clear();
	m_numel = 0;
	std::fill_n(m_size, N, 0);
	std::fill_n(m_strides, N, 0);
}

// ------------------------------------------------------------------------

/**
 * Swap contents with another ndArray.
 */
template <typename T, dimen_t N>
void ndArray<T, N>::swap(self &other)
{
	std::swap(m_numel, other.m_numel);
	for (dimen_t i = 0; i < N; ++i)
	{
		std::swap(m_size[i], other.m_size[i]);
		std::swap(m_strides[i], other.m_strides[i]);
	}
	m_data.swap(other.m_data);
}

// ------------------------------------------------------------------------

/**
 * Internal method (protected) to assign the shared pointer.
 * Dimensions are assumed to be taken care of by the public
 * assign variants (see below); only the pointer, it's length
 * and the flag 'manage' are required here.
 *
 * 'manage' allows to specify whether or not the shared pointer
 * should release the memory when the last refering instance is
 * destroyed.
 *
 * If true, the default deleter std::default_delete will be
 * assigned to the shared pointer. Note that this deleter releases
 * memory allocated USING NEW ONLY;
 * 	DO NOT use malloc/calloc or other C allocation variants.
 *
 * If false, the deleter no_delete is given instead; this will NOT
 * release the memory when the last refering instance is destroyed.
 * Use only with either externally managed (eg Matlab) or static
 * allocations.
 */
template <typename T, dimen_t N>
void ndArray<T, N>::assign_shared(pointer ptr, bool manage)
{
	if (manage)
		m_data.reset(ptr);
	else
		m_data.reset(ptr, no_delete<T>());
}

// ------------------------------------------------------------------------

// #ifdef ND_ARRAY_USING_MATLAB

// /**
//  * Assign from an mxArray.
//  * This will simply HANDLE the memory allocated by the mxArray,
//  * not manage it; no deep-copy will be performed, no dynamic
//  * allocation will occur, and no deallocation will happen when
//  * this instance (or any shallow copy) is destroyed.
//  *
//  * Use this method to handle Matlab inputs, eg:
//  *
//  * 	ndArray<const double,2> matrix( plhs[0] );
//  *  (note that T is a const type to preserve input data)
//  *
//  * or to handle Matlab outputs, eg:
//  *
//  * 	const int size[5] = {10,11,12,13,14};
//  * 	prhs[0] = mxCreateNumericArray(
//  * 		5, size, mx_type<float>::id, mxREAL );
//  * 	ndArray<float,5> matrix( prhs[0] );
//  * 	(note that the allocated array is assigned to prhs first)
//  *
//  * If the input type or number of dimensions does not correspond
//  * to the template parameters, exceptions are raised.
//  */
// template <typename T, dimen_t N>
// void ndArray<T,N>::assign( const mxArray *A )
// {
// 	// Check dimensions and type
// 	if ( ((dimen_t) mxGetNumberOfDimensions(A)) != N )
// 		throw std::domain_error("Bad dimensions.");
// 	if ( mxGetClassID(A) != mx_type<T>::id )
// 		throw std::invalid_argument("Type mismatch.");

// 	// Get input dimensions
// 	index_ptr size = (index_ptr) mxGetDimensions(A);

// 	// Call assign variant
// 	assign( (pointer) mxGetData(A), size, false );
// }

// #endif

// ------------------------------------------------------------------------

/**
 * Assign from pointer and size.
 *
 * If manage == true, the internal shared pointer m_data
 * will assume ownership of the memory pointed by ptr, and
 * try to release it using delete[] when the last refering
 * instance gets destroyed.
 * If ptr is dynamically allocated, make sure that the
 * allocation is performed with NEW, and NOT C variants
 * like malloc/calloc/etc.
 *
 * If manage == false, a dummy deleter (no_delete functor) is
 * passed to the shared pointer; nothing happens to the memory
 * pointed by ptr when the last refering instance gets destroyed.
 */
template <typename T, dimen_t N>
void ndArray<T, N>::assign(pointer ptr, index_ptr size, bool manage)
{
	if (ptr != data())
	{
		// Compute internal dimensions
		m_numel = 1;
		for (dimen_t i = 0; i < N; ++i)
		{
			m_size[i] = size[i];
			m_numel *= size[i];
			m_strides[(i + 1) % N] = m_numel;
		}
		m_strides[0] = 1;

		// Set data
		assign_shared(ptr, manage);
	}
}

// ------------------------------------------------------------------------

/**
 * Simply prints information about the dimensions of the
 * n-dimensional array (size & number of elements).
 */
template <typename T, dimen_t N>
void ndArray<T, N>::info() const
{
	if (m_data)
	{
		printf("%u-dimensional array of size (%lu", N, m_size[0]);
		for (dimen_t d = 1; d < N; ++d)
			printf(", %u", m_size[d]);
		printf(") = %u elements.\n", m_numel);
	}
	else
		printf("Empty %u-dimensional array.\n", N);
}
}

#endif