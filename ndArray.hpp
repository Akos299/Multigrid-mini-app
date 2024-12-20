// using namespace nd;
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