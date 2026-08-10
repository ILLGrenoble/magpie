/**
 * tlibs2 maths library -- containers and adapters
 * @author Tobias Weber <tobias.weber@tum.de>, <tweber@ill.fr>
 * @date 2015 - 2026
 * @license GPLv3, see 'LICENSE' file
 *
 * @note this file is based on code from my following projects:
 *         - "mathlibs" (https://github.com/t-weber/mathlibs),
 *         - "geo" (https://github.com/t-weber/geo),
 *         - "misc" (https://github.com/t-weber/misc).
 *         - "magtools" (https://github.com/t-weber/magtools).
 *         - "tlibs" (https://github.com/t-weber/tlibs).
 *
 * @desc for the references, see the 'LITERATURE' file
 *
 * ----------------------------------------------------------------------------
 * tlibs2
 * Copyright (C) 2017-2026  Tobias WEBER (Institut Laue-Langevin (ILL),
 *                          Grenoble, France).
 * tlibs1
 * Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
 *                          (TUM), Garching, Germany).
 * "magtools", "geo", "misc", and "mathlibs" projects
 * Copyright (C) 2017-2022  Tobias WEBER (privately developed).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ----------------------------------------------------------------------------
 */

#ifndef __TLIBS2_MATHS_CONTS_H__
#define __TLIBS2_MATHS_CONTS_H__

#include <cmath>
#include <vector>

#include "decls.h"
#include "operators.h"


#ifndef __TLIBS2_DEFAULT_ALLOC__
	//#define __TLIBS2_DEFAULT_ALLOC__ alloc_noinit
	#define __TLIBS2_DEFAULT_ALLOC__ std::allocator
#endif



namespace tl2 {
// ----------------------------------------------------------------------------
// adapters
// ----------------------------------------------------------------------------

/**
 * vector-like access adapter to a matrix
 */
template<class t_mat> requires is_basic_mat<t_mat>
class matvec_adapter
{
public:
	using value_type = typename t_mat::value_type;
	using size_type = decltype(t_mat{}.size1());

public:
	matvec_adapter(const t_mat &mat) : m_mat{ mat }
	{}
	~matvec_adapter() = default;

	size_type size() const { return m_mat.size1() * m_mat.size2(); }

	const value_type& operator[](size_type i) const
	{
		size_type row = i/m_mat.size2();
		size_type col = i%m_mat.size2();

		return m_mat(row, col);
	}

private:
	const t_mat& m_mat;
};


/**
 * adapter for a qvector
 */
template<typename size_t, size_t N, typename T,
	template<size_t, size_t, class...> class t_mat_base>
class qvec_adapter : public t_mat_base<1, N, T>
{
public:
	// types
	using base_type = t_mat_base<1, N, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvec_adapter(const base_type& vec) : base_type{ vec }
	{}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i)
	{
		return base_type::operator()(i,0);
	}

	const T operator[](size_t i) const
	{
		return base_type::operator()(i,0);
	}
};


/**
 * adapter for a qmatrix
 */
template<typename size_t, size_t ROWS, size_t COLS, typename T,
	template<size_t, size_t, class...> class t_mat_base>
class qmat_adapter : public t_mat_base<COLS, ROWS, T>
{
public:
	// types
	using base_type = t_mat_base<COLS, ROWS, T>;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmat_adapter(const base_type& mat) : base_type{ mat }
	{}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};


/**
 * adapter for a qvector
 */
template<typename size_t, size_t N, typename T, class t_vec_base>
class qvecN_adapter : public t_vec_base
{
public:
	// types
	using base_type = t_vec_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qvecN_adapter(const base_type& vec) : base_type{ vec }
	{}

	static constexpr size_t size() { return N; }

	T& operator[](size_t i)
	{
		return static_cast<base_type&>(*this)[i];
	}

	const T operator[](size_t i) const
	{
		return static_cast<const base_type&>(*this)[i];
	}
};


/**
 * adapter for a qmatrix
 */
template<typename size_t, size_t ROWS, size_t COLS, typename T, class t_mat_base>
class qmatNN_adapter : public t_mat_base
{
public:
	// types
	using base_type = t_mat_base;
	using size_type = size_t;
	using value_type = T;

	// constructors
	using base_type::base_type;
	qmatNN_adapter(const base_type& mat) : base_type{ mat }
	{}

	// convert from a different matrix type
	template<class t_matOther> qmatNN_adapter(const t_matOther& matOther)
		requires is_basic_mat<t_matOther>
	{
		const std::size_t minRows = std::min(
			static_cast<std::size_t>(size1()),
			static_cast<std::size_t>(matOther.size1()));
		const std::size_t minCols = std::min(
			static_cast<std::size_t>(size2()),
			static_cast<std::size_t>(matOther.size2()));

		for(std::size_t i = 0; i < minRows; ++i)
			for(std::size_t j = 0; j < minCols; ++j)
				(*this)(i, j) = static_cast<value_type>(matOther(i, j));
	}

	static constexpr size_t size1() { return ROWS; }
	static constexpr size_t size2() { return COLS; }
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// allocator
// ----------------------------------------------------------------------------
/**
 * allocator that doesn't initialise or construct its elements
 * @see https://en.cppreference.com/cpp/memory/allocator
 *
 * WARNING: without initialisation, functions like push_back also stop working
 */
template<class t_elem>
struct alloc_noinit
{
	static constexpr const bool init_elem = false;  // initialise elements?

	using value_type = t_elem;
	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;

	using propagate_on_container_move_assignment = std::true_type;
	using is_always_equal = std::true_type;


	constexpr alloc_noinit() noexcept = default;
	constexpr alloc_noinit(const alloc_noinit<t_elem>&) noexcept = default;

	constexpr ~alloc_noinit() = default;


	constexpr t_elem* allocate(std::size_t num)
	{
		if(num == 0)
			return nullptr;

		return reinterpret_cast<t_elem*>(std::malloc(sizeof(t_elem)*num));
	}


	void deallocate(t_elem *mem, [[__maybe_unused__]] std::size_t num)
	{
		if(!mem)
			return;

		std::free(reinterpret_cast<void*>(mem));
	}


	template<class t_other_elem, class ...t_args>
	void construct(t_other_elem *mem, t_args&& ...args)
	{
		if constexpr(init_elem)
		{
			if(!mem)
				return;

			// call constructor
			void *pmem = reinterpret_cast<void*>(mem);
			new(pmem) t_other_elem{ std::forward<t_args>(args)... };
		}
	}


	template<class t_other_elem>
	void destroy(t_other_elem *t)
	{
		if constexpr(init_elem)
		{
			if(!t)
				return;

			// call destructor
			t->~t_other_elem();
		}
	}
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// vector and matrix containers
// ----------------------------------------------------------------------------

/**
 * generic vector container
 */
template<class T = double,
	template<class...> class t_cont = std::vector,
	template<class...> class t_alloc = __TLIBS2_DEFAULT_ALLOC__>
requires is_basic_vec<t_cont<T, t_alloc<T>>> && is_dyn_vec<t_cont<T, t_alloc<T>>>
class vec : public t_cont<T, t_alloc<T>>
{
public:
	using value_type = T;
	using container_type = t_cont<T, t_alloc<T>>;

	using container_type::container_type;
	using container_type::size;
	using container_type::operator[];

	using typename container_type::iterator;
	using typename container_type::const_iterator;
	using typename container_type::size_type;
	using typename container_type::difference_type;
	using typename container_type::allocator_type;


/*private:
	using container_type::push_back;
	using container_type::emplace_back;*/


public:
	~vec() = default;

	vec() : container_type{}
	{}

	vec(const vec<T, t_cont, t_alloc>& other) //: container_type{ other }
	{
		this->resize(other.size());
		for(std::size_t i = 0; i < other.size(); ++i)
			(*this)[i] = other[i];
	}

	vec(vec<T, t_cont, t_alloc>&& other) noexcept : container_type{ std::move(other) }
	{}

	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other>
	vec(const vec<T_other, t_cont_other, t_alloc_other>& other)
	{
		this->operator=<T_other, t_cont_other, t_alloc_other>(other);
	}

	vec(const std::initializer_list<T>& lst)
	{
		this->resize(lst.size());
		size_type i = 0;
		for(const T& elem : lst)
		{
			(*this)[i] = elem;
			++i;
		}
	}


	vec<T, t_cont, t_alloc>& operator=(const vec<T, t_cont, t_alloc>& other)
	{
		//*static_cast<container_type*>(this) = other;
		this->resize(other.size());
		for(std::size_t i = 0; i < other.size(); ++i)
			(*this)[i] = other[i];

		return *this;
	}

	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other>
	vec<T, t_cont, t_alloc>& operator=(const vec<T_other, t_cont_other, t_alloc_other>& other)
	{
		*this = convert<vec<T, t_cont, t_alloc>, vec<T_other, t_cont_other, t_alloc_other>>(other);
		return *this;
	}


	vec(std::size_t SIZE, const T* arr = nullptr) : container_type(SIZE)
	{
		if(arr)
			from_array(arr);
	}


	const value_type& operator()(std::size_t i) const
	{
		return this->operator[](i);
	}

	value_type& operator()(std::size_t i)
	{
		return this->operator[](i);
	}


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i = 0; i < size(); ++i)
			this->operator[](i) = arr[i];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i = 0; i < size(); ++i)
			arr[i] = this->operator[](i);
	}


	friend vec operator+(const vec& vec1, const vec& vec2) { return tl2_ops::operator+(vec1, vec2); }
	friend vec operator-(const vec& vec1, const vec& vec2) { return tl2_ops::operator-(vec1, vec2); }
	friend const vec& operator+(const vec& vec1) { return tl2_ops::operator+(vec1); }
	friend vec operator-(const vec& vec1) { return tl2_ops::operator-(vec1); }

	friend vec operator*(value_type d, const vec& vec1) { return tl2_ops::operator*(d, vec1); }
	friend vec operator*(const vec& vec1, value_type d) { return tl2_ops::operator*(vec1, d); }
	friend vec operator/(const vec& vec1, value_type d) { return tl2_ops::operator/(vec1, d); }

	vec& operator*=(const vec& vec2) { return tl2_ops::operator*=(*this, vec2); }
	vec& operator+=(const vec& vec2) { return tl2_ops::operator+=(*this, vec2); }
	vec& operator-=(const vec& vec2) { return tl2_ops::operator-=(*this, vec2); }
	vec& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	vec& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }
};



/**
 * generic matrix container
 */
template<class T = double,
	template<class...> class t_cont = std::vector,
	template<class...> class t_alloc = __TLIBS2_DEFAULT_ALLOC__>
requires is_basic_vec<t_cont<T, t_alloc<T>>> && is_dyn_vec<t_cont<T, t_alloc<T>>>
class mat
{
public:
	using value_type = T;
	using container_type = t_cont<T, t_alloc<T>>;


	mat() = default;
	~mat() = default;

	mat(std::size_t ROWS, std::size_t COLS, const T* arr = nullptr)
		: m_data(ROWS*COLS), m_rowsize{ROWS}, m_colsize{COLS}
	{
		if(arr)
			from_array(arr);
	}

	mat(mat<T, t_cont, t_alloc>&& other) noexcept
		: m_data{ std::move(other.m_data) },
		  m_rowsize{ other.m_rowsize }, m_colsize{ other.m_colsize }
	{}

	mat(const mat<T, t_cont, t_alloc>& other)
		: //m_data{ other.m_data },
		  m_rowsize{ other.m_rowsize }, m_colsize{ other.m_colsize }
	{
		this->m_data.resize(other.m_data.size());
		for(std::size_t i = 0; i < other.m_data.size(); ++i)
			this->m_data[i] = other.m_data[i];
	}


	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other>
	mat(const mat<T_other, t_cont_other, t_alloc_other>& other)
	{
		this->operator=<T_other, t_cont_other, t_alloc_other>(other);
	}

	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other>
	mat<T, t_cont, t_alloc>& operator=(const mat<T_other, t_cont_other, t_alloc_other>& other)
	{
		*this = convert<mat<T, t_cont, t_alloc>, mat<T_other, t_cont_other, t_alloc_other>>(other);

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}

	mat<T, t_cont, t_alloc>& operator=(const mat<T, t_cont, t_alloc>& other)
	{
		//this->m_data = other.m_data;
		this->m_data.resize(other.m_data.size());
		for(std::size_t i = 0; i < other.m_data.size(); ++i)
			this->m_data[i] = other.m_data[i];

		this->m_rowsize = other.m_rowsize;
		this->m_colsize = other.m_colsize;

		return *this;
	}


	std::size_t size1() const { return m_rowsize; }
	std::size_t size2() const { return m_colsize; }


	// element access
	const T& operator()(std::size_t row, std::size_t col) const
	{
		return m_data[row*m_colsize + col];
	}

	T& operator()(std::size_t row, std::size_t col)
	{
		return m_data[row*m_colsize + col];
	}


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(std::size_t i = 0; i < m_rowsize; ++i)
			for(std::size_t j = 0; j < m_colsize; ++j)
				this->operator()(i, j) = arr[i*m_colsize + j];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(std::size_t i = 0; i < m_rowsize; ++i)
			for(std::size_t j = 0; j < m_colsize; ++j)
				arr[i*m_colsize + j] = this->operator()(i,j);
	}


	friend mat operator+(const mat& mat1, const mat& mat2) { return tl2_ops::operator+(mat1, mat2); }
	friend mat operator-(const mat& mat1, const mat& mat2) { return tl2_ops::operator-(mat1, mat2); }
	friend const mat& operator+(const mat& mat1) { return tl2_ops::operator+(mat1); }
	friend mat operator-(const mat& mat1) { return tl2_ops::operator-(mat1); }

	friend mat operator*(const mat& mat1, const mat& mat2) { return tl2_ops::operator*(mat1, mat2); }
	friend mat operator*(const mat& mat1, value_type d) { return tl2_ops::operator*(mat1, d); }
	friend mat operator*(value_type d, const mat& mat1) { return tl2_ops::operator*(d, mat1); }
	friend mat operator/(const mat& mat1, value_type d) { return tl2_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_cont<T, t_alloc<T>>> && is_dyn_vec<t_cont<T, t_alloc<T>>>
	friend t_vec operator*(const mat& mat1, const t_vec& vec2) { return tl2_ops::operator*(mat1, vec2); }

	mat& operator*=(const mat& mat2) { return tl2_ops::operator*=(*this, mat2); }
	mat& operator+=(const mat& mat2) { return tl2_ops::operator+=(*this, mat2); }
	mat& operator-=(const mat& mat2) { return tl2_ops::operator-=(*this, mat2); }
	mat& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	mat& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }


private:
	container_type m_data{ };
	std::size_t m_rowsize{ };
	std::size_t m_colsize{ };
};

// ----------------------------------------------------------------------------

}

#endif
