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
#include <array>
#include <type_traits>

#include "decls.h"
#include "operators.h"


#ifndef __TLIBS2_DEFAULT_ALLOC__
	#define __TLIBS2_DEFAULT_ALLOC__ tl2::alloc_noinit
	//#define __TLIBS2_DEFAULT_ALLOC__ std::allocator
#endif

#ifndef __TLIBS2_DEFAULT_STATIC_SIZE__
	#define __TLIBS2_DEFAULT_STATIC_SIZE__ 4
#endif

#ifndef __TLIBS2_DEFAULT_VEC__
	// 0: std::vector, 1: tl2::vec
	#define __TLIBS2_DEFAULT_VEC__ 1
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
		size_type mat_size2 = m_mat.size2();

		size_type row = i / mat_size2;
		size_type col = i % mat_size2;

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
		return base_type::operator()(i, 0);
	}

	const T operator[](size_t i) const
	{
		return base_type::operator()(i, 0);
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
 * WARNING: without initialisation, functions like push_back also stop working in std::vector
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
// general iterator
// ----------------------------------------------------------------------------
template<class t_vec_, bool is_const = false>
struct gen_iterator
{
	using t_vec = std::conditional<is_const, const t_vec_, t_vec_>::type;
	using value_type = typename t_vec_::value_type;
	using size_type = typename t_vec_::size_type;
	using difference_type = typename t_vec_::difference_type;
	using pointer = std::conditional<is_const, const value_type*, value_type*>::type;
	using reference = std::conditional<is_const, const value_type&, value_type&>::type;

	explicit gen_iterator(t_vec *vec = nullptr, size_type idx = 0) : m_vec{vec}, m_idx{idx}
	{ }

	gen_iterator(const gen_iterator<t_vec_, is_const>& other)
	{
		operator=(other);
	}

	gen_iterator(gen_iterator<t_vec_, is_const>&& other) noexcept
		: m_vec{ std::exchange(other.m_vec, nullptr) }, m_idx{ other.m_idx }
	{ }

	gen_iterator<t_vec_, is_const>& operator=(const gen_iterator<t_vec_, is_const>& other)
	{
		m_vec = other.m_vec;
		m_idx = other.m_idx;

		return *this;
	}

	gen_iterator<t_vec_, is_const>& operator=(gen_iterator<t_vec_, is_const>&& other)
	{
		m_vec = std::exchange(other.m_vec, nullptr);
		m_idx = other.m_idx;

		return *this;
	}

	bool operator==(gen_iterator<t_vec_, is_const> iter) const
	{
		return m_vec == iter.m_vec && m_idx == iter.m_idx;
	}

	bool operator!=(gen_iterator<t_vec_, is_const> iter) const
	{
		return !operator==(iter);
	}

	reference operator*() const
	{
		return (*m_vec)[m_idx];
	}

	template<bool is_var = !is_const, typename = std::enable_if_t<is_var>>
	value_type& operator*()
	{
		return (*m_vec)[m_idx];
	}

	gen_iterator<t_vec_, is_const>& operator++()
	{
		++m_idx;
		return *this;
	}

	gen_iterator<t_vec_, is_const>& operator+=(size_type num)
	{
		m_idx += num;
		return *this;
	}

	gen_iterator<t_vec_, is_const> operator++(int)
	{
		gen_iterator<t_vec_, is_const> prev = *this;
		operator++();
		return prev;
	}

	private:
		t_vec* m_vec{};
		size_type m_idx{};
};
// ----------------------------------------------------------------------------



// ----------------------------------------------------------------------------
// vector containers
// ----------------------------------------------------------------------------

/**
 * generic vector container with static and dynamic components
 */
template<class T = double,
	template<class...> class t_cont = std::vector,
	template<class...> class t_alloc = __TLIBS2_DEFAULT_ALLOC__,
	std::size_t STATIC_SIZE = __TLIBS2_DEFAULT_STATIC_SIZE__>
requires is_basic_vec<t_cont<T, t_alloc<T>>> && is_dyn_vec<t_cont<T, t_alloc<T>>>
class vec
{
public:
	using value_type = T;
	using container_type = t_cont<T, t_alloc<T>>;
	using static_container_type = std::array<T, STATIC_SIZE>;

	using size_type = typename container_type::size_type;
	using difference_type = typename container_type::difference_type;
	using allocator_type = typename container_type::allocator_type;

	using iterator = gen_iterator<vec<T, t_cont, t_alloc, STATIC_SIZE>, false>;
	using const_iterator = gen_iterator<vec<T, t_cont, t_alloc, STATIC_SIZE>, true>;


public:
	// constructors
	vec() = default;

	~vec()
	{
		clear();
	}

	vec(const vec<T, t_cont, t_alloc, STATIC_SIZE>& other)
	{
		const size_type other_size = other.size();
		this->resize(other_size);
		for(size_type i = 0; i < other_size; ++i)
			(*this)[i] = other[i];
	}

	vec(vec<T, t_cont, t_alloc, STATIC_SIZE>&& other) noexcept
		: m_data{ std::move(other.m_data) },
			m_static_data{ std::move(other.m_static_data) },
			m_size{ std::exchange(other.m_size, 0) }
	{ }

	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other,
		std::size_t STATIC_SIZE_OTHER>
	vec(const vec<T_other, t_cont_other, t_alloc_other, STATIC_SIZE_OTHER>& other)
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

	vec(size_type SIZE, const T* arr = nullptr)
	{
		resize(SIZE);
		if(arr)
			from_array(arr);
	}


	// assignment operators
	vec<T, t_cont, t_alloc, STATIC_SIZE>& operator=(const vec<T, t_cont, t_alloc, STATIC_SIZE>& other)
	{
		const size_type other_size = other.size();
		resize(other_size);
		for(size_type i = 0; i < other_size; ++i)
			(*this)[i] = other[i];

		return *this;
	}

	template<class T_other,
		template<class...> class t_cont_other,
		template<class...> class t_alloc_other,
		std::size_t STATIC_SIZE_OTHER>
	vec<T, t_cont, t_alloc, STATIC_SIZE>& operator=(
		const vec<T_other, t_cont_other, t_alloc_other, STATIC_SIZE_OTHER>& other)
	{
		*this = convert<vec<T, t_cont, t_alloc, STATIC_SIZE>,
			vec<T_other, t_cont_other, t_alloc_other, STATIC_SIZE_OTHER>>(other);
		return *this;
	}


	// move operator
	vec<T, t_cont, t_alloc, STATIC_SIZE>& operator=(vec<T, t_cont, t_alloc, STATIC_SIZE>&& other)
	{
		m_data = std::move(other.m_data);
		m_static_data = std::move(other.m_static_data);
		m_size = std::exchange(other.m_size, 0);

		return *this;
	}


	void push_back(const T& t)
	{
		size_type num = size();
		resize(num + 1);
		(*this)[num] = t;
	}

	void emplace_back(T&& t)
	{
		size_type num = size();
		resize(num + 1);
		(*this)[num] = std::move(t);
	}


	// conversion from / to array
	void from_array(const T* arr)
	{
		// initialise from given array data
		for(size_type i = 0; i < size(); ++i)
			this->operator[](i) = arr[i];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(size_type i = 0; i < size(); ++i)
			arr[i] = this->operator[](i);
	}


	template<class t_vec = std::vector<T>>
	t_vec to_stdvec() const
	{
		t_vec vec;
		vec.reserve(size());

		for(size_type i = 0; i < size(); ++i)
			vec.push_back(this->operator[](i));

		return vec;
	}


	// resizing
	void clear()
	{
		m_data.clear();
		m_size = 0;
	}

	void resize(size_type size)
	{
		if(size == m_size)
			return;

		m_size = size;
		if(m_size > STATIC_SIZE)
		{
			const size_t dyn_size = m_size - STATIC_SIZE;
			m_data.resize(dyn_size);
		}
	}

	void reserve(size_type size)
	{
		if(size <= m_size)
			return;

		m_size = size;
		if(m_size > STATIC_SIZE)
		{
			const size_t dyn_size = m_size - STATIC_SIZE;
			m_data.reserve(dyn_size);
		}
	}

	size_type size() const { return m_size; }


	// iterators
	iterator begin() { return iterator{this, 0}; }
	iterator end() { return iterator(this, m_size); }
	const_iterator begin() const { return const_iterator(this, 0); }
	const_iterator end() const { return const_iterator(this, m_size); }


	// calculation operators
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


	// element access
	const value_type& operator()(size_type i) const { return operator[](i); }
	value_type& operator()(size_type i) { return operator[](i); }

	const value_type& operator[](size_type i) const
	{
		if(i < STATIC_SIZE)
			return m_static_data[i];
		else if(i < m_size)
			return m_data[i - STATIC_SIZE];

		static T dummy;
		return dummy;
	}

	value_type& operator[](size_type i)
	{
		if(i < STATIC_SIZE)
			return m_static_data[i];
		else if(i < m_size)
			return m_data[i - STATIC_SIZE];

		static T dummy;
		return dummy;
	}


private:
	container_type m_data{};
	static_container_type m_static_data{};
	size_type m_size{};
};



/**
 * statically sized vector
 */
template<class T, std::size_t SIZE>
class vec_static
{
public:
	using value_type = T;

	using size_type = std::size_t;
	using difference_type = std::ptrdiff_t;

	using iterator = gen_iterator<vec_static<T, SIZE>, false>;
	using const_iterator = gen_iterator<vec_static<T, SIZE>, true>;


public:
	// constructors
	constexpr vec_static() = default;
	constexpr ~vec_static() = default;

	constexpr vec_static(const vec_static<T, SIZE>& other)
	{
		for(size_type i = 0; i < SIZE; ++i)
			(*this)[i] = other[i];
	}

	constexpr vec_static(vec_static<T, SIZE>&& other) noexcept
		: m_data{ std::move(other.m_data) }
	{ }

	template<class t_vec_other> requires is_basic_vec<t_vec_other>
	constexpr vec_static(const t_vec_other& other)
	{
		this->operator=<t_vec_other>(other);
	}

	constexpr vec_static(const std::initializer_list<T>& lst)
	{
		size_type i = 0;
		for(const T& elem : lst)
		{
			(*this)[i] = elem;
			++i;

			if(i >= SIZE)
				break;
		}
	}

	// _SIZE is a dummy size arguments as the vector is statically sized
	constexpr vec_static([[__maybe_unused__]] size_type _SIZE, const T* arr = nullptr)
	{
		assert(_SIZE <= SIZE);

		if(arr)
			from_array(arr);
	}


	// assignment operators
	constexpr vec_static<T, SIZE>& operator=(const vec_static<T, SIZE>& other)
	{
		for(size_type i = 0; i < SIZE; ++i)
			(*this)[i] = other[i];

		return *this;
	}

	template<class t_vec_other> requires is_basic_vec<t_vec_other>
	constexpr vec_static<T, SIZE>& operator=(
		const t_vec_other& other)
	{
		*this = convert<vec_static<T, SIZE>, t_vec_other>(other);
		return *this;
	}


	// move operator
	constexpr vec_static<T, SIZE>& operator=(vec_static<T, SIZE>&& other)
	{
		m_data = std::move(other.m_data);
		return *this;
	}


	// conversion from / to array
	constexpr void from_array(const T* arr)
	{
		// initialise from given array data
		for(size_type i = 0; i < SIZE; ++i)
			this->operator[](i) = arr[i];
	}

	constexpr void to_array(T* arr) const
	{
		// write elements to array
		for(size_type i = 0; i < SIZE; ++i)
			arr[i] = this->operator[](i);
	}


	constexpr size_type size() const noexcept { return SIZE; }


	// iterators
	constexpr iterator begin() { return iterator{this, 0}; }
	constexpr iterator end() { return iterator(this, SIZE); }
	constexpr const_iterator begin() const { return const_iterator(this, 0); }
	constexpr const_iterator end() const { return const_iterator(this, SIZE); }


	// calculation operators
	friend constexpr vec_static operator+(const vec_static& vec1, const vec_static& vec2) { return tl2_ops::operator+(vec1, vec2); }
	friend constexpr vec_static operator-(const vec_static& vec1, const vec_static& vec2) { return tl2_ops::operator-(vec1, vec2); }
	friend const constexpr vec_static& operator+(const vec_static& vec1) { return tl2_ops::operator+(vec1); }
	friend constexpr vec_static operator-(const vec_static& vec1) { return tl2_ops::operator-(vec1); }

	friend constexpr vec_static operator*(value_type d, const vec_static& vec1) { return tl2_ops::operator*(d, vec1); }
	friend constexpr vec_static operator*(const vec_static& vec1, value_type d) { return tl2_ops::operator*(vec1, d); }
	friend constexpr vec_static operator/(const vec_static& vec1, value_type d) { return tl2_ops::operator/(vec1, d); }

	constexpr vec_static& operator*=(const vec_static& vec2) { return tl2_ops::operator*=(*this, vec2); }
	constexpr vec_static& operator+=(const vec_static& vec2) { return tl2_ops::operator+=(*this, vec2); }
	constexpr vec_static& operator-=(const vec_static& vec2) { return tl2_ops::operator-=(*this, vec2); }
	constexpr vec_static& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	constexpr vec_static& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }


	// element access
	constexpr const value_type& operator()(size_type i) const { return operator[](i); }
	constexpr value_type& operator()(size_type i) { return operator[](i); }

	constexpr const value_type& operator[](size_type i) const
	{
		return m_data[i];
	}

	constexpr value_type& operator[](size_type i)
	{
		return m_data[i];
	}


private:
	std::array<T, SIZE> m_data{};
};




// ----------------------------------------------------------------------------
// matrix containers
// ----------------------------------------------------------------------------

/**
 * generic dynamic matrix container
 */
#if __TLIBS2_DEFAULT_VEC__ == 0
template<class T, class t_cont = class std::vector<T/*, __TLIBS2_DEFAULT_ALLOC__<T>*/>>
#else
template<class T, class t_cont = vec<T, std::vector, __TLIBS2_DEFAULT_ALLOC__,
__TLIBS2_DEFAULT_STATIC_SIZE__*__TLIBS2_DEFAULT_STATIC_SIZE__>>
#endif
requires is_basic_vec<t_cont> && is_dyn_vec<t_cont>
class mat
{
public:
	using value_type = T;
	using size_type = decltype(t_cont{}.size());


	// constructors
	mat() = default;
	~mat() = default;

	mat(size_type ROWS, size_type COLS, const T* arr = nullptr)
	: m_data(ROWS*COLS), m_rowsize{ROWS}, m_colsize{COLS}
	{
		if(arr)
			from_array(arr);
	}

	mat(const mat<T, t_cont>& other)
	{
		this->operator=(other);
	}

	template<class T_other, class t_cont_other>
	mat(const mat<T_other, t_cont_other>& other)
	{
		this->operator=<T_other, t_cont_other>(other);
	}

	mat(mat<T, t_cont>&& other) noexcept
	: m_data{ std::move(other.m_data) },
	m_rowsize{ std::exchange(other.m_rowsize, 0) },
	m_colsize{ std::exchange(other.m_colsize, 0) }
	{ }


	// assignment operators
	template<class T_other, class t_cont_other>
	mat<T, t_cont>& operator=(const mat<T_other, t_cont_other>& other)
	{
		//*this = convert<mat<T, t_cont>, mat<T_other, t_cont_other>>(other);

		const t_cont_other& dat = other.GetData();
		const size_type other_size = other.size();
		this->m_data.resize(other_size);
		for(size_type i = 0; i < other_size; ++i)
			this->m_data[i] = T(dat[i]);

		this->m_rowsize = other.size1();
		this->m_colsize = other.size2();

		return *this;
	}

	mat<T, t_cont>& operator=(const mat<T, t_cont>& other)
	{
		return operator=<T, t_cont>(other);
	}


	// move operator
	mat<T, t_cont>& operator=(mat<T, t_cont>&& other)
	{
		m_data = std::move(other.m_data);
		m_rowsize = std::exchange(other.m_rowsize, 0);
		m_colsize = std::exchange(other.m_colsize, 0);

		return *this;
	}


	// sizes
	size_type size() const { return (size_type)m_data.size(); }
	size_type size1() const { return m_rowsize; }
	size_type size2() const { return m_colsize; }


	// element access
	const T& operator()(size_type row, size_type col) const
	{
		if(row < m_rowsize && col < m_colsize)
			return m_data[row*m_colsize + col];

		static T dummy{};
		return dummy;
	}

	T& operator()(size_type row, size_type col)
	{
		if(row < m_rowsize && col < m_colsize)
			return m_data[row*m_colsize + col];

		static T dummy{};
		return dummy;
	}


	void from_array(const T* arr)
	{
		// initialise from given array data
		for(size_type i = 0; i < m_rowsize; ++i)
			for(size_type j = 0; j < m_colsize; ++j)
				this->operator()(i, j) = arr[i*m_colsize + j];
	}

	void to_array(T* arr) const
	{
		// write elements to array
		for(size_type i = 0; i < m_rowsize; ++i)
			for(size_type j = 0; j < m_colsize; ++j)
				arr[i*m_colsize + j] = this->operator()(i, j);
	}


	friend mat operator+(const mat& mat1, const mat& mat2) { return tl2_ops::operator+(mat1, mat2); }
	friend mat operator-(const mat& mat1, const mat& mat2) { return tl2_ops::operator-(mat1, mat2); }
	friend const mat& operator+(const mat& mat1) { return tl2_ops::operator+(mat1); }
	friend mat operator-(const mat& mat1) { return tl2_ops::operator-(mat1); }

	friend mat operator*(const mat& mat1, const mat& mat2) { return tl2_ops::operator*(mat1, mat2); }
	friend mat operator*(const mat& mat1, value_type d) { return tl2_ops::operator*(mat1, d); }
	friend mat operator*(value_type d, const mat& mat1) { return tl2_ops::operator*(d, mat1); }
	friend mat operator/(const mat& mat1, value_type d) { return tl2_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_cont> && is_dyn_vec<t_cont>
	friend t_vec operator*(const mat& mat1, const t_vec& vec2) { return tl2_ops::operator*(mat1, vec2); }

	mat& operator*=(const mat& mat2) { return tl2_ops::operator*=(*this, mat2); }
	mat& operator+=(const mat& mat2) { return tl2_ops::operator+=(*this, mat2); }
	mat& operator-=(const mat& mat2) { return tl2_ops::operator-=(*this, mat2); }
	mat& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	mat& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }

	const t_cont& GetData() const { return m_data; }


private:
	t_cont m_data{ };
	size_type m_rowsize{ };
	size_type m_colsize{ };
};



/**
 * generic dynamic matrix container
 */
template<class T, std::size_t ROWS, std::size_t COLS>
class mat_static
{
public:
	using value_type = T;
	using size_type = std::size_t;


	// constructors
	constexpr mat_static() = default;
	constexpr ~mat_static() = default;

	// _ROWS and _COLS are dummy size arguments as the matrix is statically sized
	constexpr mat_static(
		[[__maybe_unused__]] size_type _ROWS,
		[[__maybe_unused__]] size_type _COLS,
		const T* arr = nullptr)
	{
		assert((ROWS >= _ROWS && COLS >= _COLS));

		if(arr)
			from_array(arr);
	}

	constexpr mat_static(const T* arr)
	{
		if(arr)
			from_array(arr);
	}

	constexpr mat_static(const mat_static<T, ROWS, COLS>& other)
	{
		this->operator=(other);
	}

	template<class t_mat_other>
	requires is_basic_mat<t_mat_other>
	constexpr mat_static(const t_mat_other& other)
	{
		this->operator=<t_mat_other>(other);
	}

	constexpr mat_static(mat_static<T, ROWS, COLS>&& other) noexcept
		: m_data{ std::move(other.m_data) }
	{ }


	// assignment operators
	template<class t_mat_other>
	requires is_basic_mat<t_mat_other>
	constexpr mat_static<T, ROWS, COLS>& operator=(const t_mat_other& other)
	{
		*this = convert<mat_static<T, ROWS, COLS>, t_mat_other>(other);
		return *this;
	}

	constexpr mat_static<T, ROWS, COLS>& operator=(const mat_static<T, ROWS, COLS>& other)
	{
		return operator=<mat_static<T, ROWS, COLS>>(other);
	}


	// move operator
	constexpr mat_static<T, ROWS, COLS>& operator=(mat_static<T, ROWS, COLS>&& other)
	{
		m_data = std::move(other.m_data);
		return *this;
	}


	// sizes
	constexpr size_type size() const noexcept { return ROWS*COLS; }
	constexpr size_type size1() const noexcept { return ROWS; }
	constexpr size_type size2() const noexcept { return COLS; }


	// element access
	constexpr const T& operator()(size_type row, size_type col) const
	{
		return m_data[row*COLS + col];
	}

	constexpr T& operator()(size_type row, size_type col)
	{
		return m_data[row*COLS + col];
	}


	constexpr void from_array(const T* arr)
	{
		// initialise from given array data
		for(size_type i = 0; i < ROWS; ++i)
			for(size_type j = 0; j < COLS; ++j)
				this->operator()(i, j) = arr[i*COLS + j];
	}

	constexpr void to_array(T* arr) const
	{
		// write elements to array
		for(size_type i = 0; i < ROWS; ++i)
			for(size_type j = 0; j < COLS; ++j)
				arr[i*COLS + j] = this->operator()(i, j);
	}


	friend constexpr mat_static operator+(const mat_static& mat1, const mat_static& mat2) { return tl2_ops::operator+(mat1, mat2); }
	friend constexpr mat_static operator-(const mat_static& mat1, const mat_static& mat2) { return tl2_ops::operator-(mat1, mat2); }
	friend const constexpr mat_static& operator+(const mat_static& mat1) { return tl2_ops::operator+(mat1); }
	friend constexpr mat_static operator-(const mat_static& mat1) { return tl2_ops::operator-(mat1); }

	friend constexpr mat_static operator*(const mat_static& mat1, const mat_static& mat2) { return tl2_ops::operator*(mat1, mat2); }
	friend constexpr mat_static operator*(const mat_static& mat1, value_type d) { return tl2_ops::operator*(mat1, d); }
	friend constexpr mat_static operator*(value_type d, const mat_static& mat1) { return tl2_ops::operator*(d, mat1); }
	friend constexpr mat_static operator/(const mat_static& mat1, value_type d) { return tl2_ops::operator/(mat1, d); }

	template<class t_vec> requires is_basic_vec<t_vec>
	friend constexpr t_vec operator*(const mat_static& mat1, const t_vec& vec2) { return tl2_ops::operator*(mat1, vec2); }

	constexpr mat_static& operator*=(const mat_static& mat2) { return tl2_ops::operator*=(*this, mat2); }
	constexpr mat_static& operator+=(const mat_static& mat2) { return tl2_ops::operator+=(*this, mat2); }
	constexpr mat_static& operator-=(const mat_static& mat2) { return tl2_ops::operator-=(*this, mat2); }
	constexpr mat_static& operator*=(value_type d) { return tl2_ops::operator*=(*this, d); }
	constexpr mat_static& operator/=(value_type d) { return tl2_ops::operator/=(*this, d); }

	//constexpr const std::array<T, ROWS*COLS>& GetData() const { return m_data; }


private:
	std::array<T, ROWS*COLS> m_data{ };
};

// ----------------------------------------------------------------------------

}

#endif
