/*  Copyright (C) 2004-2015
    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
    http://uhu.es/antonio.barragan

    Collaborators:
    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es

    DPTO. DE ING. ELECTRONICA, DE SISTEMAS INFORMATICOS Y AUTOMATICA
    ETSI, UNIVERSITY OF HUELVA (SPAIN)

    For more information, please contact with authors.
    
    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _TNT_EXTENSION_HPP_
#define _TNT_EXTENSION_HPP_

/**
 * \file tnt_extension.hpp
 * \brief Implements useful functions for the management of types defined in the TNT library.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */

#include <tnt/tnt.h> // Template Numerical Toolkit (http://math.nist.gov/tnt)

namespace TNT{
	
	/**
	 * \brief Creates a TNT::Vector from a TNT::Array1D.
	 */
	template <class T>
	inline Vector<T> Array1D2Vector(Array1D<T> A)
	{
		Vector<T> V(A.dim(), &A[0]);
		return V;
	}

	/**
	 * \brief Creates a TNT::Array1D from a TNT::Vector.
	 */
	template <class T>
	inline Array1D<T> Vector2Array1D(Vector<T> V)
	{
		Array1D<T> A(V.dim(), &V[0]);
		return A;
	}

	/**
	 * \brief Creates a TNT::Matrix from a TNT::Array2D.
	 */
	template <class T>
	inline Matrix<T> Array2D2Matrix(Array2D<T> A)
	{
		Matrix<T> M(A.dim1(), A.dim2(), A[0]);
		return M;
	}

	/**
	 * \brief Creates a TNT::Array2D from a TNT::Matrix.
	 */
	template <class T>
	inline Array2D<T> Matrix2Array2D(Matrix<T> M)
	{
		Array2D<T> A(M.dim1(), M.dim2(), M[0]);
		return A;
	}

	/**
	 * \brief Creates a TNT::Matrix from a TNT::Array1D.
	 */
	template <class T>
	inline Matrix<T> Array1D2Matrix(Array1D<T> A)
	{
		Matrix<T> M(1, A.dim(), &A[0]);
		return M;
	}

	/**
	 * \brief Creates a TNT::Array1D from a TNT::Matrix.
	 * 
	 * If one dimension of the matrix is not ​​1, returns an empty array.
	 */
	template <class T>
	inline Array1D<T> Matrix2Array1D(Matrix<T> M)
	{
		if (M.dim1() == 1)
		{
			Array1D<T> A(M.dim2(), M[0]);
			return A;
		}
		else if (M.dim2() == 1)
		{
			Array1D<T> A(M.dim1(), M[0]);
			return A;
		}
		else
		{
			Array1D<T> A(0,0);
			return A;
		}	
	}

	/**
	 * \brief Creates an Identity Array2D of the given \c order.
	 */
	inline Array2D<double> makeIdentity(int order)
	{
		Array2D<double> I(order,order,0.0);
		for (int i=0;i<order;i++)
			I[i][i]=1.0;
		return I;
	}

	/**
	 * \brief Extracts the \c r row from an Array2D.
	 */
	template <class T>
	inline Array1D<T> copyrow(const Array2D<T> &Array,
	                          int r)
	{
		if (r >= Array.dim1())
		{
			Array1D<T> A;
			return A;
		}
		int m = Array.dim2();
		Array1D<T> A(m);
		for (int c=0; c<m; c++)
			A[c] = Array[r][c];
		return A;
	}

	/**
	 * \brief Extracts the \c c column from an Array2D.
	 */
	template <class T>
	inline Array1D<T> copycolumn(const Array2D<T> &Array,
	                             int c)
	{
		if (c >= Array.dim2())
		{
			Array1D<T> A;
			return A;
		}
		int n = Array.dim1();
		Array1D<T> A(n);
		for (int r=0; r<n; r++)
			A[r] = Array[r][c];
		return A;
	}

	/**
	 * \brief Extracts the \c m Array2D from an Array3D.
	 */
	template <class T>
	inline Array2D<T> copyArray2D(const Array3D<T> &Array,
	                              int m)
	{
		if (m >= Array.dim1())
		{
			Array2D<T> A;
			return A;
		}
		int n = Array.dim2();
		int o = Array.dim3();
		Array2D<T> A(n, o);
		for (int r=0; r<n; r++)
			for (int c=0; c<o; c++)
				A[r][c] = Array[m][r][c];
		return A;
	}

	/**
	 * \brief Injects a row into an Array2D in the position \c r.
	 */
	template <class T>
	inline int injectrow(Array2D<T> &A2D,
	                     const Array1D<T> &A1D,
	                     int r)
	{
		int m = A2D.dim2();
		if (m == A1D.dim() && r < A2D.dim1())
		{
			for (int c=0; c<m; c++)
				A2D[r][c] = A1D[c];
			return 0;
		}
		else
			return 1;
	}

	/**
	 * \brief Injects a column into an Array2D in the position \c c.
	 */
	template <class T>
	inline int injectcolumn(Array2D<T> &A2D,
	                        const Array1D<T> &A1D,
							int c)
	{
		int n = A2D.dim1();
		if (n == A1D.dim() && c < A2D.dim2())
		{
			for (int r=0; r<n; r++)
				A2D[r][c] = A1D[r];
			return 0;
		}
		return 1;
	}

	/**
	 * \brief Injects an Array2D into an Array3D in the position \c m.
	 */
	template <class T>
	inline int injectArray2D(Array3D<T> &A3D,
	                         const Array2D<T> &A2D,
							 int m)
	{
		int n = A3D.dim2();
		int o = A3D.dim3();
		if (n == A2D.dim1() && o == A2D.dim2() && m < A3D.dim1())
		{
			for (int r=0; r<n; r++)
				for (int c=0; c<o; c++)
					A3D[m][r][c] = A2D[r][c];
			return 0;
		}
		else
			return 1;
	}

	/**
	 * \brief Matrix-Matrix tranpose multiplication, i.e. compute A*tranpose(B).
	 * 
	 * NOTE: this is more efficient than computing the tranpose(B) explicitly,
	 * and then multiplying, as the tranpose of B is never really constructed.
	 *
	 * This function is based on transpose_mult of TNT.
	 * @param A  matrix: size N x M.
	 * @param B	matrix: size K x M.
	 * @return Returns a new matrix of size N x M.
	*/
	template <class T>
	inline Matrix<T> mult_transpose(const Matrix<T>  &A,
	                                const Matrix<T> &B)
	{

	#ifdef TNT_BOUNDS_CHECK
		assert(A.num_cols() == B.num_cols());
	#endif

		Subscript M = A.num_cols();
		Subscript N = A.num_rows();
		Subscript K = B.num_rows();

		Matrix<T> tmp(N,M);
		T sum;

		for (Subscript i=0; i<N; i++)
		for (Subscript k=0; k<K; k++)
		{
			sum = 0;
			for (Subscript j=0; j<M; j++)
				sum = sum +  A[i][j] * B[k][j];

			tmp[i][k] = sum; 
		}

		return tmp;
	}

} // TNT

#endif
