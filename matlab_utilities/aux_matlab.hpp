/*  Copyright (C) 2004-2015
    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
    http://uhu.es/antonio.barragan

    Collaborators:
    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
	MARIANO J. AZNAR, marianojose.aznar@alu.uhu.es

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

#ifndef _AUX_MATLAB_HPP_
#define _AUX_MATLAB_HPP_

/**
 * \file aux_matlab.hpp
 * \brief Defines some util functions to use with MATLAB© API.
 * 
 * \sa MATLAB© is a registered trademark of The MathWorks, Inc. http://www.mathworks.com
 */
 
#include "mex.h" // MATLAB© API
#include <flt/utilities.hpp>

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
	/**
	 * \brief Converts a column-major matrix in a row-major one.
	 * 
	 * MATLAB© uses column-major matrices, but C uses row-major.
	 */
	inline void col2row_major(const double * const colM, double *rowM, size_t num_rows, size_t num_cols)
	{
		double *p = rowM;
		for (size_t i=0; i<num_rows; i++)
			for (size_t j=0; j<num_cols; j++, p++)
				*p = colM[i+j*num_rows];
	};
	
	/**
	 * \brief Converts a column-major matrix in a row-major one.
	 * 
	 * MATLAB© uses column-major matrices, but C uses row-major.
	 */
	inline TNT::Array2D<double> col2row_major(const double * const colM, size_t num_rows, size_t num_cols)
	{
		TNT::Array2D<double> M (num_rows, num_cols);
		double *p = M[0];
		for (size_t i=0; i<num_rows; i++)
			for (size_t j=0; j<num_cols; j++, p++)
				*p = colM[i+j*num_rows];
		return M;
	};

	/**
	 * \brief Converts a row-major matrix in a column-major one.
	 *
	 * MATLAB© uses column-major matrices, but C uses row-major.
	 */
	inline void row2col_major(const double * const rowM, double *colM, size_t num_rows, size_t num_cols)
	{
		double *p = colM;
		for (size_t i=0; i<num_cols; i++)
			for (size_t j=0; j<num_rows; j++, p++)
				*p = rowM[i+j*num_cols];
/*    for (i=0;i<n;i++)
    {
        for (q=0;q<n;q++)
        {
            *J = Jac[q][i];
            J++;
        }
    }*/
	};

	/**
	 * \brief Reads a Fuzzy Model (TXT, FIS file or FIS variable)
	 * 
	 * Returns 0 if no errors.
	 */
	int readModel(const mxArray *model, System &S);

	/**
	 * \brief Reads the Fuzzy Model from a FIS variable
	 * 
	 * Returns 0 if no errors.
	 */
	int FIS2System(const mxArray *FIS, System &S);

	/**
	 * \brief Writes the Fuzzy Model in a FIS variable
	 * 
	 * Returns 0 if no errors.
	 */
	int System2FIS(System &S, mxArray *FIS);
	
} // FLT
#endif
