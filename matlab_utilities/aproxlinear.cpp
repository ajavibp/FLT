/*  Copyright (C) 2004-2022
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

/**
 * \example matlab_utilities/aproxlinear.cpp
 * MEX file that calculates the linealization of a Fuzzy Model in a point usign finite differences.
 *
 * MATLAB help:
 * \include aproxlinear.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"
#include <flt/derivatives.hpp>

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i,q,j,n,m=0;
	size_t countA,countB;
	double *X,*U=NULL,*M_A,*M_B,*M_F,h;
	System S;
	char Sist[MAX_FILE_NAME];
	char error[MAX_FILE_NAME];
	
	if(nrhs<2 || nrhs>4)
		ERRORMSG(E_NumberArgIn)

	if (nlhs>3)
		ERRORMSG(E_NumberArgOut)
	for (i=0;i<nrhs;i++)
	{
		if ((i!=2 && mxIsEmpty(prhs[i])) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	}
	
	X = mxGetPr(prhs[1]);
	if (!X)
		ERRORMSG(E_NumberArgIn)
	
	if (mxGetN(prhs[1])!=1)
		ERRORMSG(E_Column)
	n = mxGetM(prhs[1]);
	
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)
		
	if (n!=S.outputs())
		ERRORMSG(E_PointCoherent)
		
	if (nrhs>=3)
	{
		U = mxGetPr(prhs[2]);
		if (!U)
			ERRORMSG(E_NumberArgIn)
		if (mxGetN(prhs[2])!=1)
			ERRORMSG(E_Column)
		m = mxGetM(prhs[2]);
		if (m!=(S.inputs()-n))
			ERRORMSG(E_PointCoherent)
	}
	else
		U = NULL;

	plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
	M_A = mxGetPr(plhs[0]);

	if (!plhs[0] || !M_A)
	{
		if (nrhs==2 && m!=0)
			delete []U;
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_NumberArgOut)
	}

	if (nlhs>1)
	{
		plhs[1] = mxCreateDoubleMatrix(n,m,mxREAL);
		M_B = mxGetPr(plhs[1]);
		if (!plhs[1] || !M_B)
		{
			if (nrhs==2 && m!=0)
				delete []U;
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			ERRORMSG(E_NumberArgOut)
		}
	}

	if (nlhs==3)
	{
		plhs[2] = mxCreateDoubleMatrix(n,1,mxREAL);
		M_F = mxGetPr(plhs[2]);
		if (!plhs[2] || !M_F)
		{
			if (nrhs==2 && m!=0)
				delete []U;
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			mxDestroyArray(plhs[2]);
			ERRORMSG(E_NumberArgOut)
		}
	}
	Array2D<double> A(n,n), B(n,m);
	Array1D<double> F(n);
	if (nrhs==4)
	{
		h=mxGetScalar(prhs[3]);
		if (h<=0)
			ERRORMSG(E_H_GE_0)

		A = jacobianAprox (S, X, U, B, F, h);
	}
	else
		A = jacobianAprox (S, X, U, B, F);
	if(!A.dim1())
	{
		if (nrhs==2 && m!=0)
			delete []U;
		mxDestroyArray(plhs[0]);
		if (nlhs>1)
			mxDestroyArray(plhs[1]);
		if (nlhs==3)
			mxDestroyArray(plhs[2]);
		ERRORMSG(U_Overflow)
	}
	for (q=0,countA=0;q<n;q++)
	{
		for (i=0;i<n;i++,countA++)
			*(M_A + countA) = A[i][q];
	}
	if (nlhs>1)
	{
		for (j=0,countB=0;j<m;j++)
			for (i=0;i<n;i++,countB++)
				*(M_B + countB) = B[i][j];
	}
	if (nlhs==3)
		for (q=0;q<n;q++)
			M_F[q] = F[q];
}
