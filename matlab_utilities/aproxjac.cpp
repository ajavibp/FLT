/*	Copyright (C) 2004-2015
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

/**
 * \example matlab_utilities/aproxjac.cpp
 * MEX file that calculates the aproximation to the Jacobian matrix of the closed loop fuzzy system usign finite differences.
 *
 * MATLAB help:
 * \include aproxjac.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"
#include <flt/derivatives.hpp>

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i,q,n,model;
	double *X,*J,*d,h;
	System S[2];
	char Sist[MAX_FILE_NAME];
	char error[MAX_FILE_NAME];
	
	if(nrhs!=3 && nrhs!=4)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>1)
		ERRORMSG(E_NumberArgOut)
	for (i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	
	X = mxGetPr(prhs[0]);
	if (!X)
		ERRORMSG(E_NumberArgIn)
		
	n = mxGetM(prhs[0]);
	
	if (mxGetN(prhs[0])!=1)
		ERRORMSG(E_Column)
	
	for (model=1;model<3;model++)
	{
		if (readModel(prhs[model],S[model-1]))
		{
			sprintf(error,"%s %lu.\n",E_Model, model);
			ERRORMSG(error)
		}
	}
	
	if ((nrhs==3) && n!=S[0].outputs())
		ERRORMSG(E_PointCoherent)
		
	if ((nrhs==2) && n!=S[0].inputs())
		ERRORMSG(E_PointCoherent)
	
	Array2D<double> Jac(n,n);
	if (nrhs==4)
	{
		h = *mxGetPr(prhs[3]);
		if (h<=0)
			ERRORMSG(E_H_GE_0)
		Jac = jacobianAprox(S[0], S[1], X, h);
	}
	else
		Jac = jacobianAprox(S[0], S[1], X);
		
	if (!Jac.dim1())
		ERRORMSG(E_NoCoherent)
	
	plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
	J = mxGetPr(plhs[0]);
	if (!J)
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_NumberArgOut)
	}
	
	for (i=0;i<n;i++)
	{
		for (q=0;q<n;q++)
		{
			*J = Jac[q][i];
			J++;
		}
	}
}
