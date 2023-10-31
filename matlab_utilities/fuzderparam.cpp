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
 * \example matlab_utilities/fuzderparam.cpp
 * MEX file that gets the jacobian matrix of a fuzzy model with respect to its parameters.
 *
 * MATLAB help:
 * \include fuzderparam.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"
#include <flt/derivatives.hpp>

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t input, output, rule, param, m;
	int error=0;
	double *Point, *dModel_dparam;
	
	if(nrhs!=2)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>1)
		ERRORMSG(E_NumberArgOut)
	for (int i=0;i<nrhs;i++)
	{
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	}
	
	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)
	
	if (1 != mxGetN(prhs[1]))
		ERRORMSG(E_Column)
	if (S.inputs() != mxGetM(prhs[1]))
		ERRORMSG(E_PointCoherent)
	Point = mxGetPr(prhs[1]);

	m = S.outputs();
	plhs[0]=mxCreateDoubleMatrix(m,S.NumberOfAntecedents() + S.NumberOfConsequents(),mxREAL);
	if (!plhs[0])
		ERRORMSG(E_NumberArgOut)
		
	dModel_dparam = mxGetPr(plhs[0]);
	TNT::Array2D<double> dS_dparam = derfuzzy(S, Point);
	if (dS_dparam.dim1()==0 || dS_dparam.dim2()==0 )
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(O_GeneralError)
	}
	
	for (size_t j=0;j<dS_dparam.dim2();j++)
	{
		for (size_t i=0;i<dS_dparam.dim1();i++)
		{
			*dModel_dparam = dS_dparam[i][j];
			dModel_dparam++;
		}
	}
}
