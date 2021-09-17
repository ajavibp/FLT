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

/**
 * \example matlab_utilities/fuzderaffine.cpp
 * MEX file that gets the derivative of a fuzzy model with respect to its affine consequents.
 *
 * MATLAB help:
 * \include fuzderaffine.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"
#include <flt/derivatives.hpp>

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t output, rule, parameter;
	int error=0;
	double *Point, *dModel_daffine;
	
	if(nrhs!=2 && nrhs!=4)
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
	
	size_t n = S.inputs();
	size_t m = S.outputs();
		
	if (mxGetN(prhs[1])!=1)
		ERRORMSG(E_Column)
		
	if (n != mxGetM(prhs[1]))
		ERRORMSG(E_PointCoherent)
	
	Point = mxGetPr(prhs[1]);
	
	if (nrhs==4)
	{
		if ( (mxGetN(prhs[2])!=1) || (mxGetM(prhs[2])!=1))
			ERRORMSG(E_Out_No_Scalar)
		if ((mxGetN(prhs[3])!=1) || (mxGetM(prhs[3])!=1))
			ERRORMSG(E_R_No_Scalar)

		output = ((size_t)*mxGetPr(prhs[2]))-1;
		if (output >= m)
			ERRORMSG(E_OutputExceded)

		rule = ((size_t)*mxGetPr(prhs[3]))-1;
		if (rule >= S.rules(output))
			ERRORMSG(E_RulesExceded)

		plhs[0] = mxCreateDoubleScalar(1);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
			
		dModel_daffine = mxGetPr(plhs[0]);
		(*dModel_daffine) = derconseq(S, Point, output, rule, (size_t)0, error);
		if (error)
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(O_GeneralError)
		}
	}
	else // nrhs==2
	{			
		TNT::Array2D<double> dS_daffine = deraffine(S, Point);
		if (dS_daffine.dim1()==0 || dS_daffine.dim2()==0 )
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(O_GeneralError)
		}
		
		plhs[0] = mxCreateDoubleMatrix(dS_daffine.dim1(), dS_daffine.dim2(),mxREAL);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
		dModel_daffine = mxGetPr(plhs[0]);
		
		for (size_t j=0;j<dS_daffine.dim2();j++)
		{
			for (size_t i=0;i<dS_daffine.dim1();i++)
			{
				*dModel_daffine = dS_daffine[i][j];
				dModel_daffine++;
			}
		}
	}	
}
