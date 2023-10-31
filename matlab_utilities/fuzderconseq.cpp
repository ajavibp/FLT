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
 * \example matlab_utilities/fuzderconseq.cpp
 * MEX file that gets the derivative of a fuzzy model with respect to its consequents.
 *
 * MATLAB help:
 * \include fuzderconseq.m
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
	double *Point, *dModel_dconseq;
	
	if(nrhs!=2 && nrhs!=5)
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
	
	if (nrhs==5)
	{
		if ( (mxGetN(prhs[2])!=1) || (mxGetM(prhs[2])!=1))
			ERRORMSG(E_Out_No_Scalar)
		if ((mxGetN(prhs[3])!=1) || (mxGetM(prhs[3])!=1))
			ERRORMSG(E_R_No_Scalar)
		if ( (mxGetN(prhs[4])!=1) || (mxGetM(prhs[4])!=1))
			ERRORMSG(E_Param_No_Scalar)

		output = ((size_t)*mxGetPr(prhs[2]))-1;
		if (output >= m)
			ERRORMSG(E_OutputExceded)

		rule = ((size_t)*mxGetPr(prhs[3]))-1;
		if (rule >= S.rules(output))
			ERRORMSG(E_RulesExceded)

		parameter = ((size_t)*mxGetPr(prhs[4]));
		if (parameter > S.inputs())
			ERRORMSG(E_InputExceded)

		plhs[0] = mxCreateDoubleScalar(1);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
			
		dModel_dconseq = mxGetPr(plhs[0]);
		(*dModel_dconseq) = derconseq(S, Point, output, rule, parameter, error);
		if (error)
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(O_GeneralError)
		}
	}
	else // nrhs==2
	{			
		TNT::Array2D<double> dS_dconseq = derconseq(S, Point);
		if (dS_dconseq.dim1()==0 || dS_dconseq.dim2()==0 )
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(O_GeneralError)
		}
		
		plhs[0] = mxCreateDoubleMatrix(dS_dconseq.dim1(), dS_dconseq.dim2(),mxREAL);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
		dModel_dconseq = mxGetPr(plhs[0]);
		
		for (size_t j=0;j<dS_dconseq.dim2();j++)
		{
			for (size_t i=0;i<dS_dconseq.dim1();i++)
			{
				*dModel_dconseq = dS_dconseq[i][j];
				dModel_dconseq++;
			}
		}
	}	
}
