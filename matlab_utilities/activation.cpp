/*	Copyright (C) 2004-2022
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
 * \example matlab_utilities/activation.cpp
 * MEX file that calculates the fulfillment degree and the derivative of a rule.
 *
 * MATLAB help:
 * \include activation.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i, output, rule, m;
	double *Point,*W,*dW, *SumW;
	
	if(nrhs<3 || nrhs>4)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>2)
		ERRORMSG(E_NumberArgOut)
	for (i=0;i<nrhs;i++)
	{
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	}
	
	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)

	if (1!=mxGetN(prhs[1]))
		ERRORMSG(E_Column)
	if (S.inputs()!=mxGetM(prhs[1]))
		ERRORMSG(E_PointCoherent)
	if ( (mxGetN(prhs[2])!=1) || (mxGetM(prhs[2])!=1))
		ERRORMSG(E_Out_No_Scalar)
	if (nrhs==4)
	{
		if ((mxGetN(prhs[3])!=1) || (mxGetM(prhs[3])!=1))
			ERRORMSG(E_R_No_Scalar)
		rule = ((size_t)*mxGetPr(prhs[3]))-1;
		if (rule>=S.rules(output))
			ERRORMSG(E_RulesExceded)
	}

	output = ((size_t)*mxGetPr(prhs[2]))-1;
	if (output>=S.outputs())
		ERRORMSG(E_OutputExceded)
	Point = mxGetPr(prhs[1]);

	if (nrhs==4)
	{// Only one rule
		plhs[0] = mxCreateDoubleScalar(1);
		W = mxGetPr(plhs[0]);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
		if (nlhs==1)
			(*W) = S.readRule(output,rule)->activation(Point);
		else
		{
			m = S.inputs();
			plhs[1] = mxCreateDoubleMatrix(1,m,mxREAL);
			dW = mxGetPr(plhs[1]);
			(*W) = S.readRule(output,rule)->activation(Point, dW);
		}
	}
	else
	{// All rules. dW can't be returned
		plhs[0] = mxCreateDoubleMatrix(1,S.rules(output),mxREAL);
		W = mxGetPr(plhs[0]);
		if (!plhs[0])
			ERRORMSG(E_NumberArgOut)
		if (nlhs==2)
		{
			plhs[1] = mxCreateDoubleScalar(1);
			if (!plhs[1])
				ERRORMSG(E_NumberArgOut)
			SumW = mxGetPr(plhs[1]);
			*SumW = 0.0;
		}
		for (rule=0;rule<S.rules(output);rule++)
		{
			(*W) = S.readRule(output,rule)->activation(Point);
			if (nlhs==2)
				(*SumW) += (*W);
			W++;
		}
	}
}
