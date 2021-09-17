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

/**
 * \example matlab_utilities/conseq2mat.cpp
 * MEX file that extracts data from consequent of a fuzzy model.
 *
 * MATLAB help:
 * \include conseq2mat.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i, NoC, length, out, rule;
	double *data;
	
	if(nrhs<1)
		ERRORMSG(E_NumberArgIn)
	if (nrhs>1 && nrhs!=3)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>2)
		ERRORMSG(E_NumberArgOut)
	if (mxIsEmpty(prhs[0]) || mxIsNaN(*mxGetPr(prhs[0])))
		ERRORMSG(E_ArgNoValid)
	if (nrhs==3)
	{
		if (mxGetM(prhs[1])>1 || mxGetM(prhs[2])>1)
			ERRORMSG(E_RuleOutputFileVector)
		if (mxGetN(prhs[1])!=mxGetN(prhs[2]))
			ERRORMSG(E_NumberArg)
	}
	
	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)
	size_t m = S.inputs();
	if (nrhs==1)
	{
		length = S.NumberOfConsequents();
		if (length==0)
			ERRORMSG(E_BadModel)
	}
	else
		length = mxGetN(prhs[1])*(m+1);
	
	plhs[0] = mxCreateDoubleMatrix(length,1,mxREAL);
	Array1D<double> Data(length);
	
	if (nrhs==1)
		Data = S.getConsequents();
	else
	{
		length = 0;
		for (i=0;i<mxGetN(prhs[1]);i++)
		{
			out = ((size_t)*(mxGetPr(prhs[1])+i))-1;
			rule = ((size_t)*(mxGetPr(prhs[2])+i))-1;
			Array1D<double> tempData = S.readRule(out,rule)->getConsequents();
			for (size_t d=0;d<tempData.dim();d++)
				Data[length+d] = tempData[d];
			length += tempData.dim();
		}
	}

	double *MATLABData = mxGetPr(plhs[0]);
	for (i=0; i<length; i++)
		MATLABData[i] = Data[i];

	if (nlhs==2)
		plhs[1] = mxCreateDoubleScalar(length);
		
	if (length==0)
	{
		if (nlhs==1)
			ERRORMSG(E_Acquire)
		else
			plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL); // Empty matrix
	}
}
