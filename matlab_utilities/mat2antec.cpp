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
 * \example matlab_utilities/mat2antec.cpp
 * MEX file that changes the antecedent of a fuzzy model.
 *
 * MATLAB help:
 * \include mat2antec.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i, rule, out, length;
	double *data;
	
	if(nrhs<2)
		ERRORMSG(E_NumberArgIn)
	if (nrhs>2 && nrhs!=4)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>2)
		ERRORMSG(E_NumberArgOut)
	if (mxIsEmpty(prhs[0]) || mxIsNaN(*mxGetPr(prhs[0])))
		ERRORMSG(E_ArgNoValid)

	if (nrhs==4)
	{
		if (mxGetM(prhs[2])>1 || mxGetM(prhs[3])>1)
			ERRORMSG(E_RuleOutputFileVector)
		if (mxGetN(prhs[2])!=mxGetN(prhs[3]))
			ERRORMSG(E_NumberArg)
	}

	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)
	
	if (nlhs==2)
	{
		if (nrhs==2)
		{
			length = S.NumberOfAntecedents();
			if (!length)
				ERRORMSG(E_BadModel)
		}
		else
		{
			length = 0;
			for (i=0;i<mxGetN(prhs[2]);i++)
			{
				out = ((size_t)*(mxGetPr(prhs[2])+i))-1; // MATLAB is index-1 but C/C++ is index-0
				rule = ((size_t)*(mxGetPr(prhs[3])+i))-1;
				length += S.readRule(out,rule)->NumberOfAntecedents();
			}
			
			if (!length)
				mexPrintf("??? Error using ==> mat2antec\n%s, or\n%s\n",E_BadModel,E_NumberArgIn);
		}
	}

	// Save data in S
	data = mxGetPr(prhs[1]);
	int error = 0;
	
	if (nrhs==2)
		error = S.setAntecedents(data);
	else
	{
		for (i=0;i<mxGetN(prhs[2]);i++)
		{
			out = ((size_t)*(mxGetPr(prhs[2])+i))-1; // MATLAB is index-1 but C/C++ is index-0
			rule = ((size_t)*(mxGetPr(prhs[3])+i))-1;
			size_t numAntec = S.readRule(out,rule)->NumberOfAntecedents();
			error += S.readRule(out,rule)->setAntecedents(data);
			data += numAntec;
		}
	}
	
	if (error>0)
	{
		ERRORMSG(E_Store)
		plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL); // Empty matrix
		return;
	}
		
	if (nlhs==2)
		plhs[1] = mxCreateDoubleScalar(length);
		
	// Create the FIS variable
	const char* names[] = {"name","type","andMethod","orMethod","defuzzMethod","impMethod","aggMethod","input","output","rule"};
	mwSize dim[] = {1,1};
	plhs[0] = mxCreateStructArray(2,dim,10,names);

	if(System2FIS(S, plhs[0]))
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_CreateFIS)
	}
}
