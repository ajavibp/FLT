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
 * \example matlab_utilities/initializeController.cpp
 * MEX file that initializes a fuzzy Controller from a fuzzy Plant.
 *
 * MATLAB help:
 * \include initializeController.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char file[MAX_CHAR_LONG];
	
	if(nrhs!=1)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>1)
		ERRORMSG(E_NumberArgOut)
	if (mxIsEmpty(prhs[0]) || mxIsNaN(*mxGetPr(prhs[0])))
		ERRORMSG(E_ArgNoValid)
		
	System P,C;
	if (readModel(prhs[0],P))
		ERRORMSG(E_Model)

	C = initialController(P);
	if (!C.outputs())
	{
		mexErrMsgTxt("The input fuzzy model is not a valid Plant\n");
		return;
	}
	
	int numInMF = 0;
	mxArray *numOutMF = mxCreateDoubleMatrix(1,C.outputs(),mxREAL);
	for (int j=0;j<C.outputs();j++)
	{
		numInMF += C.rules(j);
		*(mxGetPr(numOutMF)+j) = C.rules(j);
	}

	mxArray *parameters[12];
	parameters[0] = mxCreateString("NumInputs");
	parameters[1] = mxCreateDoubleScalar(C.inputs());
	parameters[2] = mxCreateString("NumInputMFs");
	parameters[3] = mxCreateDoubleScalar(numInMF);
	parameters[4] = mxCreateString("NumOutputs");
	parameters[5] = mxCreateDoubleScalar(C.outputs());
	parameters[6] = mxCreateString("NumOutputMFs");
	parameters[7] = numOutMF;
	parameters[8] = mxCreateString("AddRules");
	parameters[9] = mxCreateString("none");
	parameters[10] = mxCreateString("Name");
	parameters[11] = mxCreateString(U_FISName);
	if (mexCallMATLAB(1, &plhs[0], 12, parameters, "sugfis"))
		ERRORMSG(E_FISOut)
	if(System2FIS(C,plhs[0]))
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_FISOut)
	}
}
