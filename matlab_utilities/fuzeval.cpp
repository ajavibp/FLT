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
 * \example matlab_utilities/fuzeval.cpp
 * MEX file that performs fuzzy inference calculations of open/closed loop fuzzy system.
 *
 * MATLAB help:
 * \include fuzeval.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i;
	int ode = 0;
	
	if(nrhs<2 || nrhs>4)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>1)
		ERRORMSG(E_NumberArgOut)
	for (i=0;i<nrhs;i++)
	{
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	}
	
	double *X;
	Array1D<double> dX;
	
	System S,C;
	char error[MAX_FILE_NAME];

	if (nrhs!=2 && mxGetN(prhs[0])==1 && mxGetM(prhs[0])==1 && mxIsStruct(prhs[1])==false && mxIsChar(prhs[1])==false)
		ode = 1; // If this function is called from 'ode', the 1st argument is the time

	// Read arguments
	X = mxGetPr(prhs[ode]);
	if (!X)
		ERRORMSG(E_NumberArgIn)
	if (mxGetN(prhs[ode])!=1)
		ERRORMSG(E_Column)
	if (readModel(prhs[1+ode],S))
	{
		mexPrintf("%s %s.\n",U_CheckModel,O_OfPlant);
		ERRORMSG(E_Model)
	}
	
	plhs[0] = mxCreateDoubleMatrix(S.outputs(),1,mxREAL);
	if (!plhs[0])
		ERRORMSG(E_NumberArgOut)
	
	// Run
	if (nrhs==(2+ode))
	{// Evaluate the plant
		if (mxGetM(prhs[ode])!=S.inputs())
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(E_PointCoherent)
		}
		dX = S.evaluate(X);
		if(!dX.dim())
		{
			mxDestroyArray(plhs[0]);
			mexPrintf("Fuzeval. %s\b or %s\n",U_Overflow,E_NoCoherent);
			mexErrMsgTxt("");
			return;
		}
	}
	else
	{// Evaluate the closed-loop system (plant + controller)
		if (mxGetM(prhs[ode])!=S.outputs())
		{
			mxDestroyArray(plhs[0]);
			ERRORMSG(E_PointCoherent)
		}
		if (readModel(prhs[2+ode],C))
		{
			mxDestroyArray(plhs[0]);
			mexPrintf("%s %s.\n",U_CheckModel,O_OfController);
			ERRORMSG(E_Model)
		}
		dX = evaluate(X, S, C);
		if (!dX.dim())
		{
			if (!error)
			{
				mxDestroyArray(plhs[0]);
				ERRORMSG(E_NoCoherent_BadX)
			}
			else
			{
				mxDestroyArray(plhs[0]);
				mexPrintf("Fuzeval. %s or",U_Overflow);
				ERRORMSG(E_NoCoherent)
			}
		}
	}
	
	double *MATLABOut = mxGetPr(plhs[0]);
	for (i=0;i<dX.dim();i++)
		MATLABOut[i] = dX[i];
}
