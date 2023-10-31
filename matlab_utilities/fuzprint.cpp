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
 * \example matlab_utilities/fuzprint.cpp
 * MEX file that converts a fuzzy model in a text file with its linguistic representation.
 *
 * MATLAB help:
 * \include fuzprint.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

#define ACCURACY 10 // Default accuracy to print the system

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char File[MAX_CHAR_LONG];
	size_t i,j;
	int accuracy;

	if(nrhs<2)
		ERRORMSG(E_NumberArgIn)
	if(nlhs>0)
		ERRORMSG(E_NumberArgOut)
	for (int i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	if(mxGetString(prhs[0],File,MAX_CHAR_LONG))
		ERRORMSG(E_FileOutput)
	
	System S;
	if(readModel(prhs[1],S))
		ERRORMSG(E_Model)
	
	if(nrhs>2)
	{
		if ( (nrhs<(2+S.inputs()+S.outputs())) || (nrhs>(3+S.inputs()+S.outputs())) )
			ERRORMSG(E_InNoCoherent)
		
		char **inputs=new char *[S.inputs()];
		char **outputs=new char *[S.outputs()];
		if (nrhs==(2+S.inputs()+S.outputs()))
			accuracy = ACCURACY;
		else
		{
			if (!mxIsNumeric(prhs[2+S.inputs()+S.outputs()]))
			{
				for (i=0;i<S.inputs();delete []inputs[i++]);
				for (i=0;i<S.outputs();delete []outputs[i++]);
				delete [] inputs;
				delete [] outputs;
				ERRORMSG(E_AccuracyNoNum)
			}
			accuracy = 1+(int)mxGetScalar(prhs[2+S.inputs()+S.outputs()]);
		}
		
		// Read inputs and outputs' names
		for (i=0,j=2;i<S.inputs();i++,j++)
		{
			inputs[i] = new char [MAX_CHAR_LONG];
			if(mxGetString(prhs[j],inputs[i],MAX_CHAR_LONG))
			{
				for (i=0;i<S.inputs();delete []inputs[i++]);
				for (i=0;i<S.outputs();delete []outputs[i++]);
				delete [] inputs;
				delete [] outputs;
				ERRORMSG(E_InNames)
			}
		}
		for (i=0,j=2+S.inputs();i<S.outputs();i++,j++)
		{
			outputs[i] = new char [MAX_CHAR_LONG];
			if(mxGetString(prhs[j],outputs[i],MAX_CHAR_LONG))
			{
				for (i=0;i<S.inputs();delete []inputs[i++]);
				for (i=0;i<S.outputs();delete []outputs[i++]);
				delete [] inputs;
				delete [] outputs;
				ERRORMSG(E_OutNames)
			}
		}
		if(printSystem(File,S,inputs,outputs,accuracy))
		{
			for (i=0;i<S.inputs();delete []inputs[i++]);
			for (i=0;i<S.outputs();delete []outputs[i++]);
			delete [] inputs;
			delete [] outputs;
			ERRORMSG(E_FileOutput)
		}
		for (i=0;i<S.inputs();delete []inputs[i++]);
		for (i=0;i<S.outputs();delete []outputs[i++]);
		delete [] inputs;
		delete [] outputs;
	}
	else if(printSystem(File,S))
		ERRORMSG(E_FileOutput)
}
