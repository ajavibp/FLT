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
 * \example matlab_utilities/fuzsubsystem.cpp
 * MEX file that generates a subsystem (submodel) from a fuzzy model.
 *
 * MATLAB help:
 * \include fuzsubsystem.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i,j,n1,n2,m1,m2,n;
	size_t *outputs,*rules;
	size_t nOutputs, nRules;
	double *poutput,*prule;
	System S, Snew;
	char error[MAX_CHAR_LONG];
	if(nrhs!=3)
	{
		sprintf(error,"%s %s, %s",E_NumberArg,O_OfIn,U_SeeHelp);
		ERRORMSG(error)
	}
	if (nlhs>1)
	{
		sprintf(error,"%s %s, %s",E_NumberArg,O_OfOut,U_SeeHelp);
		ERRORMSG(error)
	}
	for (int i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
			
	n1 = mxGetN(prhs[1]);
	n2 = mxGetN(prhs[2]);
	m1 = mxGetM(prhs[1]);
	m2 = mxGetM(prhs[2]);
	if (n1!=n2 || m1!=m2)
	{
		sprintf(error,"%s The sizes of Outputs and Rules must be the same, %s",E_NumberArgIn,U_SeeHelp);
		ERRORMSG(error)
	}
	if (n1==1 && n2==1 && m1==m2)
		n=m1;
	else if (m1==1 && m2==1 && n1==n2)
		n=n1;
	else
	{
		sprintf(error,"%s, %s",E_NumberArgIn,U_SeeHelp);
		ERRORMSG(error)
	}

	if(readModel(prhs[0], S))
		ERRORMSG(E_Model)
	
	poutput = mxGetPr(prhs[1]);
	prule = mxGetPr(prhs[2]);
	outputs = new size_t[n];
	rules = new size_t[n];
	for (i=0;i<n;i++)
	{
		*(outputs+i) = ((size_t)*poutput)-1; // MATLAB is index-1 but C/C++ is index-0
		*(rules+i) = ((size_t)*prule)-1;
		poutput++;
		prule++;
	}
	
	Snew = subSystem(S,n,outputs,rules);
	if(Snew.outputs()==0)
	{
		delete []outputs;
		delete []rules;
		sprintf(error,"%s, %s\nAt least, 1 rule for each output is needed.\nCheck indices.",E_NumberArgIn,U_SeeHelp);
		ERRORMSG(error)
	}
	
	delete [] outputs;
	delete [] rules;
	const char* names[]={"name","type","andMethod","orMethod","defuzzMethod","impMethod","aggMethod","input","output","rule"};
	mwSize dim[] = {1,1};
	plhs[0] = mxCreateStructArray(2,dim,10,names);
	if(System2FIS(Snew,plhs[0]))
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_FISOut)
	}
}
