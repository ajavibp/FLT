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
 * \example matlab_utilities/fuzcomb.cpp
 * MEX file that combines some fuzzy models in one.
 *
 * MATLAB help:
 * \include fuzcomb.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i,j,r,s;
	char Output[MAX_FILE_NAME];
	char error[MAX_FILE_NAME];
	size_t sysnumber,in,outputs=0,*rules;
	
	if (nlhs>1)
	{
		sprintf(error,"%s %s, %s",E_NumberArg,O_OfOut,U_SeeHelp);
		ERRORMSG(error)
	}
	if ((nlhs==0 && nrhs<3) || (nlhs==1 && nrhs<2))
	{
		sprintf(error,"%s, %s",E_NumberArg,U_SeeHelp);
		ERRORMSG(error)
	}
	int operationtype = 1 - nlhs; // 0 if output = FIS, 1 if output = TXT
	sysnumber = nrhs - operationtype;
	for (size_t i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)

	// Read models
	System *S = new System[sysnumber];
	Rule *R, *Rall;
	for (i=operationtype;i<nrhs;i++)
	{
		if(readModel(prhs[i],S[i-operationtype]))
		{
			mexPrintf("%s %d.\n",U_CheckModel,1+i-operationtype);
			delete [] S;
			ERRORMSG(E_Model)
		}
		
		// Check if the models are compatible ( = number of inputs)
		if (i==operationtype)
		{
			in = S[0].inputs();
		}
		else if (in!=S[i-operationtype].inputs())
		{
			delete [] S;
			ERRORMSG(E_SameIns)
		}
		
		outputs += S[i-operationtype].outputs();
	}

	// Combine the models
	rules = new size_t[outputs];
	j = 0;
	for (s=0;s<sysnumber;s++)
		for (i=0;i<S[s].outputs();i++,j++)
			rules[j] = S[s].rules(i);
			
	System *Sall = new System(in,outputs,rules);
	outputs = 0;
	for (s=0;s<sysnumber;s++)
	{
		for (j=0;j<S[s].outputs();j++,outputs++)
		{
			
			for (r=0;r<S[s].rules(j);r++)	
			{
				Rall = Sall->readRule(outputs,r);
				R = S[s].readRule(j,r);
				for (i=0;i<in;i++)
				{
					Rall->changeFunction(R->readFunction(i)->type(),i);
					(*Rall->readFunction(i)) = (*R->readFunction(i));
					Rall->changeTSK(R->readTSK(i),i);
				}
				Rall->changeTSK(R->readTSK(in),in);
			}
		}
	}
	
	double *limitsmin = new double[S[0].inputs()];
	double *limitsmax = new double[S[0].inputs()];
	for (i=0;i<S[0].inputs();i++)
	{
		limitsmin[i] = S[0].in_low(i);
		limitsmax[i] = S[0].in_high(i);
	}
	
	Sall->in_low(limitsmin);
	Sall->in_high(limitsmax);
	delete [] limitsmin;
	delete [] limitsmax;
	limitsmin = new double[outputs];
	limitsmax = new double[outputs];
	outputs=0;
	for (s=0;s<sysnumber;s++)
	{
		for (j=0;j<S[s].outputs();j++,outputs++)
		{
			limitsmin[outputs] = S[s].out_low(j);
			limitsmax[outputs] = S[s].out_high(j);
		}
	}
	
	Sall->out_low(limitsmin);
	Sall->out_high(limitsmax);
	delete [] limitsmin;
	delete [] limitsmax;

	// Ending
	if (nlhs==0)
	{
		if(mxGetString(prhs[0],Output,MAX_FILE_NAME))
		{
			delete [] S;
			delete Sall;
			delete [] rules;
			ERRORMSG(E_FileOutput)
		}
		if(System2TXT(*Sall,Output))
		{
			delete [] S;
			delete Sall;
			delete [] rules;
			mexPrintf("??? Error using ==> fuzcomb\n%s\n%s %s\n",E_FileOutput,U_OR,E_NoValidFP);
		}
	}
	else
	{
		const char* names[]={"name","type","andMethod","orMethod","defuzzMethod","impMethod","aggMethod","input","output","rule"};
		mwSize dim[]={1,1};
		plhs[0]=mxCreateStructArray(2,dim,10,names);
		if(System2FIS(*Sall, plhs[0]))
		{
			mxDestroyArray(plhs[0]);
			delete [] S;
			delete Sall;
			delete [] rules;
			ERRORMSG(E_CreateFIS)
		}
	}
	delete [] S;
	delete Sall;
	delete [] rules;
}
