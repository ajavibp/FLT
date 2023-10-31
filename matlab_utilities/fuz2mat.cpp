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
 * \example matlab_utilities/fuz2mat.cpp
 * MEX file that extracts all data from a fuzzy model.
 *
 * MATLAB help:
 * \include fuz2mat.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	size_t i, length, out, rule;
	double *data;
	
	if(nrhs<1)
		ERRORMSG(E_NumberArgIn)
	if (nrhs>1)
	{
		if (nrhs!=3)
			ERRORMSG(E_NumberArgIn)
		if (mxGetM(prhs[1])>1 || mxGetM(prhs[2])>1)
			ERRORMSG(E_RuleOutputFileVector)
		if (mxGetN(prhs[1])!=mxGetN(prhs[2]))
			ERRORMSG(E_NumberArg)
	}
	if (nlhs>2)
		ERRORMSG(E_NumberArgOut)
	if (mxIsEmpty(prhs[0]) || mxIsNaN(*mxGetPr(prhs[0])))
		ERRORMSG(E_ArgNoValid)
		
	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)
	
	size_t NoC = S.inputs() + 1;
	if (nrhs==1)
	{
		length = S.NumberOfAntecedents() + S.NumberOfConsequents();
	}
	else
	{
		length = 0;
		for (i=0;i<mxGetN(prhs[2]);i++)
		{
			out = ((size_t)*(mxGetPr(prhs[1])))-1; // MATLAB is index-1 but C/C++ is index-0
			rule = ((size_t)*(mxGetPr(prhs[2])))-1;
			size_t numAntec = S.readRule(out, rule)->NumberOfAntecedents();
			length += numAntec + NoC;
		}
	}

	Array1D<double> Data(length);
	
	if (nrhs==1)
	{
		Array1D<double> AntecData = S.getAntecedents();
		Array1D<double> ConseqData = S.getConsequents();
		for (size_t d=0;d<AntecData.dim();d++)
			Data[d] = AntecData[d];
		for (size_t d=0;d<ConseqData.dim();d++)
			Data[AntecData.dim()+d] = ConseqData[d];		
	}
	else
	{
		size_t dataIdx = 0;
		for (i=0;i<mxGetN(prhs[2]);i++)
		{
			out = ((size_t)*(mxGetPr(prhs[1])))-1; // MATLAB is index-1 but C/C++ is index-0
			rule = ((size_t)*(mxGetPr(prhs[2])))-1;
			Rule *R = S.readRule(out,rule);

			Array1D<double> AntecData = R->getAntecedents();
			for (size_t d=0;d<AntecData.dim();d++, dataIdx++)
				Data[dataIdx] = AntecData[d];

			Array1D<double> ConseqData = R->getConsequents();
			for (size_t d=0;d<ConseqData.dim();d++, dataIdx++)
				Data[dataIdx] = ConseqData[d];
		}
	}

	plhs[0] = mxCreateDoubleMatrix(length,1,mxREAL);
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
