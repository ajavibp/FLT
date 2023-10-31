/*  Copyright (C) 2004-2022
	ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
	http://uhu.es/antonio.barragan

	Collaborator:
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
 * \example matlab_utilities/fis2txt.cpp
 * MEX file that writes the fuzzy model in an ASCII text file.
 *
 * MATLAB help:
 * \include fis2txt.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char Fichero[MAX_CHAR_LONG];
	char error[MAX_CHAR_LONG];
	if(nrhs!=2)
	{
		sprintf(error,"%s %s, %s",E_NumberArg,O_OfIn,U_SeeHelp);
		ERRORMSG(error)
	}
	if (nlhs>0)
	{
		sprintf(error,"%s %s, %s",E_NumberArg,O_OfOut,U_SeeHelp);
		ERRORMSG(error)
	}
	for (int i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)
	if(mxGetString(prhs[1],Fichero,MAX_CHAR_LONG))
		ERRORMSG(E_FileInput)

	System S;
	if(readModel(prhs[0],S))
		ERRORMSG(E_Model)

	if(System2TXT(S,Fichero))
		mexPrintf("Error using ==> fis2txt\n%s\n%s %s\n",E_FileOutput,U_OR,E_NoValidFP);
}
