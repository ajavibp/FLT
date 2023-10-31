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
 * \example matlab_utilities/extractPoints.cpp
 * MEX file that extracts the significative points of a fuzzy model.
 *
 * MATLAB help:
 * \include extractPoints.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"

using namespace FLT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	char file[MAX_CHAR_LONG];
	size_t i,j;
	int error = 0;
	size_t numMeshPoints = 5;
	double precision = 1E-3;
	bool onlyStatePoints = true;
	double *d;
	
	if(nrhs<1 || nrhs>4)
		ERRORMSG(E_NumberArgIn)
	if (nlhs>1)
		ERRORMSG(E_NumberArgOut)
	if (mxIsEmpty(prhs[0]) || mxIsNaN(*mxGetPr(prhs[0])))
		ERRORMSG(E_ArgNoValid)

	System S;
	if (readModel(prhs[0],S))
		ERRORMSG(E_Model)

	if (nrhs>1)
	{
		d = mxGetPr(prhs[1]);
		if (!d)
			ERRORMSG(E_NumberArgIn)
			
		numMeshPoints = (size_t)(*d);
		if (numMeshPoints<2)
			ERRORMSG(E_N_MESH_P)
		if (nrhs>2)
		{
			d = mxGetPr(prhs[2]);
			if (!d)
				ERRORMSG(E_NumberArgIn)
			precision = *d;
			if (precision<=0.0)
				ERRORMSG(E_PRECISION)
			if (nrhs==4)
			{
				d=mxGetPr(prhs[3]);
				if (!d)
					ERRORMSG(E_NumberArgIn)
				onlyStatePoints = (*d<=0.0);
			}
		}
	}

	TNT::Array2D<double> aPoints = extractPoints(S,numMeshPoints,precision,onlyStatePoints);
	plhs[0] = mxCreateDoubleMatrix(aPoints.dim1(),aPoints.dim2(),mxREAL);
	double *Points = mxGetPr(plhs[0]);
	if (!Points)
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_NumberArgOut)
	}
	for (j=0;j<aPoints.dim2();j++)
	{
		for (i=0;i<aPoints.dim1();i++)
		{
			*Points = aPoints[i][j];
			Points++;
		}
	}
}
