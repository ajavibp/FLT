/*  Copyright (C) 2004-2015
	ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
	http://uhu.es/antonio.barragan

	Collaborators:
	JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
	AGUSTIN JIMENEZ AVELLO, agustin.jimenez@upm.es
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
 * \example matlab_utilities/Kalmanconseq.cpp
 * MEX file that calculates an iteration of EKF only for consequents.
 *
 * MATLAB help:
 * \include Kalmanconseq.m
 *
 * Source code:
*/

#include "aux_matlab.hpp"
#include <flt/Kalman.hpp>

using namespace FLT;
using namespace TNT;

void mexFunction (int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

	// Reads and checks arguments
	if(nrhs<5 | nrhs>6)
		ERRORMSG(E_NumberArgIn)

	if (nlhs>3)
		ERRORMSG(E_NumberArgOut)

	for (size_t i=0;i<nrhs;i++)
		if (mxIsEmpty(prhs[i]) || mxIsNaN(*mxGetPr(prhs[i])))
			ERRORMSG(E_ArgNoValid)

	System Model;
	if (readModel(prhs[0],Model))
		ERRORMSG(E_Model)

	int m = Model.outputs();

	double *p_input = mxGetPr(prhs[1]);
	if (!p_input)
		ERRORMSG(E_NumberArgIn)
	int n = mxGetM(prhs[1]);
	if (n==1)
		n = mxGetN(prhs[1]);
	Array1D<double> input(n, p_input);

	double *p_output = mxGetPr(prhs[2]);
	if (!p_output)
		ERRORMSG(E_NumberArgIn)
	n = mxGetM(prhs[2]);
	if (n==1)
		n = mxGetN(prhs[2]);
	Array1D<double> output(n, p_output);

	double *p_covariance = mxGetPr(prhs[3]);
	if (!p_covariance)
		ERRORMSG(E_NumberArgIn)
	Array2D<double> covariance = col2row_major(p_covariance, m, m);

	double *p_P = mxGetPr(prhs[4]);
	if (!p_P)
		ERRORMSG(E_NumberArgIn)
	n = mxGetM(prhs[4]);
	Array2D<double> P = col2row_major(p_P, n, n);

	// Does the iteration
	Array1D<double> error;
	if (nrhs==5)
		error = KalmanConseq(Model, input, output, covariance, P);
	else
	{
		double *p_Phi = mxGetPr(prhs[5]);
		if (!p_Phi)
			ERRORMSG(E_NumberArgIn)
		Array2D<double> Phi = col2row_major(p_Phi, n, n);

		error = KalmanConseq(Model, input, output, covariance, P, Phi);
	}
	
	if (error.dim() == 0)
		ERRORMSG(E_NumberArgIn)

	// Generates outputs
	const char* names[]={"name","type","andMethod","orMethod","defuzzMethod","impMethod","aggMethod","input","output","rule"};
	mwSize dim[] = {1,1};

	plhs[0] = mxCreateStructArray(2,dim,10,names);

	if(System2FIS(Model, plhs[0]))
	{
		mxDestroyArray(plhs[0]);
		ERRORMSG(E_CreateFIS)
	}

	if (nlhs>1)
	{
		plhs[1] = mxCreateDoubleMatrix(n,n,mxREAL);
		double *p_newP = mxGetPr(plhs[1]);
		if (!p_newP)
		{
			mxDestroyArray(plhs[0]);
			mxDestroyArray(plhs[1]);
			ERRORMSG(E_NumberArgOut)
		}
		row2col_major(P[0], p_newP, n, n);

		if (nlhs>2)
		{
			plhs[2] = mxCreateDoubleMatrix(m,1,mxREAL);
			double *p_error = mxGetPr(plhs[2]);
			if (!p_error)
			{
				mxDestroyArray(plhs[0]);
				mxDestroyArray(plhs[1]);
				mxDestroyArray(plhs[2]);
				ERRORMSG(E_NumberArgOut)
			}
			for (size_t i=0;i<m;i++)
			{
				*p_error = error[i];
				p_error++;
			}
		}
	}
}
