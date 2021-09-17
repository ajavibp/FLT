/*  Copyright (C) 2004-2015
    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
    http://uhu.es/antonio.barragan

    Collaborators:
    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
    AGUSTIN JIMENEZ AVELLO, agustin.jimenez@upm.es
    MARIANO J. AZNAR, marianojose.aznar@alu.uhu.es
    BASIL M. AL-HADITHI, basil.alhadithi@upm.es

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

#ifndef _KALMAN_HPP_
#define _KALMAN_HPP_

/**
 * \file Kalman.hpp
 * \brief Implements the adaptation of a fuzzy model by the extended Kalman filter.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt \n
 * For GNU Scientific Library (GSL) documentation see http://gnu.org/software/gsl
 */
 
#include <flt/derivatives.hpp>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

namespace FLT // Use the FTL (Fuzzy Logic Tool) namespace
{
	/**
	 * \brief Computes antecedets adjust by the discrete extended Kalman filter.
	 * 
	 * @param Model Is the fuzzy system whose antecedents will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration.
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @param Phi Is the Jacobian matrix that relates the parameters to be set with the
	 * following value of these parameters. Usually this will be an identity matrix
	 * of appropriate dimensions.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.
	 */
	TNT::Array1D<double> KalmanAntec(System &Model,
											  TNT::Array1D<double> &input,
											  TNT::Array1D<double> &output,
											  TNT::Array2D<double> &covariance,
											  TNT::Array2D<double> &P,
											  TNT::Array2D<double> &Phi);

	/**
	 * \brief Computes antecedets adjust by the discrete extended Kalman filter
	 * where Phi is assumed to be the identity matrix.
	 *
	 * @param Model Is the fuzzy system whose antecedents will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration.
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.
	 */
	TNT::Array1D<double> KalmanAntec(System &Model,
											  TNT::Array1D<double> &input,
											  TNT::Array1D<double> &output,
											  TNT::Array2D<double> &covariance,
											  TNT::Array2D<double> &P);

    /**
	 * \brief Computes consequents adjust by the discrete extended Kalman filter.
	 * 
	 * @param Model Is the fuzzy system whose consequents will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration (current State vector, x(k)).
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @param Phi Is the Jacobian matrix that relates the parameters to be set with the
	 * following value of these parameters. Usually this will be an identity matrix
	 * of appropriate dimensions.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.
	 */
	TNT::Array1D<double> KalmanConseq(System &Model,
											   TNT::Array1D<double> &input,
											   TNT::Array1D<double> &output,
											   TNT::Array2D<double> &covariance,
											   TNT::Array2D<double> &P,
											   TNT::Array2D<double> &Phi);
											   
	/**
	 * \brief Computes consequents adjust by the discrete extended Kalman filter
	 * where Phi is assumed to be the identity matrix.
	 * 
	 * @param Model Is the fuzzy system whose consequents will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration (current State vector, x(k)).
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.
	 */
	TNT::Array1D<double> KalmanConseq(System &Model,
											  TNT::Array1D<double> &input,
											  TNT::Array1D<double> &output,
											  TNT::Array2D<double> &covariance,
											  TNT::Array2D<double> &P);

	/**
	 * \brief Computes the simultaneous adjustment of antecedets and consequents by
	 * the discrete extended Kalman filter.
	 * 
	 * @param Model Is the fuzzy system whose parameters will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration.
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @param Phi Is the Jacobian matrix that relates the parameters to be set with the
	 * following value of these parameters. Usually this will be an identity matrix
	 * of appropriate dimensions.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.
	 */
	TNT::Array1D<double> KalmanFuz(System &Model,
											  TNT::Array1D<double> &input,
											  TNT::Array1D<double> &output,
											  TNT::Array2D<double> &covariance,
											  TNT::Array2D<double> &P,
											  TNT::Array2D<double> &Phi);
											   
	/**
	 * \brief Computes the simultaneous adjustment of antecedets and consequents by
	 * the discrete extended Kalman filter where Phi is assumed to be the identity matrix.
	 * 
	 * @param Model Is the fuzzy system whose parameters will be adjusted.
	 * 
	 * @param input Is the input applied to the system in the current iteration.
	 * 
	 * @param output Is the real output of the system in the current iteration.
	 * 
	 * @param covariance Is the noise covariance matrix estimated from the hope operator.
	 * 
	 * @param P Is the covariance matrix of the filter.
	 * 
	 * @return Returns the error vector of the iteration, or an empty vector if there was some error.
	 *
	 * Model and P will be changed according to the Kalman discrete filte algorithm.	 */
	TNT::Array1D<double> KalmanFuz(System &Model,
											  TNT::Array1D<double> &input,
											  TNT::Array1D<double> &output,
											  TNT::Array2D<double> &covariance,
											  TNT::Array2D<double> &P);	
} // FLT

#endif
