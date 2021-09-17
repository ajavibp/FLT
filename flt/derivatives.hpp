/*  Copyright (C) 2004-2015
    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
    http://uhu.es/antonio.barragan

    Collaborators:
    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
    AGUSTIN JIMENEZ AVELLO, agustin.jimenez@upm.es
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

#ifndef _DERIVATIVES_HPP_
#define _DERIVATIVES_HPP_

/**
 * \file derivatives.hpp
 * \brief Calculates several derivatives for a fuzzy system.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */
 
#include <flt/utilities.hpp>

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
		/**
		 * \brief Computes the open-loop jacobian matrix in \c x for the Plant \c S.
		 * 
		 * Obtains the linearization \f$\mathbf{f}(\mathbf{x},\mathbf{u})
		 * \approx \mathbf{F} + \mathbf{A x} + \mathbf{Bu}\f$.
		 * 
		 * @return Returns the dynamic matrix, \f$\mathbf{A}\f$, or an empty array if error.
		 */
		TNT::Array2D<double> jacobian(System &S,
                                      const double *const x,
                                      const double *const u,
                                      TNT::Array2D<double> &B,
                                      TNT::Array1D<double> &F);

		/**
		 * \brief Computes the closed-loop jacobian matrix in \c x for the Plant \c S and the Controller \c C.
		 * 
		 * @return Returns the jacobian matrix or an empty array if error.
		 */
		TNT::Array2D<double> jacobian(System &S,
							          System &C,
							          const double *const x);

		/**	
		 * \brief Computes the approximation of the closed-loop jacobian matrix in \c x for the Plant
		 * \c S and the Controller \c C using finite differences with a step \c h.
		 * 
		 * @return Returns an empty array if error.
		 */
		TNT::Array2D<double> jacobianAprox(System &S,
							               System &C,
							               const double * const x,
							               double h=0.001);

		/**
		 * \brief Computes the approximation of the open-loop jacobian matrix in \c x for the Plant
		 * \c S using finite differences with a step \c h.
		 * 
		 * Obtains the linearization \f$\mathbf{f}(\mathbf{x},\mathbf{u})
		 * \approx \mathbf{F} + \mathbf{A x} + \mathbf{Bu}\f$.
		 * 
		 * @return Returns the dynamic matrix, \f$\mathbf{A}\f$, or an empty array if error.
		 */
		TNT::Array2D<double> jacobianAprox(System &S,
							               const double *const x,
							               const double *const u,
							               TNT::Array2D<double> &B,
							               TNT::Array1D<double> &F,
							               double h=0.001);
		
		/**
		 * \brief Obtains the derivative of a fuzzy model with respect to a consequent.
		 * 
		 * It is assumed that the output of the system is the same as the parameter, otherwise it should be considered derivative equal to 0.
		 * 
		 * @return If \c output, \c rule or \c parameter are out of System limits, the function returns 0 and \c error = 1.
		 */
		double derconseq(System &S,
						 double *x,
						 size_t output,
						 size_t rule,
						 size_t parameter,
						 int &error);

		/**
		 * \brief Obtains the jacobian matrix of a fuzzy model with respect to its consequents.
		 * 
		 * @return Returns an empty array if error.
		 */
		TNT::Array2D<double> derconseq(System &S,
						               double *x);
						               
						               		/**
		 * \brief Obtains the jacobian matrix of a fuzzy model with respect to its affine consequents.
		 * 
		 * @return Returns an empty array if error.
		 */
		TNT::Array2D<double> deraffine(System &S,
						               double *x);
		
		/**
		 * \brief Obtains the derivative of a fuzzy model with respect to a parameter of an antecedent.
		 * 
		 * It is assumed that the output of the system is the same as the parameter, otherwise it should be considered derivative equal to 0.
		 * 
		 * @return If \c output, \c rule or \c parameter are out of System limits, the function returns 0 and \c error = 1.
		 */
		double derantec(System &S,
						double *x,
						size_t input,
						size_t output,
						size_t rule,
						size_t parameter,
						int &error);

		/**
		 * \brief Obtains the jacobian matrix of a fuzzy model with respect to its antecedents.
		 * 
		 * @return Returns an empty array if error.
		 */
		TNT::Array2D<double> derantec(System &S,
						              double *x);
		
		/**
		 * \brief Obtains the jacobian matrix of a fuzzy model with respect to all of its parameters.
		 * 
		 * @return Returns an empty array if error.
		 */		
		TNT::Array2D<double> derfuzzy(System &S,
						              double *x);
} // FLT

#endif
