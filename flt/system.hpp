/*  Copyright (C) 2004-2015
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

#ifndef _SYSTEM_HPP_
#define _SYSTEM_HPP_

/**
 * \file system.hpp
 * \brief Defines a general Takagi-Sugeno fuzzy System.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */

#include <flt/rule.hpp>

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
	#define LOW_DEFAULT_LIMIT -100.0 ///< Default minimum of the universe of discourse.
	#define HIGH_DEFAULT_LIMIT 100.0 ///< Default maximum of the universe of discourse.
	
	/**
	 * \class System system.hpp
	 * \brief This class contains methods and attributes common to fuzzy Systems.
	 * 
	 * This class stores a completely general Takagi-Sugeno fuzzy System.
	 * 
	 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
	 */
	class System
	{
		size_t in;  ///< Number of inputs of the System.
		size_t out; ///< Number of outputs of the System.
		size_t *N;  ///< Number of rules of the System for each output.
		double *low_in; ///< Low limits for the inputs of the System.
		double *low_out; ///< Low limits for the outputs of the System.
		double *high_in; ///< High limits for the inputs of the System.
		double *high_out; ///< High limits for outputs of the System.
		Rule **R; ///< System Rules.
	public:
		size_t inputs(void) const; ///< Returns the number of inputs of the System.
		size_t outputs(void) const; ///< Returns the number of outputs of the System.
		size_t rules(size_t output) const; ///< Returns the number of rules for the output \c output of the System.
		
		TNT::Array1D<size_t> rules(void) const; ///< Returns the number of rules for each output of the System.
		
		TNT::Array1D<double> in_low(void) const; ///< Returns the low limits of the System inputs.
		double in_low(size_t input) const; ///< Returns the low limit of the input \c input.
		
		TNT::Array1D<double> in_high(void) const; ///< Returns the high limits of the System inputs.
		double in_high(size_t input) const; ///< Returns the high limit of the input \c input.
		
		TNT::Array1D<double> out_low(void) const; ///< Returns the low limits of the System outputs.
		double out_low(size_t output) const; ///< Returns the low limit of the output \c output.
		
		TNT::Array1D<double> out_high(void) const; ///< Returns the high limits of the System outputs.
		double out_high(size_t output) const; ///< Returns the high limit of the output \c output.

		void in_low(double *limits); ///< Changes the low limits of the System inputs.
		void in_high(double *limits); ///< Changes the high limits of the System inputs.
		void out_low(double *limits); ///< Changes the low limits of the System outputs.
		void out_high(double *limits); ///< Changes the high limits of the System outputs.

		/**	
		 * \brief Checks if \c point meets the input limits of the system.
		 * 
		 * @return Returns 0 if <c>low_in[i] <= point[i] <= high_in[i]</c>, for <c>i = 1..in</c>.\n
		 * If <c>point[i] < low_in[i]</c> or <c>point[i] > high_in[i]</c>, returns 1 + the lowest value of \c i that breaks limits.
		 */
		int checkInputLimits(const double * const point);

		/**
		 * \brief Checks if \c output meets the output limits of the system.
		 * 
		 * @return Returns 0 if <c>low_out[j] <= output[j] <= high_out[j]</c>, for <c>i = j..out</c>.\n
		 * If <c>output[j] < low_out[j]</c> or <c>output[j] > high_out[j]</c>, returns 1 + the lowest value of \c j that breaks limits.
		 */
		int checkOutputLimits(const double * const output);

		/**
		 * \brief This function checks the parameters of all Membership functions of the System,
		 * and corrects them if possible.
		 * 
		 * @return Returns 0 if no errors, 1 if an error occurs, and -1 if the parameters were corrected.
		 */
		int test(void) const;
		
		/**
		 * \brief Changes a Rule in the System.
		 * 
		 * @return Returns 0 if no errors.\n
		 * An error means that the rule or the output indicated are beyond the limits of the system.
		 */
		int changeRule(const Rule &R,
					   size_t r,
					   size_t output);
		
		/**
		 * \brief Reads a Rule from the System.
		 * 
		 * @return Returns 0 if no errors.\n
		 * An error means that the rule or the output indicated are beyond the limits of the system.
		 */
		int readRule(Rule &R,
		             size_t output,
					 size_t r) const;
		
		/**
		 * \brief Extracts a pointer to a Rule from the System.
		 * 
		 * @return Returns NULL if the rule or the output indicated are beyond the limits of the system.
		 */
		Rule* readRule(size_t output,
		               size_t r);
		
		/**
		 * \brief Initializes the System with \c in inputs, \c out outputs and \c N[out] rules for each output.
		 */
		void initialize(size_t in,
		                size_t out,
						size_t *N);

		/**
		 * \brief Gets the number of parameters stored in the antecedents of the System.
		 */
		size_t NumberOfAntecedents(void);

		/**
		 * \brief Sets all the parameters of the antecedents of the System.
		 * 
		 * The antecedents must be sorted for each output, for each rule.
		 * 
		 * @return Returns the result of running the System::test() method. 
		 */
		int setAntecedents(const double *const parameters);
	
		/**
		 * \brief Gets all the parameters of the antecedents of the System.
		 * 
		 * The antecedents must be sorted for each output, for each rule.
		 */
		TNT::Array1D<double> getAntecedents(void);

		/**
		 * \brief Gets the number of parameters stored in the consequents of the System.
		 */
		size_t NumberOfConsequents(void);
		
		/**
		 * \brief Sets all the consequets of the System.
		 * 
		 * The consequets must be sorted for each output, for each rule, the affine term and the consequents of each input.
		 */
		void setConsequents(const double *const parameters);
	
		/**
		 * \brief Gets the consequets of the System.
		 * 
		 * The consequets are sorted for each output, for each rule, the affine term and the consequents of each input.
		 */
		TNT::Array1D<double> getConsequents(void);

		/**
		 * \brief Evaluates the System for input \c point.
		 * 
		 * If \c point or \c output are beyond the limits of the system,
		 * a warning message is sent to the standard error stream.
		 */
		TNT::Array1D<double> evaluate(const double * const point);

		System & operator=(const System &S);
		System();
		System(const System &S);
		System(size_t in,
		       size_t out,
			   size_t *N);
		~System();
	};
}  // FLT namespace

#endif
