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

#ifndef _RULE_HPP_
#define _RULE_HPP_

/**
 * \file rule.hpp
 * \brief Defines a Takagi-Sugeno fuzzy Rule.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */

#include <stdio.h>

#include <flt/membership.hpp>

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{

	/**
	 * \class Rule rule.hpp
	 * \brief This class contains methods and attributes common to fuzzy Rules.
	 * 
	 * This class stores a completely general Takagi-Sugeno Rule.
	 * 
	 * \sa You can use the Anymf Membership function to remove a variable from the antecedent of a Rule.
	 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
	 */
	class Rule
	{
		size_t in; ///< Number of inputs of the Rule.
		double *TSK; ///< TSK consequent of the Rule.
		Membership **membership; ///< Membership functions of the Rule.
		inline void clean(void); ///< Cleans the dynamic memory.
	public:
		size_t n_inputs(void) const; ///< Extract the number of inputs of the Rule.

		/**
		 * \brief This function checks the parameters of all Membership functions in the Rule,
		 * and corrects them if possible.
		 * 
		 * @return Returns 0 if no errors, 1 if an error occurs, and -1 if the parameters were corrected.
		 */
		int test(void) const;
		
		/**
		 * \brief Changes a parameter in the TSK consequent.
		 * 
		 * @return Returns 0 if no errors occurred.\n
		 * An error means the input is greater than the number of parameter of the consequent.
		 * 
		 * @note TSK(0) is the affine term, TSK(1) is the term associated to the 1st input variable, ...
		 */
		int changeTSK(double value,
		              size_t input);

		/**
		 * \brief Changes the TSK consequent.
		 */
		void changeTSK(const double* const TSK);
		
		/**
		 * \brief Reads a TSK parameter.
		 * 
		 * @return Returns the selected parameter. If \c input is greater than the number of Rule inputs,
		 * an error is sent to the standard error stream and return 0.
		 */
		double readTSK(size_t input) const;
		
		/**
		 * \brief Reads the TSK consequent.
		 */
		TNT::Array1D<double> readTSK(void) const;
	
		/**
		 * \brief Extracts a pointer to a Membership function.
		 * 
		 * @return Returns NULL if the input is greater than the number of Rule inputs.\n
		 */
		Membership *readFunction(size_t input);
		
		/**
		 * \brief Changes a Membership function of the Rule.
		 * 
		 * @return Returns 0 if no errors occurred.\n
		 * An error means the input is greater than the number of Rule inputs.
		 */
		int changeFunction(const Membership &P,
		                   size_t input);

		/**
		 * \brief Changes the type af a Membership function.
		 * 
		 * @return Returns 0 if no errors occurred.\n
		 * An error means the type of Membership function is not recognized (1) or
		 * \c input is greater than the number of Rule inputs (2).
		 * 
		 * \note The parameters of the new Membership function are not assigned,
		 * you must do so by Membership::edit method.
		 */
		int changeFunction(TYPE_MF type,
		                   size_t input);

		/**
		 * \brief Initializes the Rule for \c in inputs.
		 * 
		 * Empty membership functions are used.
		 * You can change them with \c changeFunction method or using the Membership class methods.
		 */
		void initialize(size_t in);

		/**
		 * \brief  Gets the number of parameters stored in the antecedents of the Rule.
		 * 
		 * @note The function <c>NumberOfConsequents(void)</c> is not necessary because this value is <c>in+1</c>.
		 */
		size_t NumberOfAntecedents(void);

		/**
		 * \brief Sets the parameters of the antecedents of the Rule.
		 * 
		 * @return Returns the result of running the Rule::test() method.
		 */
		int setAntecedents(const double *const parameters);

		/**
		 * \brief Gets the parameters of the antecedents of the Rule.
		 */
		TNT::Array1D<double> getAntecedents(void);

		/**
		 * \brief Sets the consequents of the Rule.
		 */
		void setConsequents(const double *const parameters);

		/**
		 * \brief Gets the consequents of the Rule.
		 */
		TNT::Array1D<double> getConsequents(void);

		/**
		 * \brief Calculates the matching degree of the rule in \c point and, optionally, its derivative.
		 * 
		 * Optionally, this function also calculates the derivative of the degree of
		 * activation of the rule with respect to the input at \c point.\n
		 * dW represents \f$dW / dX\f$, and its length is equal to the number of entries in the rule.
		 * 
		 * The length of \c point must match the number of inputs in the rule.
		 */
		double activation(const double * const point,
		                  double *dW = NULL) const;
			 
		Rule &operator=(const Rule &R);
		Rule();
		Rule(const Rule &R);
		Rule(size_t in);
		~Rule();
	};
}  // FLT namespace

#endif
