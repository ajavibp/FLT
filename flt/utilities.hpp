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

#ifndef _UTILITIES_HPP_
#define _UTILITIES_HPP_

/**
 * \file utilities.hpp
 * \brief Implements useful functions for fuzzy systems.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */

#include <stdlib.h>
#include <string.h>

#include <flt/system.hpp>

#ifndef strcmpi
    #define strcmpi strcasecmp // Make the change to compile in Windows(MinGW) & GNU/Linux(gcc/g++)
#endif

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
	/**
	 * \brief Saves the fuzzy System \c S in a text file.
	 * 
	 * @return Returns 0 if no errors.\n
	 * \sa See TXT2System to see the template for the TXT file.
	 */
	int System2TXT(System &S,
				   const char *file);
	
	/**
	 * \brief Reads a fuzzy System from a text file.
	 * 
	 * @return Returns an empty System if error.
	 *
	 * \par Template for the TXT file
	 * The text in italics is only helpful, should not appear in the template.\n
	 * \n
	 * <c><b>'Number of inputs' 'Number of outputs'\n
	 * 'Number of rules for the 1st ouput' 'Number of rules for the 2nd ouput'</b></c> ...\n
	 * \n
	 * <i>Rules: For each output, for each rule in this output:\n
	 * Antecedents. For each input:</i>\n
	 * <c><b>'Type membership function' 'Parameters separated by spaces'</b></c>\n
	 * <i>Consequents:</i>\n
	 * <c><b>'Affine term'</b></c> <i>and for each input</i> <c><b>'consequents terms separated by spaces'\n
	 * \n
	 * 'lower limit for input 1' 'upper limit for input 1'\n
	 * 'lower limit for input 2' 'upper limit for input 2'</b></c>\n
	 * ... <i>for all the inputs</i>\n
	 * \n
	 * <c><b>'lower limit for output 1' 'upper limit for output 1'\n
	 * 'lower limit for output 2' 'upper limit for output 2'</b></c>\n
	 * ... <i>for all the outputs</i></c>
	 */
	System TXT2System(const char *file);

	/**
	 * \brief Writes a fuzzy system in its linguistic form.
	 * 
	 * "IF temperature is high and pressure is ... THEN valve is ..."
	 * 
	 * You can specify the name of the variables and the number of decimals
	 * that will be used to represent the system.
	 */	
	int printSystem(const char *file,
					System &S,
					char *inputs[] = NULL,
					char *outputs[] = NULL,
					int accuracy = 10);

	/**
	 * \brief Evaluates a closed loop fuzzy system in the given \c point.
	 * 
	 * \c S is the plant and \c C is a fuzzy controller.
	 * 
	 * @return Returns an empty array if error.
	 */
	TNT::Array1D<double> evaluate(const double * const point,
								  System &S,
								  System &C);
	
	/**
	 * \brief Creates a fuzzy Controller from a given plant.
	 * 
	 * The fuzzy controller will have all the rules that the plant owns in each of its outputs.\n
	 * Consequents all initialized to 0.
	 * 
	 * @return Returns an empty System if error.
	 */
	System initialController(System &Plant);
	
	/**
	 * \brief Extracts a subsystem from a fuzzy System.
	 * 
	 * @return Returns an empty System if error.
	 */
	System subSystem(System &S,
					 size_t nrules,
					 size_t *outputs,
					 size_t *rules);

	/**
	 * \brief Extracts the representative points of a fuzzy system.
	 * 
	 * @param S Is a fuzzy System.\n 
	 * @param numMeshPoints Is used in Memberships functions without significative points, like ANYMF, CONSTMF, SMF, ...\n
	 * @param precision Is used to remove similar points.\n
	 * @param addMesh If \c addMesh is true, a mesh of \c numMeshPoints for each dimmension will be added to System points.\n
	 * @param onlyStatePoints If \c onlyStatePoints is true, this function only extracts the state vector points, not the control vector points.
	 */
	TNT::Array2D<double> extractPoints(System &S,
	                                   unsigned int numMeshPoints = 5,
									   double precision = 1E-3,
									   bool addMesh = false,
									   bool onlyStatePoints = true);

}  // FLT namespace

#endif
