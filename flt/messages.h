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

#ifndef _MESSAGES_H_
#define _MESSAGES_H_

#include <cstdlib>
#include <iostream>

#ifdef MATLAB_MEX_FILE
	#include <mex.h>
#endif

/**
 * \file messages.h
 * \brief Defines output messages (informational messages, warnings, errors, ...).
 * 
 * Some macros for generating error messages are also defined to facilitate generation
 * console or MATLAB© applications.
 * 
 * This file defines the messages used in all functions according to the following schedule:
 *
 * \b E_* -> Error messages.\n 
 * \b U_* -> User Messages.\n
 * \b O_* -> Other messages.\n
 */

/**
 * \namespace FLT
 * \brief Fuzzy Logic Tools (FLT) namespace. 
 */
namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
/**
 * \def WARNINGMSG(menssage)
 * \brief Generates the given wanring message
 * 
 * This macro is used to generate a warning message.\n
 * The standard error stream is used, but if the code is being compiled in MATLAB©,
 * uses the \c mexErrMsgTxt function from the MATLAB© API.
 */
	#ifdef MATLAB_MEX_FILE
		#define WARNINGMSG(menssage) {mexWarnMsgTxt(menssage);}
	#else
		#define WARNINGMSG(menssage) {std::cerr << menssage << std::endl;}
	#endif

/**
 * \def ERRORMSG(menssage)
 * \brief Generates the given error and return.
 * 
 * This macro is used to generate an error.\n
 * The standard error stream is used, but if the code is being compiled in MATLAB©,
 * uses the \c mexErrMsgTxt function from the MATLAB© API.
 */
	#ifdef MATLAB_MEX_FILE
		#define ERRORMSG(menssage) {mexErrMsgTxt(menssage); return;}
	#else
		#define ERRORMSG(menssage) {std::cerr << menssage << std::endl; return;}
	#endif
	
/**
 * \def ERRORMSGVAL(menssage, value)
 * \brief Generates the given error and return \c value.
 * 
 * This macro is used to generate an error and return a value.\n
 * The standard error stream is used, but if the code is being compiled in MATLAB©,
 * uses the \c mexErrMsgTxt function from the MATLAB© API.
 */
	#ifdef MATLAB_MEX_FILE
		#define ERRORMSGVAL(menssage, value) {mexErrMsgTxt(menssage); return value;}
	#else
		#define ERRORMSGVAL(menssage, value) {std::cerr << menssage << std::endl; return value;}
	#endif

	// Some constants
	
	#define MAX_CHAR_LONG 256 ///< Maximum length of messages strings.
	#define MAX_FILE_NAME 256 ///< Maximum length of file names.
	#define PRECISION 10 ///< Precision for data extraction (for use with f/i/o/stream).
	
	// Warning Messages( W_... ) ###############################################
		#define W_InputsLimitsExceeded "Warning, the input exceeds the limits of the system."
		#define W_OutputsLimitsExceeded "Warning, the output exceeds the limits of the system."
	
	// Error Messages( E_... ) #################################################
		// Model
			#define E_Model "Error reading model."
			#define E_NoSugeno "Model 'type' should be 'sugeno'."
			#define E_NoOutputs "Error, the model has not any Output."
			#define E_NoInputs "Error, the model has not any Input."
			#define E_NoRules "Error, some Outputs of the model have not any Rule."
			#define E_No1Conseq "Only one consequent for rule is allowed."
			#define E_SameIns "All models must have the same inputs."
			#define E_NoCoherent "The Controller is not coherent with the Plant."
			#define E_InNoCoherent "Inputs arguments are not coherent with fuzzy model."
			#define E_BadModel "Undefined or wrong fuzzy model.\nCheck rules, membership functions and theirs parameters"
			#define E_NoCoherent_BadX "The Controller is not coherent with the Plant, or the point is not coherent with the closed loop system."
			#define E_ArgNoValid "Empty and NaN matrices are not valid inputs arguments."
			#define E_AccuracyNoNum "The accuracy should be a number."
			#define E_Acquire "Error reading model or incorrect initial model."
			#define E_Store "Error storing data. The vector size is wrong, or the data is invalid for the membership functions, for example, a zero width in Gaussmf function."
			
		// FIS
			#define E_NoFIS "Fuzzy Model must be a Sugeno FIS structure."
			#define E_CreateFIS "Error at create FIS."
			#define E_NoProd "Implication method must be 'prod'."
			#define E_FISOut "Error writing FIS."
			#define E_InputExceded "Input limit exceeded."
			#define E_OutputExceded "Output limit exceeded."
			#define E_RulesExceded "Rule limit exceeded."
			#define E_ParamExceded "Parameter number exceeded the number of parameters of the membership function."
			
		// Membership functions
			#define E_MFType "Error reading membership function type:"
			#define E_BadMF "Membership function was not identified."
			#define E_MFTypeDef "Error in membership function type."
			#define E_ParamMF "Error reading membership function parameters."
			#define E_NoValidFP "Any parameter of membership function is not valid."
			#define E_MFDef "Error in membership function definition."
			#define E_AnymfUse "ANYMF function must not be evaluated. It is possible that there is an error in the sources, contact with authors, please."

		// Consequents
			#define E_Conseq "It is not posible read consequent type in rule:"
			#define E_ParamConseq "Error reading TSK's consequent parameters."
			#define E_BadConseq "Not valid cosequent."
			
		// Files
			#define E_FileInput "Input file error."
			#define E_FileOutput "Output file error."
			
		// Design
			#define E_Syntax "Syntax error."
			#define E_NumberParam "Number of parameters incorrect."
			#define E_NumberArg "Error in number of arguments."
			#define E_NumberArgIn "Error in inputs arguments"
			#define E_NumberArgOut "Error in outputs arguments."
			#define E_NumberRulesOut "Error, the number of rules must be a vector of size equal to the number of system outputs."
			#define E_Column "The state vector must be a column vector."
			#define E_In_Out_Column "Input and output vector must be column vectors."
			#define E_InNames "Error reading names of inputs."
			#define E_OutNames "Error reading names of outputs."
			#define E_Controller_Init "Controller initialization error."
			#define E_No_Points_Extracted "No points were extracted."
			#define E_Jacobian_Calculation "It happened an error during the calculation of the Jacobian."
			
		// Others
			#define E_ArrayExceded "Array limit exceded."
			#define E_ChangingParameter "Error changing the parameter."
			#define E_WritingRule "Error writing the rule"
			#define E_NoSumProm "Error, 'defuzzMethod' should be 'wtaver'."
			#define E_RuleOutputFileVector "Rules and Outputs should be a file vector."
			#define E_PointCoherent "The point/s should be coherent with fuzzy models."
			#define E_In_No_Scalar "The number of input must be a scalar."
			#define E_Out_No_Scalar "The number of output must be a scalar."
			#define E_R_No_Scalar "The number of rule must be a scalar."
			#define E_Param_No_Scalar "The number of the parameter must be a scalar."
			#define E_Epoch_No_Scalar "The number of epoch must be a scalar."
			#define E_Alfa_No_Scalar "Alfa must be a scalar."
			#define E_Error_No_Scalar "The error must be a scalar."
			#define E_ModelError "Modeling error."
			#define E_H_GE_0 "'h' must be >= 0."
			#define E_N_MESH_P "Number of point to make the mesh must be >1."
			#define E_PRECISION "Precision must be >0."

	// Messages to users ( U_... ) ###############################################
		// FIS
			#define U_FISName "Unnamed"
			#define U_FISInput "Input"
			#define U_FISOutput "Output"
			#define U_RuleAbbreviation "R"
			#define U_InputAbbreviation "-In"
			#define U_OutputAbreviation "-Out"
			#define U_If "IF"
			#define U_Is "is"
			#define U_And "and"
			#define U_OR "or"
			#define U_Then "THEN"
			
		// Others
			#define U_SeeManual "See also Fuzzy Logic Toolbox Manual."
			#define U_SeeHelp "see the help."
			#define U_CheckModel "Check model's definition."
			#define U_MathException "Warning. Math exception generated."
			#define U_Overflow "Happened an overflow."
			#define U_FewData "There are few data for modeling the system."

	// Others messages ( O_... ) ####################################################
		// Inputs
			#define O_Input "Input:"
			#define O_OfIn "of in"
			
		// Outputs
			#define O_OfOut "of out"
			
		// Rules
			#define O_RuleWeigh "The rule's weigh:"
			#define O_RuleConnection "The rule's connection:"
			#define O_AntecedentSize "The size of rule's antecedent"
			#define O_ConsequentSize "The size of rule's consequent"
			#define O_Rule "Rule:"
			
		// Others
			#define O_DifferentOf "is different of"
			#define O_IsNotCorrect "is not correct."
			#define O_MF "Membership Function:"
			#define O_OfPlant "of the Plant"
			#define O_OfController "of Controller"
			#define O_GeneralError "An error occurred. Please, check all parameters or contact with authors."
} // Namespace FLT
#endif
