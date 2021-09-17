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

#ifndef _MEMBERSHIP_HPP_
#define _MEMBERSHIP_HPP_

/**
 * \file membership.hpp
 * \brief Defines Membership functions.
 * 
 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
 */

#include <math.h>

#include <tnt/tnt_array1d.h> // Template Numerical Toolkit, http://math.nist.gov/tnt
#include <tnt/tnt_array2d.h>

#include <flt/messages.h>

namespace FLT // Use the FTL (Fuzzy Logic Tools) namespace
{
	#define M_EPS 1e-16 ///< Epsilon value
	#define MAX_SIZE_TYPE_NAME 16 ///< Maximum size for the names of the Membership functions.

	/**
	 * \fn int sign(double x)
	 * \brief Implementation of the sign function.
	 * 
	 * @return Returns 0 if x=0, 1 if x>0 and -1 if x<0.
	 */
	inline int sign(double x)
	{
		if (x > 0.0)
			return 1;
		else if (x < 0.0)
			return -1;
		else
			return 0;
	}

	/**
	 * \enum TYPE_MF
	 * \brief Enumeration with the implemented Membership functions.
	 * 
	 * FLT::TYPE_MF is used to determinate the Membership function type, FLT::MF_NAMES are
	 * the names of the Membership functions, and FLT::MF_PARAM_NUMBER are the initial number
	 * of parameters of the Membership functions.

	 * \remark FLT::TYPE_MF, FLT::MF_NAMES and FLT::MF_PARAM_NUMBER must be consistent.
	 * 
	 * \remark \b Important: If the first or last element of FLT::TYPE_MF are changed, you needd to check
	 * membership.cpp (and perhaps others files), because these elements are used to make some comparisons.
	 * 
	 * \var TYPE_MF ANYMF
	 * Any Membership function (FLT::Anymf).
	 * 
	 * \var TYPE_MF CONSTMF
	 * Constant Membership function (FLT::Constmf).
	 *
	 * \var TYPE_MF GAUSSMF
	 * Gaussian Membership function (FLT::Gaussmf).
	 * 
	 * \var TYPE_MF GAUSS2MF
	 * Double gaussian Membership function (FLT::Gauss2mf).
	 *
	 * \var TYPE_MF GBELLMF
	 * Bell Membership function (FLT::GBellmf).
	 * 
	 * \var TYPE_MF PIMF
	 * Pi (S-Z) Membership function (FLT::Pimf).
	 * 
	 * \var TYPE_MF PSIGMF
	 * Product of sigmoidals Membership function (FLT::PSigmf).
	 * 
	 * \var TYPE_MF SMF
	 * S Membership function (FLT::Smf).
	 * 
	 * \var TYPE_MF SIGMF
	 * Sigmoidal Membership function (FLT::Sigmf).
	 * 
	 * \var TYPE_MF SIG2MF
	 * Difference of sigmoidals Membership function (FLT::Sig2mf).
	 * 
	 * \var TYPE_MF TRAPMF
	 * Trapezoidal Membership function (FLT::Trapmf).
	 * 
	 * \var TYPE_MF TRIMF
	 * Triangular Membership function (FLT::Trimf).
	 * 
	 * \var TYPE_MF ZMF
	 * Z Membership function (FLT::Zmf).
	 */
	enum TYPE_MF // Implemented Membership functions
	{
		ANYMF = 0,
		CONSTMF,
		GAUSSMF,
		GAUSS2MF,
		GBELLMF,
		PIMF,
		PSIGMF,
		SMF,
		SIGMF,
		SIG2MF,
		TRAPMF,
		TRIMF,
		ZMF
	};

	/**
	 * \var MF_NAMES
	 * \brief Names of the Membership functions.
	 * 
	 * \remark FLT::TYPE_MF, FLT::MF_NAMES and FLT::MF_PARAM_NUMBER must be consistent.
	 */
	static const char *MF_NAMES[]=
	{// Names of TYPE_MF Membership Functions - UPPERCASE
		"ANYMF",
		"CONSTMF",
		"GAUSSMF",
		"GAUSS2MF",
		"GBELLMF",
		"PIMF",
		"PSIGMF",
		"SMF",
		"SIGMF",
		"SIG2MF",
		"TRAPMF",
		"TRIMF",
		"ZMF"
	};

	/**
	 * \var MF_PARAM_NUMBER
	 * \brief Initial number of parameters of each Membership function in FLT::TYPE_MF.
	 *
	 * \remark These values are used for initialization, in the code is better use the
	 * \c num_params function from Membership class.\n
	 * 
	 * \remark FLT::TYPE_MF, FLT::MF_NAMES and FLT::MF_PARAM_NUMBER must be consistent.
	 */
	static const size_t MF_PARAM_NUMBER[] =
	{
		0, // \var ANYMF
		1, // CONSTMF
		2, // GAUSSMF
		4, // GAUSS2MF
		3, // GBELLMF
		4, // PIMF
		4, // PSIGMF
		2, // SMF
		2, // SIGMF
		4, // SIG2MF
		4, // TRAPMF
		3, // TRIMF
		2  // ZMF
	};

	// Declarations ################################################################

	/**
	 * \class Membership membership.hpp
	 * \brief This class contains methods and attributes common to all Membership functions.
	 * 
	 * This class can store the type of Membership function and all its parameters and evaluates
	 * the function and its derivatives.
	 * 
	 * \sa FLT::TYPE_MF for Membership functions details.
	 * \sa FLT::MF_NAMES for Membership functions names.
	 * \sa FLT::MF_PARAM_NUMBER for initial number of parameters.
	 * \sa For Template Numerical Toolkit (TNT) documentation see http://math.nist.gov/tnt
	 * \sa FLT::createMF virtual constructor.
	 */
	class Membership // -Parent class-
	{
	protected:
		double *param; ///< Vector of parameters that define the Membership function.
		size_t n;  ///< Number of parameters of the Membership function.
		TYPE_MF type_mf; ///< Type of Membership function.
	public:
		/**
		 * \brief Evaluates the Membership function.
		 * 
		 * \remark You must first ensure that the function is not Anymf.
		 */
		virtual double eval(double x) const = 0;
		
		/**
		 * \brief Evaluates the derivative of the Membership function with respect to x, \f$ d\mu(x) / dx \f$
		 * 
		 * \remark You must first ensure that the function is not Anymf.
		 */
		virtual double evalder(double x) const = 0;
		
	/**
	 * \brief Evaluates the derivative of the Membership function with respect to a parameter, \f$ d\mu(x) / dparam \f$
	 * 
	 * \remark You must first ensure that the function is not Anymf.
	 */
		virtual double paramder(size_t parameter,
		                        double x) const = 0;
		
		/**
		 * \brief This function checks parameters in the Membership function, and corrects them if possible.
		 * 
		 * @return Returns 0 if no errors, 1 if an error occurs, and -1 if the parameters were corrected.
		 */
		virtual int test(void) const;
		
		/**
		 * \brief Reads the number of parameters that define the Membership function.
		 * 
		 * We recommend using this function instead of the constant FLT::MF_PARAM_NUMBER,
		 * since it is possible to create Membership functions with a variable number of parameters.
		 */
		size_t num_params(void) const;
		
		/**
		 * \brief Returns the type of Membership function
		 * 
		 * The type of Membership function is represented with the enumeration FLT::TYPE_MF.
		 */
		TYPE_MF type(void) const;
		
		/**
		 * \brief Changes the Membership function type.
		 * 
		 * @return Returns 0 if no errors occurred.\n
		 * An error means the type of Membership function is not recognized.
		 * 
		 * \note The parameters of the new Membership function are not assigned,
		 * you must do so by Membership::edit method.
		 */
		int type(TYPE_MF type_mf);
		
		/**
		 * \brief Changes a parameter.
		 * 
		 * @return Returns 0 if no errors occurred.\n
		 * An error means the index value is greater than the number of function parameters.
		 * 
		 * \note The parameter value is not checked. If you want to check it,
		 * you must use Membership::test method.
		 */
		int edit(size_t index,
		         double value);

		/**
		 * \brief Changes the vector of parameters.
		 * 
		 * \note The values of the parameters are not checked.
		 * If you want to check it, you must use Membership::test method.
		 */
		void edit(const double* const parameters);
		
		/**
		 * \brief Reads a parameter of the Membership function.
		 * 
		 * @return If \c index is greater than the number of function parameters,
		 * an error is sent to the standard error stream and return 0.
		 */
		double read(size_t index) const;
		
		/**
		 * \brief Returns the vector of parameters.
		 */
		TNT::Array1D<double> read(void) const;

		Membership &operator=(const Membership &P);
		Membership(const Membership &P);
		Membership(); 
		~Membership();
};

	/**
	 * \fn Membership *createMF(TYPE_MF t);
	 * \brief Virtual constructor for Membership functions.
	 * 
	 * This function create and return a Membership function of the selected type.
	 */
	Membership *createMF(TYPE_MF t);
	
	// Anymf ########################################################################
	/**
	 * \class Anymf
	 * \brief Any Membership function.
	 * 
	 * Anymf function is used to remove a variable from the antecedent of a Rule. It doesn't have parameters.
	 * 
	 * For example:
	 * <i>If x1 is A1 and x3 is A2 ...</i> is created using <i>x2 is ANYMF</i>.
	 * \remark <b>Anymf never must be evaluated</b>, use \c <i>if</i> to avoid its evaluation.
	 */
	class Anymf:public Membership // This Membership function don't exist in MATLAB©.
	{ // This type is used to initializate Rules
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Anymf(void);
	};

	// Constmf ######################################################################
	/**
	 * \class Constmf
	 * \brief Constant Membership function.
	 * 
	 * Constmf is a constant Membership function with 1 parameter.\n
	 * This function always return its parameter.\n
	 * \f[ \mu [c](x) = c, \; 0 \leq c \leq 1 \f]
	 */
	class Constmf:public Membership // This Membership function don't exist in MATLAB© by default. It is 'constmf.m'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Constmf(double c = 1.0);
	};

	// Gaussmf #####################################################################
	/**
	 * \class Gaussmf
	 * \brief Gaussian Membership function.
	 * 
	 * Gaussmf is a Gaussian Membership function with 2 parameters.\n
	 * \f[ \mu [c,\beta](x) = e^{-\left( \frac{x-c}{\beta}\right)^2 }, \; \beta > 0 \f]	 
	 */
	class Gaussmf:public Membership // In MATLAB© is 'gaussmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		int test(void) const; // If the parameters of the Membership function are valid returns 0
		Gaussmf(double Center = 0.5,
				double Beta = 0.3); // Sigma_MATLAB = Beta/sqrt(2);
	};

	// Gauss2mf ####################################################################
	/**
	 * \class Gauss2mf
	 * \brief Double gaussian Membership function.
	 * 
	 * Gauss2mf is a Double Gaussian Membership function with 4 parameters.\n
	 * It is obtained by the product of 2 Gaussmf functions.\n
	 * \image html gauss2mf.png
	 * \latexonly
	 * \begin{displaymath}
	 * \mu _{Gauss2mf}(c_1, \beta _1, c_2, \beta _2) =
	 * \left\{
	 * 	\begin{array}{ll}
	 * 		\mu _{Gaussmf}(c_1, \beta _1) & \hspace{-11em} \text{if } x < c_1 \leq c_2\vspace{0.3em}\\
	 * 		\mu _{Gaussmf}(c_2, \beta _2) & \hspace{-11em} \text{if } x > c_2, c_1 \leq c_2\vspace{0.3em}\\
	 * 			1 & \hspace{-11em} \text{if } c_1 \leq x \leq c_2 \text{ and } c_1 \leq c_2\vspace{0.3em}\\
	 * 		\mu _{Gaussmf}(c_1, \beta _1) \mu _{Gaussmf}(c_2, \beta _2) \hspace{1em} \text{if } c_1 > c_2
	 * 	\end{array}
	 * \right.
	 * \end{displaymath}
	 * \endlatexonly
	 */
	 class Gauss2mf:public Membership // In MATLAB© is 'gauss2mf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		int test(void) const; // Modify the parameters to be correct
		Gauss2mf(double Center1 = 0.45,
				 double Beta1 = 0.25,
				 double Center2 = 0.55,
				 double Beta2 = 0.25); // Sigma_MATLAB = Beta/sqrt(2);
	};

	// GBellmf ######################################################################
	/**
	 * \class GBellmf
	 * \brief Bell Membership function.
	 * 
	 * GBellmf is a bell Membership function with 3 parameters.\n
	 * \f[ \mu [\alpha ,\beta, c](x) = \frac {1}{1+ \left| \frac{x-c}{\alpha}\right| ^{2 \beta}}, \; \alpha \neq 0 \f]
	 */
	 class GBellmf:public Membership // In MATLAB© is 'gbellmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		int test(void) const; // Modify the parameters to be correct
		GBellmf(double a = 0.25,
		        double b = 2.5,
				double c = 0.5);
	};

	// Pimf ########################################################################
	/**
	 * \class Pimf
	 * \brief Pi (S-Z) Membership function.
	 * 
	 * Pimf is a S-Z Membership function with 4 parameters.\n
	 * It is obtained by the product of a Smf and a Zmf Membership functions.\n
	 * \f[ \mu [a,b,c,d](x) = \mu_S[a,b](x) \cdot \mu_Z[c,d](x) \f]
	 */
	 class Pimf:public Membership // In MATLAB© is 'pimf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Pimf(double a = 0.05,
			 double b = 0.45,
			 double c = 0.55,
			 double d = 0.95);
	};

	// PSigmf ################################################################
	/**
	 * \class PSigmf
	 * \brief Product of sigmoidals Membership function.
	 * 
	 * PSigmf is the product of 2 sigmoidal Membership functions (Sigmf), and it has 4 parameters.\n
	 * \f[ \mu [\alpha_1,c_1,\alpha_2,c_2](x)=\mu_{Sigm}[\alpha_1,c_1](x) \cdot \mu_{Sigm}[\alpha_2,c_2](x) \f]
	 */
	 class PSigmf:public Membership // In MATLAB© is 'psigmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		PSigmf(double a1 = 11,
		       double c1 = 0.25,
			   double a2 = -11,
			   double c2 = 0.75);
	};

	// Smf #########################################################################
	/**
	 * \class Smf
	 * \brief S Membership function.
	 * 
	 * Smf is a S Membership function with 2 parameters.\n
	 * \image html smf.png
	 * \latexonly
	 * \begin{displaymath}
	 * \mu _{Smf}[a, b](x) = 
	 * \left\{
	 * 	\begin{array}{ll}
	 * 		0 & \text{if } x \leq a\\
	 * 		2\left(\frac{x-a}{b-a}\right)^2 & \text{if } a < x \leq \frac{a+b}{2}\vspace{0.3em}\\
	 * 		1 - 2\left(\frac{b-x}{b-a}\right)^2 & \text{if } \frac{a+b}{2} < x < b\vspace{0.3em}\\
	 * 		1 & \text{if } x \geq b
	 * 	\end{array}
	 * \right.
	 * \end{displaymath}
	 * \endlatexonly
	 */
	 class Smf:public Membership // In MATLAB© is 'smf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Smf(double a = 0.05,
		    double b = 0.45);
	};

	// Sigmf #################################################################
	/**
	 * \class Sigmf
	 * \brief Sigmoidal Membership function.
	 * 
	 * Sigmf is a sigmoidal Membership function with 2 parameters.\n
	 * \f[ \mu [\alpha,c](x) = \frac{1}{1+e^{\alpha(c-x)}} \f]
	 */
	 class Sigmf:public Membership // In MATLAB© is 'sigmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Sigmf(double a = 14,
		      double c = 0.25);
	};

	// Sig2mf ################################################################
	/**
	 * \class Sig2mf
	 * \brief Difference of sigmoidals Membership function.
	 * 
	 * Sig2mf is a double sigmoidal Membership function with 4 parameters.\n
	 * It is obtained by the difference between 2 Sigmf functions.\n
	 * \f[ \mu[\alpha_1,c_1,\alpha_2,c_2](x)=\mu_{Sigm}[\alpha_1,c_1](x)-\mu_{Sigm}[\alpha_2,c_2](x) \f]
	 */
	 class Sig2mf:public Membership // In MATLAB© is 'dsigmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Sig2mf(double a1 = 11,
		       double c1 = 0.25,
			   double a2 = 11,
			   double c2 = 0.75);
	};

	// Trapmf ###############################################################
	/**
	 * \class Trapmf
	 * \brief Trapezoidal Membership function.
	 * 
	 * Trapmf is a trapeziodal Membership function with 4 parameters.\n
	 * \image html trapmf.png
	 * \latexonly
	 * \begin{displaymath}
	 * 	\mu _{Trapmf}[a,b,c,d](x) = 
	 * 	\left\{
	 * 		\begin{array}{ll}
	 * 			\frac{x-a}{b-a} & \text{if } a<x<b\vspace{0.3em}\\
	 * 			1 & \text{if } b \leq x \leq c\vspace{0.3em}\\
	 * 			\frac{d-x}{d-c} & \text{if } c < x < d\vspace{0.3em}\\
	 * 			1 & \text{if } x \geq d
	 * 		\end{array}
	 * 	\right.
	 * \end{displaymath}
	 * \endlatexonly
	 * \f[ a < b < c < d \f]
	 */
	 class Trapmf:public Membership // In MATLAB© is 'trapmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		int test(void) const; // Modify the parameters to be correct
		Trapmf(double a = 0.05,
		       double b = 0.45,
			   double c = 0.55,
			   double d = 0.95);
	};

	// Trimf ################################################################
	/**
	 * \class Trimf
	 * \brief Triangular Membership function.
	 * 
	 * Trimf is a triangular Membership function with 3 parameters.\n
	 * \image html trimf.png
	 * \latexonly
	 * \begin{displaymath}
	 * 	\mu _{Trimf}[a,b,c](x) = 
	 * 	\left\{
	 * 		\begin{array}{ll}
	 * 			\frac{x-a}{b-a} & \text{if } a<x \leq b\vspace{0.3em}\\
	 * 			\frac{c-x}{c-b} & \text{if } b < x < c\vspace{0.3em}\\
	 * 			0 & \text{if } x \geq c
	 * 		\end{array}
	 * 	\right.
	 * \end{displaymath}
	 * \endlatexonly
	 * \f[ a < b < c \f]
	 */
	 class Trimf:public Membership // In MATLAB© is 'trimf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		int test(void) const; // Modify the parameters to be correct
		Trimf(double a = 0.0,
		      double b = 0.5,
			  double c = 1.0);
	};

	// Zmf #########################################################################
	/**
	 * \class Zmf
	 * \brief Z Membership function.
	 * 
	 * Zmf is a Z Membership function with 2 parameters.\n
	 * \image html zmf.png
	 * \latexonly
	 * \begin{displaymath}
	 * 	\mu _{Zmf}[a, b](x) = 
	 * 	\left\{
	 * 		\begin{array}{ll}
	 * 			1 & \text{if } x \leq a\\
	 * 			1-2\left(\frac{x-a}{a-b}\right)^2 & \text{if } a < x \leq \frac{a+b}{2}\\
	 * 			2\left(\frac{b-x}{a-b}\right)^2 & \text{if } \frac{a+b}{2} < x < b\\
	 * 			0 & \text{if } x \geq b
	 * 		\end{array}
	 * 	\right.
	 * \end{displaymath}
	 * \endlatexonly
	 */
	 class Zmf:public Membership // In MATLAB© is 'zmf'
	{
	 public:
		double eval(double x) const;
		double evalder(double x) const;
		double paramder(size_t parameter,
		                double x) const;
		Zmf(double a = 0.55,
		    double b = 0.95);
	};
	
} // FLT namespace

#endif
