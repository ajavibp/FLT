#include <flt/membership.hpp>

using namespace FLT;
using namespace TNT;
using namespace std;

// Implementations #############################################################

// Membership  -Parent Class- ##################################################

int Membership::test(void) const
{
	return 0;
}

size_t Membership::num_params(void) const
{
	return n;
}

TYPE_MF Membership::type(void) const
{
	return type_mf;
}

int Membership::type(TYPE_MF type_mf)
{
    if (type_mf != this->type_mf)
    {
        if (type_mf>ZMF || type_mf<ANYMF) // The fisrt and last value of TYPE_MF
            return 1;
        if (param != NULL)
            delete [] param;
        n = MF_PARAM_NUMBER[type_mf];
        this->type_mf = type_mf;
        if (n==0)
        	param = NULL;
    	else
	        param = new double[n];
    }
    return 0;
}

int Membership::edit(size_t index,double value)
{
    if (index >= n)
        return 1;
    param[index] = value;
    return 0;
}

void Membership::edit(const double* const parameters)
{
	for (size_t i=0; i<n; i++)
		param[i] = parameters[i];
}

double Membership::read(size_t index) const
{
    if (index >= n)
		ERRORMSGVAL(E_ArrayExceded, 0.0)
    return param[index];
}

Array1D<double> Membership::read(void) const
{
	Array1D<double> parameters(n, param);
	return parameters;
}

Membership & Membership::operator=(const Membership &P)
{
    if (this == &P)
        return *this;

    if (param != NULL)
        delete [] param;

    n = P.num_params();
	type_mf = P.type();

	if (n)
	{
		param = new double[n];
		for (size_t i=0; i<n; i++)
		    param[i] = P.read(i);
	}

    return *this;
}

Membership::Membership(const Membership &P)
{
    n = P.num_params();
    type_mf = P.type();

	if (n)
	{
		param = new double[n];
		for (size_t i=0; i<n; i++)
		    param[i] = P.read(i);
	}
}

Membership::Membership()
{
    n = MF_PARAM_NUMBER[ANYMF];
    type_mf = ANYMF;
    param = NULL;
}

Membership::~Membership()
{
    if (param != NULL)
    {
        delete [] param;
        param = NULL;
    }
}

// Anymf ########################################################################

double Anymf::eval(double x) const
{
	ERRORMSGVAL(E_AnymfUse, 0.0)
}

double Anymf::evalder(double x) const
{
	ERRORMSGVAL(E_AnymfUse, 0.0)
}

double Anymf::paramder(size_t param, double x) const
{
	ERRORMSGVAL(E_AnymfUse, 0.0)
}

Anymf::Anymf(void)
{
    n = MF_PARAM_NUMBER[ANYMF];
    type_mf = ANYMF;
    param = NULL;
}

// Constmf ####################################################################

double Constmf::eval(double x) const
{
	return param[0];
}

double Constmf::evalder(double x) const
{
	return 0.0;
}

double Constmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	return 1.0;
}

Constmf::Constmf(double C)
{
    n = MF_PARAM_NUMBER[CONSTMF];
    type_mf = CONSTMF;
    param = new double[n];
    param[0] = C;
}

// Gaussmf #####################################################################

double Gaussmf::eval(double x) const
{
    double temp = (x - param[0]) / param[1];
    return exp(-temp * temp);
}

double Gaussmf::evalder(double x) const
{
    return -2.0 * eval(x) * (x - param[0]) / pow(param[1], 2.0);
}

double Gaussmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double df_p = 2.0 * eval(x) * (x - param[0]) / pow(param[1], 2.0);
	if (parameter==1)
		df_p *= (x-param[0])/param[1];

	return df_p;
}

int Gaussmf::test(void) const
{
    if (param[1] == 0.0)
        return 1;
	int error = 0;
    if (param[1] < 0.0)
	{
        param[1] = -param[1];
		error = -1;
	}

    return error;
}

Gaussmf::Gaussmf(double Center, double Beta)
{
    n = MF_PARAM_NUMBER[GAUSSMF];;
    type_mf = GAUSSMF;
    param = new double[n];
    param[0] = Center;
    param[1] = Beta;
}

// Gauss2mf ####################################################################

double Gauss2mf::eval(double x) const
{
    double y1, y2;

    if (x<param[0])
	{
		Gaussmf G1(param[0], param[1]);
        y1 = G1.eval(x);
	}
    else
        y1 = 1.0;

    if (x>param[2])
	{
		Gaussmf G2(param[2], param[3]);
		y2 = G2.eval(x);
	}
    else
        y2 = 1.0;

    return y1 * y2;
}

double Gauss2mf::evalder(double x) const
{
	Gaussmf G1(param[0], param[1]);
	Gaussmf G2(param[2], param[3]);

	double y1, y2, dy1, dy2;

	if (x < param[0])
	{
		y1 = G1.eval(x);
		dy1 = G1.evalder(x);
	}
	else
	{
		y1 = 1.0;
		if (x == param[0])
			dy1 = G1.evalder(x) / 2.0;
		else
			dy1 = 0.0;
	}

	if (x > param[2])
	{
		y2 = G2.eval(x);
		dy2 = G2.evalder(x);
	}
	else
	{
		y2 = 1.0;
		if (x == param[2])
			dy2 = G2.evalder(x) / 2.0;
		else
			dy2 = 0.0;
	}

	return y1*dy2 + y2*dy1;
}

double Gauss2mf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	Gaussmf G1(param[0], param[1]);
	Gaussmf G2(param[2], param[3]);

	if (parameter<2)
	{
		if (x < param[0])
			return G1.paramder(parameter, x) * G2.eval(x);
		else if (x > param[0])
			return 0.0;
		else
			return G1.paramder(parameter, x) * G2.eval(x) / 2.0;
	}
	else
	{
		if (x > param[2])
			return G1.eval(x) * G2.paramder(parameter-2, x);
		else if (x < param[2])
			return 0.0;
		else
			return G1.eval(x) * G2.paramder(parameter-2, x) / 2.0;
	}
}

int Gauss2mf::test(void) const
{
    if (param[1] == 0.0 || param[3] == 0.0)
        return 1;

	int error = 0;
    if (param[1] < 0.0)
	{
        param[1] = -param[1];
		error = -1;
	}

    if (param[3] < 0.0)
	{
        param[3] = -param[3];
		error = -1;
	}

    return error;
}

Gauss2mf::Gauss2mf(double Center1, double Beta1, double Center2, double Beta2)
{
    n = MF_PARAM_NUMBER[GAUSS2MF];
    type_mf = GAUSS2MF;
    param = new double[n];
    param[0] = Center1;
    param[1] = Beta1;
    param[2] = Center2;
    param[3] = Beta2;
}

// GBellmf ######################################################################

double GBellmf::eval(double x) const
{
    return 1.0/(1.0+pow(fabs((x-param[2])/param[0]),2.0*param[1]));
}

double GBellmf::evalder(double x) const
{
    return -2.0*param[1]*eval(x)*eval(x)*pow(fabs((x-param[2])/param[0]),(2.0*param[1]-1.0))*sign((x-param[2])/param[0])/param[0];
}

double GBellmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double df_p=0;
	if (parameter == 0)
		df_p = 2.0*pow(eval(x),2.0)*param[1]*sign(param[0])*fabs(x-param[2])*pow(fabs((x-param[2])/param[0]),2.0*param[1]-1.0)/pow(param[0],2.0);
	else if(parameter==1)
		df_p = -2.0*pow(eval(x),2.0)*log(fabs((x-param[2])/param[0]))*pow(fabs((x-param[2])/param[0]),2.0*param[1]);
	else
		df_p = 2.0*pow(eval(x),2.0)*param[1]*pow(fabs(x-param[2]),2.0*param[1]-2.0)*(x-param[2])/pow(fabs(param[0]),2.0*param[1]);

	return df_p;
}

int GBellmf::test(void) const
{
    if (param[0] == 0.0)
        return 1;
    return 0;
}

GBellmf::GBellmf(double a, double b, double c)
{
    n = MF_PARAM_NUMBER[GBELLMF];
    type_mf = GBELLMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
    param[2] = c;
}

// Pimf ########################################################################

double Pimf::eval(double x) const
{
	Smf S(param[0], param[1]);
	Zmf Z(param[2], param[3]);

	return S.eval(x) * Z.eval(x);
}

double Pimf::evalder(double x) const
{
	Smf S(param[0], param[1]);
	Zmf Z(param[2], param[3]);

	return S.eval(x)*Z.evalder(x) + S.evalder(x)*Z.eval(x);
}

double Pimf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	Smf S(param[0], param[1]);
	Zmf Z(param[2], param[3]);

	if (parameter<2)
		return S.paramder(parameter, x) * Z.eval(x);
	else
		return S.eval(x) * Z.paramder(parameter-2, x);
}

Pimf::Pimf(double a, double b, double c, double d)
{
    n = MF_PARAM_NUMBER[PIMF];
    type_mf = PIMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
    param[2] = c;
    param[3] = d;
}

// PSigmf ################################################################

double PSigmf::eval(double x) const
{
    double s1 = 1.0/(1.0+exp(param[0]*(param[1]-x)));
    double s2 = 1.0/(1.0+exp(param[2]*(param[3]-x)));

    return s1*s2;
}

double PSigmf::evalder(double x) const
{
	Sigmf S1(param[0], param[1]);
	Sigmf S2(param[2], param[3]);

	return S1.eval(x)*S2.evalder(x) + S1.evalder(x)*S2.eval(x);
}

double PSigmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	Sigmf S1(param[0], param[1]);
	Sigmf S2(param[2], param[3]);

	if (parameter<2)
		return S1.paramder(parameter, x) * S2.eval(x);
	else
		return S1.eval(x) * S2.paramder(parameter-2, x);
}

PSigmf::PSigmf(double a1, double c1, double a2, double c2)
{
    n = MF_PARAM_NUMBER[PSIGMF];
    type_mf = PSIGMF;
    param = new double[n];
    param[0] = a1;
    param[1] = c1;
    param[2] = a2;
    param[3] = c2;
}

// Smf #########################################################################

double Smf::eval(double x) const
{
    if (param[0] >= param[1])
        return (x >= (param[0]+param[1])/2.0);
    else if (x <= param[0])
        return 0;
    else if (x > param[0] && x <= ((param[0]+param[1])/2.0))
        return 2.0*pow(((x-param[0])/(param[1]-param[0])),2.0);
    else if (x >= param[1])
        return 1.0;
    else
        return 1.0-2.0*pow(((param[1]-x)/(param[1]-param[0])),2.0);
}

double Smf::evalder(double x) const
{
    if (x <= param[0] || x >= param[1])
        return 0.0;
    else if (x > param[0] && x <= ((param[0]+param[1])/2.0))
        return 4.0*(x-param[0])/pow(param[1]-param[0],2.0);
    else
        return 4.0*(param[1]-x)/pow(param[1]-param[0],2.0);
}

double Smf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double a = param[0];
	double b = param[1];
	if (x<a || x>b)
		return 0.0;

	double p = (a+b)/2.0;
	if (parameter == 0)
	{
		if (x<p)
			return 4.0*(x-a)*(x-b)/pow(b-a,3.0);
		else if (x>p)
			return 4.0*pow(x-b,2.0)/pow(a-b,3.0);
		else if (x==a)
			return 2.0*(x-a)*(x-b)/pow(b-a,3.0);
		else if (x==b)
			return 2.0*pow(x-b,2.0)/pow(a-b,3.0);
		else // x==p
			return 2.0*(x-b)/pow(a-b,2.0);
	}
	else
	{
		if (x<p)
			return 4.0*pow(a-x,2.0)/pow(a-b,3.0);
		else if (x>p)
			return 4.0*(x-a)*(b-x)/pow(a-b,3.0);
		else if (x==a)
			return 2.0*pow(a-x,2.0)/pow(a-b,3.0);
		else if (x==b)
			return 2.0*(x-a)*(b-x)/pow(a-b,3.0);
		else // x==p
			return 2.0*(a-x)/pow(a-b,2.0);
	}
}

Smf::Smf(double a, double b)
{
    n = MF_PARAM_NUMBER[SMF];
    type_mf = SMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
}

// Sigmf #################################################################

double Sigmf::eval(double x) const
{
    return 1.0/(1.0+exp(param[0]*(param[1]-x)));
}

double Sigmf::evalder(double x) const
{
    return param[0]*eval(x)*eval(x)*exp(param[0]*(param[1]-x));
}

double Sigmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	if (parameter==0)
		return (x-param[1])*eval(x)*eval(x)*exp(param[0]*(param[1]-x));
	else
		return -param[0]*eval(x)*eval(x)*exp(param[0]*(param[1]-x));
}

Sigmf::Sigmf(double a, double c)
{
    n = MF_PARAM_NUMBER[SIGMF];
    type_mf = SIGMF;
    param = new double[n];
    param[0] = a;
    param[1] = c;
}

// Sig2mf ################################################################

double Sig2mf::eval(double x) const
{
    double s1 = 1.0/(1.0+exp(param[0]*(param[1]-x)));
    double s2 = 1.0/(1.0+exp(param[2]*(param[3]-x)));
    return s1-s2;
}

double Sig2mf::evalder(double x) const
{
    double ds1 = param[0]*eval(x)*eval(x)*exp(param[0]*(param[1]-x));
    double ds2 = param[2]*eval(x)*eval(x)*exp(param[2]*(param[3]-x));
    return ds1-ds2;
}

double Sig2mf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	Sigmf S1(param[0], param[1]);
	Sigmf S2(param[2], param[3]);

	return S1.paramder(parameter, x) - S2.paramder(parameter-2, x);
}

Sig2mf::Sig2mf(double a1, double c1, double a2, double c2)
{
    n = MF_PARAM_NUMBER[SIG2MF];
    type_mf = SIG2MF;
    param = new double[n];
    param[0] = a1;
    param[1] = c1;
    param[2] = a2;
    param[3] = c2;
}

// TRAPMF ###############################################################

double Trapmf::eval(double x) const
{
    if (x>param[0] && x<param[1])
        return (x-param[0])/(param[1]-param[0]);
    else if (x>param[2] && x<param[3])
        return (param[3]-x)/(param[3]-param[2]);
    else if (x>=param[1] && x<=param[2])
        return 1.0;
    else
        return 0.0;
}

double Trapmf::evalder(double x) const
{ // On singular points df(x) = [df(x+) + dx(x-)] / 2
    if (x>param[0] && x<param[1])
        return 1.0/(param[1]-param[0]);
    else if (x>param[2] && x<param[3])
        return 1.0/(param[2]-param[3]);
    else if (x==param[0] || x==param[1])
        return 0.5/(param[1]-param[0]);
    else if (x==param[2] || x==param[3])
        return 0.5/(param[2]-param[3]);
    else
        return 0.0;
}

double Trapmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double df_p=0.0;
    if (x>=param[0] && x<=param[1])
	{
		if (parameter==0)
			df_p = (x-param[1]) / pow(param[1]-param[0],2.0);
		else if (parameter==1)
			df_p = (param[0]-x) / pow(param[1]-param[0],2.0);
		if (x==param[0] || x==param[1])
			df_p /= 2.0;
	}
    else if (x>=param[2] && x<=param[3])
	{
		if (parameter==2)
			df_p = (param[3]-x) / pow(param[3]-param[2],2.0);
		else if (parameter==3)
			df_p = (x-param[2]) / pow(param[3]-param[2],2.0);
		if (x==param[2] || x==param[3])
			df_p /= 2.0;
	}

	return df_p;
}

int Trapmf::test(void) const
{
	int error = 0;
	for (size_t i=1;i<n;i++)
	{
		if (param[i] <= param[i-1])
		{
			param[i] = param[i-1] + M_EPS;
			error = -1;
		}
	}

    return error;
}

Trapmf::Trapmf(double a, double b, double c, double d)
{
    n = MF_PARAM_NUMBER[TRAPMF];
    type_mf = TRAPMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
    param[2] = c;
    param[3] = d;
}

// Trimf #######################################################################

double Trimf::eval(double x) const
{
    if (x>param[0] && x<=param[1])
        return (x-param[0])/(param[1]-param[0]);
    else if (x>param[1] && x<param[2])
        return (param[2]-x)/(param[2]-param[1]);
    else
        return 0.0;
}

double Trimf::evalder(double x) const
{// On singular points df(x) = [df(x+) + dx(x-)] / 2
    if (x>param[0] && x<param[1])
        return 1.0/(param[1]-param[0]);
    else if (x>param[1] && x<param[2])
        return 1.0/(param[1]-param[2]);
    else if (x==param[0])
        return 0.5/(param[1]-param[0]);
    else if (x==param[1])
        return (1.0/(param[1]-param[0])+1/(param[1]-param[2]))/2.0;
    else if (x==param[2])
        return 0.5/(param[1]-param[2]);
    else
        return 0.0;
}

double Trimf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double df_p=0.0;
	if (x>=param[0] && x<param[1])
	{
		if (parameter == 0)
			df_p = (x-param[1]) / pow(param[1]-param[0],2.0);
		else if (parameter == 1)
			df_p = (param[0]-x) / pow(param[1]-param[0],2.0);
		if (x==param[0])
			df_p /= 2.0;
	}
    else if (x>param[1] && x<=param[2])
	{
		if (parameter == 1)
			df_p = (param[2]-x) / pow(param[2]-param[1],2.0);
		else if (parameter == 2)
			df_p = (x-param[1]) / pow(param[2]-param[1],2.0);
		if (x==param[2])
			df_p /= 2.0;
	}
	else if (x==param[1])
	{
		if (parameter==1)
			df_p = ( (param[0]-x) / pow(param[1]-param[0],2) + (param[2]-x) / pow(param[2]-param[1],2.0) ) / 2.0;
		else
			df_p /= 2.0;
	}

	return df_p;
}

int Trimf::test(void) const
{
	int error = 0;
	for (size_t i=1;i<n;i++)
	{
		if (param[i] <= param[i-1])
		{
			param[i] = param[i-1] + M_EPS;
			error = -1;
		}
	}

    return error;
}

Trimf::Trimf(double a, double b, double c)
{
    n = MF_PARAM_NUMBER[TRIMF];
    type_mf = TRIMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
    param[2] = c;
}

// Zmf #########################################################################

double Zmf::eval(double x) const
{
    if (param[0]>=param[1])
        return (x<=(param[0]+param[1])/2.0);
    else if (x<=param[0])
        return 1.0;
    else if (x>=param[1])
        return 0.0;
    else if (x>param[0] && x<=((param[0]+param[1])/2.0))
        return 1.0-2.0*pow((x-param[0])/(param[0]-param[1]),2.0);
    else
        return 2.0*pow((param[1]-x)/(param[0]-param[1]),2.0);
}

double Zmf::evalder(double x) const
{
    if (x<=param[0] || x>=param[1])
        return 0.0;
    else if (x>param[0] && x<=((param[0]+param[1])/2.0))
        return -4.0*(x-param[0])/pow(param[0]-param[1],2.0);
    else
        return -4.0*(param[1]-x)/pow(param[0]-param[1],2.0);
}

double Zmf::paramder(size_t parameter, double x) const
{
	if (parameter > n)
		ERRORMSGVAL(E_ParamExceded, 0.0);

	double a = param[0];
	double b = param[1];
	if (x<a || x>b)
		return 0.0;

	double p = (a+b)/2.0;
	if (parameter == 0)
	{
		if (x<p)
			return 4.0*(x-a)*(x-b)/pow(a-b,3.0);
		else if (x>p)
			return 4.0*pow(b-x,2.0)/pow(b-a,3.0);
		else if (x==a)
			return 2.0*(x-a)*(x-b)/pow(a-b,3.0);
		else if (x==b)
			return 2.0*pow(b-x,2.0)/pow(b-a,3.0);
		else // x==p
			return 2.0*(b-x)/pow(a-b,2.0);
	}
	else
	{
		if (x<p)
			return 4.0*pow(a-x,2.0)/pow(b-a,3.0);
		else if (x>p)
			return 4.0*(a-x)*(b-x)/pow(a-b,3.0);
		else if (x==a)
			return 2.0*pow(a-x,2.0)/pow(b-a,3.0);
		else if (x==b)
			return 2.0*(a-x)*(b-x)/pow(a-b,3.0);
		else // x==p
			return 2.0*(x-a)/pow(a-b,2.0);
	}
}

Zmf::Zmf(double a, double b)
{
    n = MF_PARAM_NUMBER[ZMF];
    type_mf = ZMF;
    param = new double[n];
    param[0] = a;
    param[1] = b;
}

// Virtual Constructor #########################################################

Membership *FLT::createMF(TYPE_MF t)
{
    switch(t)
    {
        case ANYMF:
            return new Anymf;
        case CONSTMF:
        	return new Constmf;
        case GAUSSMF:
            return new Gaussmf;
        case GAUSS2MF:
            return new Gauss2mf;
        case GBELLMF:
            return new GBellmf;
        case PIMF:
            return new Pimf;
        case PSIGMF:
            return new PSigmf;
        case SMF:
            return new Smf;
        case SIGMF:
            return new Sigmf;
        case SIG2MF:
            return new Sig2mf;
        case TRAPMF:
            return new Trapmf;
        case TRIMF:
            return new Trimf;
        case ZMF:
            return new Zmf;
        default:
		{
			cerr << E_BadMF << endl;
            return NULL;
		}
    }
}
