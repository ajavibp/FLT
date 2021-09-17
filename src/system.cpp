#include <flt/system.hpp>

using namespace FLT;
using namespace TNT;
using namespace std;

size_t System::inputs(void) const
{
	return in;
}

size_t System::outputs(void) const
{
	return out;
}

size_t System::rules(size_t output) const
{
    if (output >= out)
		ERRORMSGVAL(E_OutputExceded, 0)
    return N[output];
}

Array1D<size_t> System::rules(void) const
{
	Array1D<size_t> rules(out);
	for (size_t j=0;j<out;j++)
		rules[j] = N[j];
	return rules;
}

Array1D<double> System::in_low(void) const
{
	Array1D<double> limits(in);
	for (size_t i=0;i<in;i++)
		limits[i] = low_in[i];
	return limits;
}

double System::in_low(size_t input) const
{
    if (input >= in)
		ERRORMSGVAL(E_InputExceded, LOW_DEFAULT_LIMIT)
    return low_in[input];
}

Array1D<double> System::in_high(void) const
{
	Array1D<double> limits(in);
	for (size_t i=0;i<in;i++)
		limits[i] = high_in[i];
	return limits;
}

double System::in_high(size_t input) const
{
    if (input >= in)
		ERRORMSGVAL(E_InputExceded, HIGH_DEFAULT_LIMIT)
    return high_in[input];
}

Array1D<double> System::out_low(void) const
{
	Array1D<double> limits(out);
	for (size_t j=0;j<out;j++)
		limits[j] = low_out[j];
	return limits;
}

double System::out_low(size_t output) const
{
    if (output >= out)
		ERRORMSGVAL(E_OutputExceded, LOW_DEFAULT_LIMIT)
    return low_out[output];
}

Array1D<double> System::out_high(void) const
{
	Array1D<double> limits(out);
	for (size_t j=0;j<out;j++)
		limits[j] = high_out[j];
	return limits;
}

double System::out_high(size_t output) const
{
    if (output >= out)
		ERRORMSGVAL(E_OutputExceded, HIGH_DEFAULT_LIMIT)
    return high_out[output];
}

void System::in_low(double *limits)
{
	for (size_t i=0;i<in;i++)
		low_in[i] = limits[i];
}

void System::in_high(double *limits)
{
	for (size_t i=0;i<in;i++)
		high_in[i] = limits[i];
}

void System::out_low(double *limits)
{
	for (size_t j=0;j<out;j++)
		low_out[j] = limits[j];
}

void System::out_high(double *limits)
{
	for (size_t j=0;j<out;j++)
		high_out[j] = limits[j];
}

int System::checkInputLimits(const double * const point)
{
	for (int i=0; i<in; i++)
		if (point[i] < low_in[i] || point[i] > high_in[i])
			return i+1;
	return 0;
}

int System::checkOutputLimits(const double * const output)
{
	for (int j=0; j<out; j++)
		if (output[j] < low_out[j] || output[j] > high_out[j])
			return j+1;
	return 0;
}

int System::test(void) const
{
	int error = 0;	
    if (N==NULL || R==NULL || out==0 || in==0)
        return 1;
    else
    {
        for (size_t j=0;j<out;j++)
        {
            if (R[j] == NULL)
                return 1;
			int test = R[j]->test();
			if (test > 0)
				return 1;
            if (test < 0)
                error = test;
        }
    }
    return error;
}

int System::changeRule(const Rule &R, size_t r, size_t output)
{
    if (output >= out || r >= N[output])
        return 1;
    this->R[output][r] = R;
    return 0;
}

int System::readRule(Rule &R, size_t output, size_t r) const
{
    if (output >= out || r >= N[output])
        return 1;
    R = this->R[output][r];
    return 0;
}

Rule * System::readRule(size_t output, size_t r)
{
    if (output >= out)
		ERRORMSGVAL(E_OutputExceded, NULL)
    if (r >= N[output])
		ERRORMSGVAL(E_RulesExceded, NULL)
    return &R[output][r];
}

void System::initialize(size_t in, size_t out, size_t *N)
{
    size_t r,i,j;
	if (in==0 || out ==0)
		return;
	if (this->N)
		delete [] (this->N);
	this->N = new size_t[out];
	if (R != NULL)
	{
		for (j=0;j<(this->out);j++)
			if (R[j] != NULL)
				delete [] R[j];
		delete [] R;
	}
	if (low_in != NULL)
		delete [] low_in;
	if (low_out != NULL)
		delete [] low_out;
	if (high_in != NULL)
		delete [] high_in;
	if (high_out != NULL)
		delete [] high_out;
    this->in = in;
    this->out = out;
    R = new Rule*[out];
    low_in = new double [in];
	high_in = new double [in];
	low_out = new double [out];
	high_out = new double [out];
	for (i=0;i<in;i++)
	{
		low_in[i] = LOW_DEFAULT_LIMIT;
		high_in[i] = HIGH_DEFAULT_LIMIT;
	}
    for (j=0;j<out;j++)
    {
		low_out[j] = LOW_DEFAULT_LIMIT;
		high_out[j] = HIGH_DEFAULT_LIMIT;
        this->N[j] = N[j];
        R[j] = new Rule[N[j]];
        for (r=0;r<N[j];r++)
            R[j][r].initialize(in);
    }
}

size_t System::NumberOfAntecedents(void)
{
    size_t size_data = 0;

    for (size_t j=0;j<out;j++)
        for (size_t r=0;r<N[j];r++)
            for (size_t i=0;i<in;i++)
                size_data += R[j][r].readFunction(i)->num_params();
				
    return size_data;
}

int System::setAntecedents(const double *const parameters)
{
    Membership *P;
    size_t index = 0;
    for (size_t j=0;j<out;j++)
    {
        for (size_t r=0;r<N[j];r++)
        {
            for (size_t i=0;i<in;i++)
            {
                P = R[j][r].readFunction(i);
                size_t size_data = P->num_params();
                for (size_t p=0;p<size_data;p++, index++)
                    P->edit(p,parameters[index]);
            }
        }
    }
    return test();
}

Array1D<double> System::getAntecedents(void)
{
	Array1D<double> parameters(NumberOfAntecedents());
	
    size_t index = 0;
    Membership *P;
    for (size_t j=0;j<out;j++)
    {
        for (size_t r=0;r<N[j];r++)
        {
            for (size_t i=0;i<in;i++)
            {
                P = R[j][r].readFunction(i);
                size_t size_data = P->num_params();
                for (size_t p=0;p<size_data;p++, index++)
                    parameters[index] = P->read(p);
            }
        }
    }
    return parameters;
}

size_t System::NumberOfConsequents(void)
{
	size_t rules = 0;
    for (size_t j=0;j<out;j++)
        rules += N[j];
		
    return rules*(in+1);
}

void System::setConsequents(const double *const parameters)
{
    size_t index = 0;
    for (size_t j=0;j<out;j++)
        for (size_t r=0;r<N[j];r++)
            for (size_t i=0;i<=in;i++, index++)
                R[j][r].changeTSK(parameters[index], i);

    return;
}

Array1D<double> System::getConsequents(void)
{	
	Array1D<double> consequents(NumberOfConsequents());
	
	size_t index = 0;
    for (size_t j=0;j<out;j++)
        for (size_t r=0;r<N[j];r++)
            for (size_t i=0;i<=in;i++, index++)
                consequents[index] = R[j][r].readTSK(i);

    return consequents;
}

Array1D<double> System::evaluate(const double * const point)
{
    const size_t m = in-out; // Number of control inputs
	Array1D<double> output(out, 0.0);
	
	if (checkInputLimits(point))
		WARNINGMSG(W_InputsLimitsExceeded)
	
    int N_max = 0; // Indices
    for (size_t i=0;i<out;i++)
    {
        if (N[i] > N_max)
            N_max = N[i];
    }

    double *X = new double [out+1];
    X[0] = 1.0;
    for (size_t i=0;i<out;i++)
        *(X+i+1) = *(point+i); // State vector extension

    // Auxiliary values
    double temp, Sum_W, at, bt;
    double *W = new double [N_max]; // Matching degree of the plant

    for (size_t i=0;i<out;i++)
    {

        Sum_W = 0;
        for (size_t l=0;l<N[i];l++)
        {
            W[l] = R[i][l].activation(point); // W[l][i]
            Sum_W += W[l];
        }

        for (size_t k=0;k<=out;k++)
        {
            at = 0.0;
            for (size_t l=0;l<N[i];l++)
                at += ( W[l] * R[i][l].readTSK(k) );

//          at=at/Sum_W; // at[k][i] // The division is made at the end
            output[i] = output[i] + at*X[k];
        }
	
        for (size_t j=0;j<m;j++)
        {
            bt = 0.0;
            for (size_t l=0;l<N[i];l++)
                bt += ( W[l] * R[i][l].readTSK(j+out+1) );
				
//          bt=bt/Sum_W; // bt[j][i] // The division is made at the end
            output[i] = output[i] + bt * point[j+out]; // point is actually u here
        }

        output[i] = output[i] / (Sum_W + M_EPS);
	}

//	if (checkOutputLimits(&output[0]))
//		WARNINGMSG(W_OutputsLimitsExceeded)

    // Free dynamic memory
    delete [] W;
    delete [] X;

    return output;
}

System & System::operator=(const System &S)
{
    size_t r,j;
    if (this != &S)
    {
        if (N != NULL)
            delete [] N;
        N = new size_t[S.outputs()];
        if (R != NULL)
        {
            for (j=0;j<out;j++)
                if (R[j] != NULL)
                    delete [] R[j];
            delete [] R;
        }
		if (low_in != NULL)
			delete [] low_in;
		if (low_out != NULL)
			delete [] low_out;
		if (high_in != NULL)
			delete [] high_in;
		if (high_out != NULL)
			delete [] high_out;

        in = S.inputs();
        out = S.outputs();
        low_in = new double [in];
        high_in = new double [in];
        low_out = new double [out];
        high_out = new double [out];
        R = new Rule*[out];
		for (j=0;j<in;j++)
		{
			low_in[j] = S.in_low(j);
			high_in[j] = S.in_high(j);
		}
        for (j=0;j<out;j++)
        {
			low_out[j] = S.out_low(j);
			high_out[j] = S.out_high(j);
            N[j] = S.rules(j);
            R[j] = new Rule[N[j]];
            for (r=0;r<N[j];r++)
                R[j][r] = S.R[j][r];
        }
    }
    return *this;
}

System::System()
{
    in = 0;
    out = 0;
    N = NULL;
    R = NULL;
	low_in = NULL;
	high_in = NULL;
	low_out = NULL;
	high_out = NULL;
}

System::System(const System &S)
{
    size_t r,i,j;
    in = S.inputs();
    out = S.outputs();
    low_in = new double [in];
	high_in = new double [in];
	low_out = new double [out];
	high_out = new double [out];
	N = new size_t[out];
    R = new Rule*[out];
	for (i=0;i<in;i++)
	{
        low_in[i] = S.in_low(i);
        high_in[i] = S.in_high(i);
	}
    for (j=0;j<out;j++)
    {
		low_out[j] = S.out_low(j);
		high_out[j] = S.out_high(j);
        N[j] = S.rules(j);
        R[j] = new Rule[N[j]];
        for (r=0;r<N[j];r++)
            R[j][r] = S.R[j][r];
    }
}

System::System(size_t in, size_t out, size_t *N)
{
    size_t r,j,i;
    this->in = in;
    this->out = out;
    low_in = new double [in];
	high_in = new double [in];
	low_out = new double [out];
	high_out = new double [out];
    this->N = new size_t[out];
	for (i=0;i<in;i++)
	{
		low_in[i] = LOW_DEFAULT_LIMIT;
		high_in[i] = HIGH_DEFAULT_LIMIT;
	}
    R = new Rule*[out];
    for (j=0;j<out;j++)
    {
		low_out[j] = LOW_DEFAULT_LIMIT;
		high_out[j] = HIGH_DEFAULT_LIMIT;
        this->N[j] = N[j];
        R[j] = new Rule[N[j]];
        for (r=0;r<N[j];r++)
            R[j][r].initialize(in);
    }
}

System::~System()
{
    size_t j;
	if (low_in != NULL)
	{
		delete [] low_in;
		low_in = NULL;
	}
	if (high_in != NULL)
	{
		delete [] high_in;
		high_in = NULL;
	}
	if (low_out != NULL)
	{
		delete [] low_out;
		low_out = NULL;
	}
	if (high_out != NULL)
	{
		delete [] high_out;
		high_out = NULL;
	}
    if (N != NULL)
    {
        delete [] N;
        N = NULL;
    }
    if (R != NULL)
    {
        for (j=0;j<out;j++)
        {
            if (R[j] != NULL)
                delete [] R[j];
            R[j] = NULL;
        }
        delete [] R;
        R = NULL;
    }
}
