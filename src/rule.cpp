#include <flt/rule.hpp>

using namespace FLT;
using namespace TNT;
using namespace std;

// Implementations ############################################################

void Rule::clean(void)
{
    if (membership != NULL && in>0)
    {
        for (size_t v=0;v<in;v++)
        {
            if (membership[v] != NULL)
            {
                delete membership[v];
                membership[v] = NULL;
            }
        }
        delete [] membership;
    }
    membership = NULL;
    if (TSK != NULL)
    {
        delete [] TSK;
        TSK=NULL;
    }
}

size_t Rule::n_inputs(void) const
{
	return in;
}

int Rule::test(void) const
{
    if (membership==NULL || in==0)
        return 1;
	int error = 0;
	for (size_t v=0;v<in;v++)
	{
		if (membership[v]==NULL)
			return 1;
		int test = membership[v]->test();
		if (test > 0)
			return 1;
		if (test < 0)
			error = test;
	}
    return error;
}

int Rule::changeTSK(double value, size_t input)
{
    if (input>in)
        return 1;
    TSK[input] = value;
    return 0;
}

void Rule::changeTSK(const double* const TSK)
{
	for (int i=0; i<in; i++)
		this->TSK[i] = TSK[i];
}

double Rule::readTSK(size_t input) const
{
    if (input>in)
		ERRORMSGVAL(E_ArrayExceded, 0.0)
    return TSK[input];
}

Array1D<double> Rule::readTSK(void) const
{
	Array1D<double> consequent(in, TSK);
    return consequent;
}

Membership *Rule::readFunction(size_t input)
{
    if (input>=in)
		ERRORMSGVAL(E_ArrayExceded, NULL)
    return membership[input];
}

int Rule::changeFunction(const Membership &P, size_t input)
{
    if (input >= in)
        return 1;
    if (membership[input])
		delete membership[input]; // Ensure memory clean 
	membership[input] = createMF(P.type());
	*(membership[input]) = P; // Make the change
	return 0;
}

int Rule::changeFunction(TYPE_MF type, size_t input)
{
    if (input >= in)
        return 2;
	delete membership[input];
    membership[input] = createMF(type);
	if (membership[input] == NULL)
		return 1;
    return 0;
}

void Rule::initialize(size_t in)
{
    if (this->in != 0)
        clean();
    this->in = in;
    TSK = new double [(in+1)];
    membership = new Membership *[in];
    for (size_t v=0;v<in;v++)
        membership[v] = NULL;
}

size_t Rule::NumberOfAntecedents(void)
{
    size_t size_data = 0;

    for (size_t i=0;i<in;i++)
        size_data += membership[i]->num_params();
		
    return size_data;
}

int Rule::setAntecedents(const double *const parameters)
{
    Membership *P;
    size_t index = 0;

    for (size_t i=0;i<in;i++)
    {
        P = membership[i];
        size_t size_data = P->num_params();
        for (size_t p=0;p<size_data;p++, index++)
            P->edit(p,parameters[index]);
    }
    return test();
}

Array1D<double> Rule::getAntecedents(void)
{
	Array1D<double> parameters(NumberOfAntecedents());

    size_t index = 0;
    Membership *P;
    for (size_t i=0;i<in;i++)
    {
        P = membership[i];
        size_t size_data = P->num_params();
        for (size_t p=0;p<size_data;p++, index++)
            parameters[index] = P->read(p);
    }
    return parameters;
}

void Rule::setConsequents(const double *const parameters)
{  
    for (size_t i=0;i<=in;i++)
        TSK[i] = parameters[i];
		
    return;
}

Array1D<double> Rule::getConsequents(void)
{
	Array1D<double> parameters (in+1);

    for (size_t i=0;i<=in;i++)
        parameters[i] = TSK[i];
	
    return parameters;
}

double Rule::activation(const double * const point, double *dW) const
{
	Membership *P;
	double W = 1.0;
    for (size_t k=0;k<in;k++)
	{
		P = membership[k];
		if (P->type() != ANYMF) // ANYMF doesn't affect to the system
			W *= P->eval(point[k]);
		if (dW)
		{
			dW[k] = P->evalder(point[k]);
			for (int q=0;q<in;q++)
			{
				if (q!=k)
				{
					P = membership[q];
					if (P->type() != ANYMF) // ANYMF doesn't affect to the system
						dW[k] *= P->eval(point[q]);
				}
			}
		}
	}

	return W;
}

Rule& Rule::operator=(const Rule &R)
{
    if(this==&R)
        return *this;
    initialize(R.n_inputs());
    for (size_t v=0;v<in;v++)
    {
        membership[v]=createMF(R.membership[v]->type());
        (*membership[v])=(*R.membership[v]);
    }
    for (size_t v=0;v<=in;v++)
        TSK[v]=R.readTSK(v);
    return *this;
}

Rule::Rule()
{
    TSK=NULL;
    membership=NULL;
    in=0;
}

Rule::Rule(size_t in)
{
    this->in=0;
    initialize(in);
}

Rule::Rule(const Rule &R)
{
    in=0;
    initialize(R.n_inputs());
    for (size_t v=0;v<in;v++)
    {
        membership[v]=createMF(R.membership[v]->type());
        (*membership[v])=(*R.membership[v]);
    }
    for (size_t v=0;v<=in;v++)
        TSK[v]=R.readTSK(v);
}

Rule::~Rule(){clean();}
