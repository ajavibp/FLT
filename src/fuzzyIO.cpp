#include <flt/fuzzyIO.hpp>

using namespace FLT;
using namespace std;

// TYPE_MF
ostream &std::operator<<(ostream &F, TYPE_MF &T)
{
    F << MF_NAMES[T];
    return F;
}

istream &std::operator>>(istream &F, TYPE_MF &T)
{
	int iT;
    char type[MAX_SIZE_TYPE_NAME];

    F >> type;
    for (iT=(int)ANYMF;iT<=(int)ZMF;iT++)
        if (!strcasecmp(MF_NAMES[(TYPE_MF)iT],type))
            break;
    if (iT>(int)ZMF)
    {
    	cerr << E_BadMF << endl;
        exit(1);
    }
    T=(TYPE_MF)iT;
    return F;
}

// Membership
ostream &std::operator<<(ostream &F, Membership &P)
{
    F.precision(PRECISION);

    TYPE_MF T = P.type();
    F << T << "\t";

    for (size_t i=0;i<P.num_params();i++)
        F << P.read(i) << "\t";
    
    return F;
}

istream &std::operator>>(istream &F, Membership &P)
{
    for (size_t i=0;i<P.num_params();i++)
    {
    	double d;
        F >> d;
        if(P.edit(i,d))
        {
        	cerr << E_ChangingParameter << endl;
            exit(1);
        }
    }

    return F;
}

// Rule
ostream &std::operator<<(ostream &F, Rule &R)
{
    F.precision(PRECISION);

    for (size_t i=0; i<R.n_inputs(); i++)
        F << *R.readFunction(i);
    F << "\n";

    for (size_t i=0; i<=R.n_inputs(); i++)
        F << R.readTSK(i) << "\t";
    F << "\n";

    return F;
}

istream &std::operator>>(istream &F, Rule &R)
{
    double temp;
    TYPE_MF t;

    for (size_t i=0; i<R.n_inputs(); i++)
    {
        F >> t;
        if (R.changeFunction(t, i))
            exit(1);

        Membership *P = createMF(t);
        F >> *P;

        if(R.changeFunction(*P, i))
        {
        	cerr << E_WritingRule << endl;
            exit(1);
        }

        delete P;
    }

    for (size_t i=0; i<=R.n_inputs(); i++)
    {
        F >> temp;
        if(R.changeTSK(temp, i))
        {
        	cerr << E_WritingRule << endl;
            exit(1);
        }
    }

    return F;
}


// System
ostream &std::operator<<(ostream &F, System &S)
{
    size_t *rules = new size_t[S.outputs()];
    const size_t m = S.outputs();
	const size_t in = S.inputs();

    for (size_t j=0; j<m; j++)
        *(rules + j) = S.rules(j);

    F.precision(PRECISION);

    F << in << "\t" << m << "\n";
    for (size_t j=0;j<m; j++)
        F << S.rules(j) << "\t";
    F << "\n\n";

	Rule *R;
    for (size_t j=0; j<m; j++)
    {
        for (size_t r=0;r<rules[j];r++)
        {
            R = S.readRule(j,r);
            F << *R << "\n";
        }
    }
    delete [] rules;

	for (size_t i=0;i<in;i++)
		F << S.in_low(i) << "\t" << S.in_high(i) << "\n";
	F<<"\n";

	for (size_t j=0;j<m;j++)
		F << S.out_low(j) << "\t" << S.out_high(j) << "\n";

    return F;
}

istream &std::operator>>(istream &F, System &S)
{
    size_t in,out;
    F >> in;
    F >> out;

    size_t *rules = new size_t[out];
    for (size_t i=0;i<out;i++)
        F >> rules[i];

    S.initialize(in, out, rules);

    for (size_t i=0;i<out;i++)
    {
        for (size_t r=0;r<rules[i];r++)
        {
        	Rule R(in);
            F >> R;

            if(S.changeRule(R,r,i))
            {
                cerr << E_RulesExceded << endl;
                exit(1);
            }
        }
    }
    delete [] rules;

	double *limits_min = new double[in];
	double *limits_max = new double[in];

	for (size_t i=0;i<in;i++)
	{
		F >> limits_min[i];
		F >> limits_max[i];
	}

	S.in_low(limits_min);
	S.in_high(limits_max);

	delete [] limits_min;
	delete [] limits_max;

	limits_min = new double[out];
	limits_max = new double[out];

	for (size_t i=0;i<out;i++)
	{
		F >> limits_min[i];
		F >> limits_max[i];
	}

	S.out_low(limits_min);
	S.out_high(limits_max);

	delete [] limits_min;
	delete [] limits_max;

    return F;
}
