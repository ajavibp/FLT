#include <flt/derivatives.hpp>

using namespace FLT;
using namespace TNT;

Array2D<double> FLT::jacobian(System &S, const double *const x, const double *const u, TNT::Array2D<double> &B, TNT::Array1D<double> &F)
{
    const size_t n = S.outputs(); // Order
    const size_t m = S.inputs() - n; // Number of control inputs
	
	Array2D<double> A(n,n);
	if (B.dim1()!=n || B.dim2()!=m)
	{
		Array2D<double> Bnew(n,m,0.0);
		B = Bnew.copy();
	}
	if (F.dim()!=n)
	{
		Array1D<double> Fnew(n);
		F = Fnew.copy();
	}
	
    size_t *M = new size_t[n]; // Number of rules of the plant for each output
	size_t M_max = 0;
    for (size_t i=0;i<n;i++)
    {
        M[i] = S.rules(i);
        if (M[i]>M_max)
            M_max = M[i];
    }

    double temp1,temp2;
    double a1,a2,b1,b2,at,bt;
    double *X = new double [n+1]; // State vector extension
    X[0] = 1.0;
    for (size_t i=0;i<n;i++)
        X[i+1] = x[i];
    double **W = new double *[M_max]; // Matching degree of the rules of the plant
    double ***der_wx = new double**[M_max]; // dW/dx
    double ***der_wu = new double**[M_max]; // dW/du
    for (size_t l=0;l<M_max;l++)
    {
        W[l] = new double[n];
        der_wx[l] = new double*[n];
        der_wu[l] = new double*[n];
        for (size_t i=0;i<n;i++)
        {
            der_wx[l][i] = new double[n];
            der_wu[l][i] = new double[m];
        }
    }
    double *Sum_W = new double[n];
    Membership *P;
    Rule *R;
	
    // Algorithm
    for (size_t i=0;i<n;i++)
    {
        Sum_W[i] = 0.0;
        for (size_t l=0;l<M[i];l++)
        {
            R = S.readRule(i,l);
            W[l][i] = 1.0;
            for (size_t k=0;k<n;k++)
            {
                P = R->readFunction(k);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(x[k]);
            }
            for (size_t j=0;j<m;j++)
            {
                P = R->readFunction(j+n);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(u[j]);
            }
            // W[l][i]
            Sum_W[i] += W[l][i];
			
            for (size_t q=0;q<n;q++)
            {
				P = R->readFunction(q);
				if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
				{
					der_wx[l][i][q] = P->evalder(x[q]);
					for (size_t k=0;k<n;k++)
					{
						if (k!=q)
						{
							P = R->readFunction(k);
							if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
								der_wx[l][i][q] *= P->eval(x[k]);
						}
					}
					for (size_t j=0;j<m;j++)
					{
						P = R->readFunction(j+n);
						if (P->type() != ANYMF) // ANYMF doesn't affect to the system
							der_wx[l][i][q] *= P->eval(u[j]);
					}
				}
				else
					der_wx[l][i][q] = 0.0;
            } // der_wx[l][i][q]
			
            for (size_t t=0;t<m;t++)
            {
				P = R->readFunction(t+n);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
				{
					der_wu[l][i][t] = P->evalder(u[t]);
					for (size_t k=0;k<n;k++)
					{
						P = R->readFunction(k);
						if (P->type() != ANYMF) // ANYMF doesn't affect to the system
							der_wu[l][i][t] *= P->eval(x[k]);
					}
					for (size_t j=0;j<m;j++)
					{
						if (j!=t)
						{
							P = R->readFunction(j+n);
							if (P->type() != ANYMF) // ANYMF doesn't affect to the system
								der_wu[l][i][t] *= P->eval(u[j]);
						}
					}
				}
				else
					der_wu[l][i][t] = 0.0;
            }
        }// Sum_W[i]
    }
	
    for (size_t i=0;i<n;i++)
    {
        for (size_t t=0;t<m;t++) // To obtain B
        {
            bt = 0.0;
            for (size_t l=0;l<M[i];l++)
            {
                R = S.readRule(i,l);
                bt += (W[l][i]*R->readTSK(t+n+1));
            }
            bt /= (Sum_W[i]+M_EPS); // bt[t][i]
			
            temp2 = 0.0;
            for (size_t k=0;k<=n;k++)
            {
                temp1 = 0.0;
                for (size_t l=0;l<M[i];l++)
                {
                    a1 = S.readRule(i,l)->readTSK(k);
                    for (size_t p=0;p<M[i];p++)
                    {
                        a2 = S.readRule(i,p)->readTSK(k);
                        temp1 += (der_wu[l][i][t]*W[p][i]*(a1-a2));
                    }
                }
                temp1 *= (X[k]/(Sum_W[i] * Sum_W[i] + M_EPS));
                temp2 += temp1;
            }
            for (size_t j=0;j<m;j++)
            {
                temp1 = 0.0;
                for (size_t l=0;l<M[i];l++)
                {
                    b1 = S.readRule(i,l)->readTSK(j+n+1);
                    for (size_t p=0;p<M[i];p++)
                    {
                        b2 = S.readRule(i,p)->readTSK(j+n+1);
                        temp1 += (der_wu[l][i][t]*W[p][i]*(b1-b2));
                    }
                }
                temp1 *= (u[j] / (Sum_W[i]*Sum_W[i]+M_EPS));
                temp2 += temp1;
            }
            B[i][t] = bt + temp2;
        }
        for (size_t q=0;q<n;q++) // To obtain A
        {
            at = 0.0;
            for (size_t l=0;l<M[i];l++)
            {
                R = S.readRule(i,l);
                at += (W[l][i]*R->readTSK(q+1));
            }
            at /= (Sum_W[i]+M_EPS); // at[q][i]
			
            temp2 = 0.0;
            for (size_t k=0;k<=n;k++)
            {
                temp1 = 0.0;
                for (size_t l=0;l<M[i];l++)
                {
                    a1 = S.readRule(i,l)->readTSK(k);
                    for (size_t p=0;p<M[i];p++)
                    {
                        a2 = S.readRule(i,p)->readTSK(k);
                        temp1 += (der_wx[l][i][q]*W[p][i]*(a1-a2));
                    }
                }
                temp1 *= (X[k]/(Sum_W[i]*Sum_W[i]+M_EPS));
                temp2 += temp1;
            }
            for (size_t j=0;j<m;j++)
            {
                temp1 = 0.0;
                for (size_t l=0;l<M[i];l++)
                {
                    b1 = S.readRule(i,l)->readTSK(j+n+1);
                    for (size_t p=0;p<M[i];p++)
                    {
                        b2 = S.readRule(i,p)->readTSK(j+n+1);
                        temp1 += (der_wx[l][i][q]*W[p][i]*(b1-b2));
                    }
                }
                temp1 *= (u[j]/(Sum_W[i]*Sum_W[i]+M_EPS));
                temp2 += temp1;
            }
            A[i][q] = at + temp2;
        }
    }
	
    // To obtain f(x,u)
	Array1D<double> xu(n+m);
	for (size_t i=0;i<n;i++)
		xu[i] = x[i];
	for (size_t j=0;j<m;j++)
		xu[j+n] = u[j];
	F = S.evaluate(&xu[0]);

    // Free the dynamic memory
    delete [] M;
    for (size_t l=0;l<M_max;l++)
    {
        delete [] W[l];
        for (size_t i=0;i<n;i++)
        {
            delete [] der_wu[l][i];
            delete [] der_wx[l][i];
        }
        delete [] der_wu[l];
        delete [] der_wx[l];
    }
    delete [] W;
    delete [] der_wx;
    delete [] der_wu;
    delete [] Sum_W;
    delete [] X;

    return A;
}

TNT::Array2D<double> FLT::jacobian(System &S, System &C, const double *const x)
{
	/* ### To do ###########################################################
	 * Change this function to make use of activation and other new methods.
	 * ##################################################################### */
	
	// Indices
	size_t i, k, q, j, r, l, s, p, v; 
	
    // System and Controller constants
    const size_t n = S.outputs(); // System order
    if (n != C.inputs() || n != S.inputs()-C.outputs())
        ERRORMSGVAL(E_NoCoherent, Array2D<double>())
		
    const size_t m = C.outputs(); // Number of control inputs
	
    size_t *M = new size_t[n];
    size_t M_max = 0; // The maximum number of rules of the plant
    for (i=0;i<n;i++)
    {
        M[i] = S.rules(i); // Number of rules of the plant in the ith output
        if (M[i] > M_max)
            M_max = M[i];
    }
    size_t *N = new size_t[m];
    size_t N_max = 0; // The maximum number of rules of the controller
    for (j=0;j<m;j++)
    {
        N[j] = C.rules(j); // Number of rules of the plant in the jth output
        if (N[j] > N_max)
            N_max = N[j];
    }

    // Auxiliary values
    // using M_max and N_max to make the memory reservation
    Membership *P;
    Rule *R;

    double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, temp4 = 0.0;
    double Sum2 = 0.0, Sum3 = 0.0;
    double a1, a2, b1, b2, c1, c2;
    double at;  // Variable coefficient 'a' for the plant. 

    double *X = new double [n+1]; // Extension of the state vector
    X[0] = 1.0;
    for (i=0;i<n;i++)
        X[i+1] = x[i];
	
    double *Sum_activ = new double[m];
    double *Sum_W = new double[n];
    double **activ = new double *[N_max]; // Matching degree of the rules of the controller (omega)
    double **ct = new double *[n+1]; // Variable coefficient for the controller
    double *u = new double [m]; // Control signal
    double **W = new double *[M_max]; // Matching degree of the rules of the plant
    double *bt = new double [m]; // Variable coefficient 'b' for the plant 
    double **der_u = new double *[m]; // du/dx
    double ***der_w = new double **[M_max]; // dw/dx
    double ***der_activ = new double **[N_max]; //domega/dx

    for (j=0;j<m;j++)
        der_u[j] = new double [n];
    for (l=0;l<M_max;l++)
    {
        W[l] = new double[n];
        der_w[l] = new double *[n];
        for (i=0;i<n;i++)
            der_w[l][i] = new double[n];
    }
    for (r=0;r<N_max;r++)
    {
        activ[r] = new double[m];
        der_activ[r] = new double *[n];
        for (j=0;j<m;j++)
            der_activ[r][j] = new double[n];
    }
    for (i=0;i<=n;i++)
        ct[i] = new double[m];
		
	Array2D<double> J(n,n);

    // Calculation procedure
    for (j=0;j<m;j++)
    {
        Sum_activ[j] = 0.0;
        for (r=0;r<N[j];r++)
        {
            R = C.readRule(j,r);
            activ[r][j] = 1.0;
            for (k=0;k<n;k++)
            {
                P = R->readFunction(k);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
					activ[r][j] *= P->eval(x[k]);
            }
            Sum_activ[j] += activ[r][j];
        } // activ[r][j] is calculated here
        u[j] = 0.0;
        for (k=0;k<=n;k++)
        {
            ct[k][j] = 0.0;
            for (r=0;r<N[j];r++)
                ct[k][j] += (activ[r][j] * C.readRule(j,r)->readTSK(k));
            ct[k][j] /= Sum_activ[j]; // ct[k][j] is calculated here
            u[j] += (X[k] * ct[k][j]);
        } // u[j] is calculated here
    }
    for (i=0;i<n;i++)
    {
        Sum_W[i] = 0.0;
        for (l=0;l<M[i];l++)
        {
            R = S.readRule(i,l);
            W[l][i] = 1.0;
            for (k=0;k<n;k++)
            {
                P = R->readFunction(k);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(x[k]);
            }
            for (j=0;j<m;j++)
            {
                P = R->readFunction(j+n);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(u[j]);
            }
            // W[l][i] is calculated here
            Sum_W[i] += W[l][i];
        }
    }
    for (j=0;j<m;j++)
    {
        for (q=0;q<n;q++)
        {
            for (r=0;r<N[j];r++)
            {
                R = C.readRule(j,r);
				P = R->readFunction(q);
				if (P->type() != ANYMF) // ANYMF doesn't affect to the system
				{
					der_activ[r][j][q] = P->evalder(x[q]);
					temp1 = 1.0;
					for (k=0;k<n;k++)
					{
						if (k != q)
						{
							P = R->readFunction(k);
							if (P->type() != ANYMF) // ANYMF doesn't affect to the system
								temp1 *= P->eval(x[k]);
						}
					}
					der_activ[r][j][q] *= temp1;
				}
				else
					der_activ[r][j][q] = 0.0;
            }
            // der_activ[r][j][q] is calculated here
            temp1 = 0;
            for (k=0;k<=n;k++)
            {
                temp2 = 0;
                for (r=0;r<N[j];r++)
                {
                    c1 = C.readRule(j,r)->readTSK(k);
                    for (s=0;s<N[j];s++)
                    {
                        c2 = C.readRule(j,s)->readTSK(k);
                        temp2 += (der_activ[r][j][q] * activ[s][j] * (c1-c2));
                    }
                }
                temp2 *= X[k];
                temp1 += temp2;
            }
            der_u[j][q] = ct[q+1][j] + temp1 / pow(Sum_activ[j],2);
        } // der_u[j][q] is calculated here
    }
    for (i=0;i<n;i++)
    {
        for (j=0;j<m;j++)
        {
            bt[j] = 0.0;
            for (l=0;l<M[i];l++)
                bt[j] += (W[l][i] * S.readRule(i,l)->readTSK(j+n+1));
            bt[j] /= Sum_W[i]; // bt[j][i] is calculated here
        }
        for (q=0;q<n;q++)
        {
            at = 0.0;
            Sum2 = 0.0;
            for (j=0;j<m;j++)
                Sum2 += (bt[j] * ct[q+1][j]); // Sum2[i][q] is calculated here
            for (l=0;l<M[i];l++)
            {
                R = S.readRule(i,l);
                at += (W[l][i] * R->readTSK(q+1));
				P = R->readFunction(q);
				if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
				{
					temp1 = P->evalder(x[q]);
					for(j=0;j<m;j++)
					{
						P = R->readFunction(j+n);
						if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
							temp1 *= P->eval(u[j]);
					}
					temp3 = 1.0;
					for(k=0;k<n;k++)
					{
						if (k!=q)
						{
							P = R->readFunction(k);
							if (P->type() != ANYMF) // ANYMF doesn't affect to the system
								temp3 *= P->eval(x[k]);
						}
					}
					temp1 *= temp3;
					temp3 *= R->readFunction(q)->eval(x[q]); // R->readFunction(q) != ANYMF (see if)
					temp2 = 0.0;
					for (v=0;v<m;v++)
					{
						P = R->readFunction(v+n);
						if (P->type() != ANYMF) // ANYMF doesn't affect to the system
						{
							temp4 = der_u[v][q] * P->evalder(u[v]);
							for (j=0;j<m;j++)
							{
								if (j != v)
								{
									P = R->readFunction(j+n);
									if (P->type() != ANYMF) // ANYMF doesn't affect to the system
										temp4 *= P->eval(u[j]);
								}
							}
							temp2 += temp4;
						} // else temp4 = 0 -> temp2 don't change						
					}
					der_w[l][i][q] = temp1 + temp2 * temp3;
				}
				else
					der_w[l][i][q] = 0.0;
				// der_w[l][i][q] is calculated here
            }
            at /= Sum_W[i]; // at[q][i] is calculated here
            Sum3 = 0.0;
            for (k=0;k<=n;k++)
            {
                temp1 = 0.0;
                for (l=0;l<M[i];l++)
                {
                    a1 = S.readRule(i,l)->readTSK(k);
                    for (p=0;p<M[i];p++)
                    {
                        a2 = S.readRule(i,p)->readTSK(k);
                        temp1 += (der_w[l][i][q] * W[p][i] * (a1-a2));
                    }
                }
                temp1 /= pow(Sum_W[i],2.0);
                temp2 = 0.0;
                for (j=0;j<m;j++)
                {
                    temp3 = 0.0;
                    for (l=0;l<M[i];l++)
                    {
                        b1 = S.readRule(i,l)->readTSK(j+n+1);
                        for (p=0;p<M[i];p++)
                        {
                            b2 = S.readRule(i,p)->readTSK(j+n+1);
                            temp3 += (der_w[l][i][q] * W[p][i] * (b1-b2));
                        }
                    }
                    temp3 *= ct[k][j]/pow(Sum_W[i],2.0);
                    temp4 = 0.0;
                    for (r=0;r<N[j];r++)
                    {
                        c1 = C.readRule(j,r)->readTSK(k);
                        for (s=0;s<N[j];s++)
                        {
                            c2 = C.readRule(j,s)->readTSK(k);
                            temp4 += (der_activ[r][j][q] * activ[s][j] * (c1-c2));
                        }
                    }
                    temp4 *= bt[j] / pow(Sum_activ[j],2.0);
                    temp2 += temp3 + temp4;
                }
                Sum3 += ((temp1 + temp2) * X[k]);
            }
            J[i][q] = at + Sum2 + Sum3; // J[i][q] is calculated here
        }
    }

    // Free dynamic memory
    for (j=0;j<m;j++)
        delete [] der_u[j];
    delete [] der_u;
    der_u = NULL;
    delete [] bt;
    bt = NULL;
    for (i=0;i<=n;i++)
        delete [] ct[i];
    delete [] ct;
    ct = NULL;
    for (l=0;l<M_max;l++)
    {
        delete [] W[l];
        for (i=0;i<n;i++)
            delete [] der_w[l][i];
        delete [] der_w[l];
    }
    delete [] W;
    W = NULL;
    delete [] der_w;
    der_w = NULL;
    for (r=0;r<N_max;r++)
    {
        delete [] activ[r];
        for (j=0;j<m;j++)
            delete [] der_activ[r][j];
        delete [] der_activ[r];
    }
    delete [] activ;
    activ = NULL;
    delete [] der_activ;
    der_activ = NULL;
    delete [] Sum_W;
    Sum_W = NULL;
    delete [] Sum_activ;
    Sum_activ = NULL;
    delete [] u;
    u = NULL;
    delete [] X;
    X = NULL;
    delete [] N;
    N = NULL;
    delete [] M;
    M = NULL;
	
    return J;
}

TNT::Array2D<double> FLT::jacobianAprox(System &S, const double *const x, const double *const u, TNT::Array2D<double> &B, TNT::Array1D<double> &F, double h)
{
	double *xu, *xu2;

    const size_t n = S.outputs();
	const size_t m = S.inputs() - n;
	
	Array2D<double> A(n,n);
	if (B.dim1()!=m || B.dim2()!=n)
	{
		Array2D<double> Bnew(m,n);
		B = Bnew.copy();
	}

	xu = new double[n+m];
	xu2 = new double[n+m];

	for (size_t i=0;i<n;i++)
	{
		xu[i] = x[i];
		xu2[i] = x[i];
	}
	for (size_t i=0;i<m;i++)
	{
		xu[n+i] = u[i];
		xu2[n+i] = u[i];
	}

	Array1D<double> dx = S.evaluate(xu);
	Array1D<double> dx2;

	for (size_t j=0;j<n;j++)
	{
		xu2[j] = xu2[j]+h;

		dx2 = S.evaluate(xu2);
		for (size_t i=0;i<n;i++)
			A[i][j] = (dx2[i] - dx[i]) / h;

		xu2[j] = xu[j];
	}
	for (size_t j=0;j<m;j++)
	{
	    size_t index = n + j;
		xu2[index] = xu2[index] + h;
		dx2 = S.evaluate(xu2);
		for (size_t i=0;i<n;i++)
			B[i][j] = (dx2[i] - dx[i]) / h;
		xu2[index] = xu[index];
	}
	// To obtain f(x0,u0)
	F = S.evaluate(xu);
	
	delete []xu2;
	delete []xu;
	return A;
}

TNT::Array2D<double> FLT::jacobianAprox(System &S, System &C, const double * const x, double h)
{
	size_t i,j,k;
	double *x2;

    // System and controller constants
    // n = System order
    const size_t n = S.outputs();

    if (n != C.inputs() || n != (S.inputs() - C.outputs()))
        ERRORMSGVAL(E_NoCoherent, Array2D<double>())
	
    const size_t m = C.outputs(); // Number of control inputs

	Array2D<double> J(n,n);
	Array1D<double> dx = evaluate(x,S,C);
	if (!dx.dim())
		ERRORMSGVAL(E_NoCoherent, Array2D<double>())
	
	x2 = new double[n];
	Array1D<double> dx2(n);

	for (j=0;j<n;j++)
	{
		for (k=0;k<n;k++)
			x2[k] = x[k];
		x2[j] = x2[j]+h;

		dx2 = evaluate(x2,S,C);
		if (!dx2.dim())
		{
			delete []x2;
			ERRORMSGVAL(E_NoCoherent, Array2D<double>())
		}
		for (i=0;i<n;i++)
			J[i][j] = (dx2[i] - dx[i])/h; // Aproximation of the jacobian matrix
	}
	
	// Free dynamic memory
	delete []x2;
	return J;
}


double FLT::derconseq(System &S, double *x, size_t output, size_t rule, size_t parameter, int &error)
{
	if (S.outputs()<=output || S.rules(output)<=rule || S.inputs()<parameter)
	{
		error = 1;
		return 0.0;
	}
	double sumW = 0.0;
	double W_LI;
	
	for (size_t l=0;l<S.rules(output);l++)
	{
		double w = S.readRule(output, l)->activation(x);
		sumW += w;
		if (l==rule)
			W_LI = w;
	}

	// State vector extension X=[1;x]
	if (parameter!=0)
		return W_LI*x[parameter-1]/sumW;
	else
		return W_LI/sumW;
}


TNT::Array2D<double> FLT::derconseq(System &S, double *x)
{
	const size_t n = S.inputs();
	const size_t m = S.outputs();
	int error = 0;
	const size_t num_conseq4rule = n + 1;
	Array2D<double> dS_dconseq(m,S.NumberOfConsequents(),0.0);
	
	for (size_t i=0;i<m;i++) // Output for the consequent
    {
		size_t n_data=0;
		for (size_t I=0;I<m;I++) // Outputs of the Model
		{
			if (i!=I)
			{
				// dS_dconseq[I][n_data] = 0.0 because dS_dconseq has been initilized to 0.0
				n_data += (S.rules(I)*num_conseq4rule);
				continue;
			}
			
			// i==I
			size_t N = S.rules(i);
			for (size_t L=0;L<N;L++) // Rules of the Model
			{
				for (size_t J=0;J<=n;J++) //Inputs of the Model
				{
					dS_dconseq[I][n_data] = derconseq(S, x, I, L, J, error);
					if (error)
					{
						Array2D<double> empty;
						return empty;
					}
					n_data++;
				}
			}
		}
	}
	return dS_dconseq;
}


TNT::Array2D<double> FLT::deraffine(System &S, double *x)
{
	const size_t n = S.inputs();
	const size_t m = S.outputs();
	int error = 0;
	size_t numberofRules = 0;
	for (size_t j=0;j<m;j++)
		numberofRules += S.rules(j);
	Array2D<double> dS_daffine(m,numberofRules,0.0);
	
	for (size_t i=0;i<m;i++) // Output for the consequent
    {
		size_t n_data=0;
		for (size_t I=0;I<m;I++) // Outputs of the Model
		{
			if (i!=I)
			{
				// dS_dconseq[I][n_data] = 0.0 because dS_dconseq has been initilized to 0.0
				n_data += S.rules(I);
				continue;
			}
			
			// i==I
			size_t N = S.rules(i);
			for (size_t L=0;L<N;L++) // Rules of the Model
			{
				dS_daffine[I][n_data] = derconseq(S, x, I, L, (size_t)0, error);
				if (error)
				{
					Array2D<double> empty;
					return empty;
				}
				n_data++;
			}
		}
	}
	return dS_daffine;
}

double FLT::derantec(System &S, double *x, size_t input, size_t output, size_t rule, size_t parameter, int &error)
{
	const size_t n = S.inputs();
	const size_t m = S.outputs();
	const size_t N = S.rules(output);
	if (m<=output || N<=rule || n<=input)
	{
		error = 1;
		return 0.0;
	}
	
	Rule *R = S.readRule(output, rule);
	if (!R)
	{
		error = 1;
		return 0.0;
	}
	
	Membership *M;
	double dW_dantec = 1.0;
	for (size_t q=0;q<n;q++)
	{
		if (q == input)
			continue;
		M = R->readFunction(q);
		if (!M)
		{
			error = 1;
			return 0.0;
		}
		if (M->type() == ANYMF) // ANYMF doesn't affect to the system
			continue;
		dW_dantec *= M->eval(x[q]);
	}
	
	M = R->readFunction(input);
	if (M->type() == ANYMF) // ANYMF doesn't affect to the system
		return 0.0;
	dW_dantec *= M->paramder(parameter, x[input]);

	double *X = new double [n+1]; // Extended state vector
    X[0]=1.0;
    for (size_t i=0;i<n;i++)
        X[i+1] = x[i];
	
	double sum1=0.0;
	for (size_t j=0;j<=n;j++)
	{
		double sum2 = 0.0;
		double sumW = 0.0;
		for (size_t l=0;l<N;l++)
		{
			double W_lI = S.readRule(output, l)->activation(x);
			sumW += W_lI;
			sum2 += ( W_lI*( S.readRule(output, rule)->readTSK(j) - S.readRule(output, l)->readTSK(j) ) );
		}
		sum1 += ( sum2*X[j]/pow(sumW,2) );
	}
	delete []X;
	
	return dW_dantec*sum1;
}


TNT::Array2D<double> FLT::derantec(System &S, double *x)
{	
	const size_t n = S.inputs();
	const size_t m = S.outputs();
	Array2D<double> dS_dantec(m,S.NumberOfAntecedents(),0.0);

	for (size_t i=0;i<m;i++) // Output for the antecedent
    {
		size_t n_data=0;
		for (size_t I=0;I<m;I++) // Outputs of the Model
		{
			if (i!=I)
			{
				// dS_dantec[I][n_data] = 0.0 because dS_dantec has been initilized to 0.0
				for (size_t l=0;l<S.rules(I);l++)
					n_data += S.readRule(I,l)->NumberOfAntecedents();
			}
			else
			{
				for (size_t L=0;L<S.rules(I);L++) // Rules of the Model
				{
					Rule *R = S.readRule(I, L);
					if (!R)
					{
						Array2D<double> empty;
						return empty;
					}
		
					for (size_t J=0;J<n;J++) // Inputs of the Model
					{
						Membership *M = R->readFunction(J);
						if (!M)
						{
							Array2D<double> empty;
							return empty;
						}					
						size_t num_params = M->num_params();
						for (size_t p=0;p<num_params;p++) // Parameters of the antecedent
						{
							int error = 0;
							dS_dantec[i][n_data] = derantec(S, x, J, I, L, p, error);
							if (error)
							{
								Array2D<double> empty;
								return empty;
							}
							n_data++;
						}
					}
				}
			}
		}
	}
	
	return dS_dantec;
}

TNT::Array2D<double> FLT::derfuzzy(System &S, double *x)
{	
	const size_t n = S.inputs();
	const size_t m = S.outputs();
	const size_t num_conseq4rule = n + 1;
	
	Array2D<double> dS_dparam(m, S.NumberOfAntecedents() + S.NumberOfConsequents(), 0.0);

	for (size_t i=0;i<m;i++) // Output for the antecedent
    {
		size_t n_data = 0;
		for (size_t I=0;I<m;I++) // Outputs of the Model
		{
			if (i!=I)
			{
				// dS_dparam[I][n_data] = 0.0 because dS_dantec has been initilized to 0.0
				for (size_t l=0;l<S.rules(I);l++)
					n_data += S.readRule(I,l)->NumberOfAntecedents();
				n_data += (S.rules(I) * num_conseq4rule);
			}
			else
			{
				// Extract data from antecedents
				size_t N = S.rules(I);
				for (size_t L=0;L<N;L++) // Rules of the Model
				{
					Rule *R = S.readRule(I, L);
					if (!R)
					{
						Array2D<double> empty;
						return empty;
					}
		
					for (size_t J=0;J<n;J++) // Inputs of the Model
					{
						Membership *M = R->readFunction(J);
						if (!M)
						{
							Array2D<double> empty;
							return empty;
						}					
						size_t num_params = M->num_params();
						for (size_t p=0;p<num_params;p++) // Parameters of the antecedent
						{
							int error = 0;
							dS_dparam[i][n_data] = derantec(S, x, J, I, L, p, error);
							if (error)
							{
								Array2D<double> empty;
								return empty;
							}
							n_data++;
						}
					}
				}
				
				// Extract data from conequents
				for (size_t L=0;L<N;L++) // Rules of the Model
				{
					for (size_t J=0;J<=n;J++) //Inputs of the Model
					{
						int error=0;
						dS_dparam[I][n_data] = derconseq(S, x, I, L, J, error);
						if (error)
						{
							Array2D<double> empty;
							return empty;
						}
						n_data++;
					}
				}
			}
		}
	}
	
	return dS_dparam;
}

