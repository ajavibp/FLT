#include <flt/utilities.hpp>

#include <fstream>
#include <flt/fuzzyIO.hpp>

using namespace FLT;
using namespace TNT;
using namespace std;

int FLT::System2TXT(System &S, const char *file)
{
    ofstream fs(file);
    if(!fs.good() || S.test()>0)
    	return 1;
    
    fs << S;
    
    if(fs.good())
    {
    	fs.close();
    	return 0;
    }
    else
    {
    	fs.close();
    	return 1;
    }
}

System FLT::TXT2System(const char *file)
{
    ifstream fs(file);
	System emptyS;
	
    if(!fs.good())
    	return emptyS;
    
    System S;
	fs >> S;    

	fs.close();

    if (S.test()>0)
        return emptyS;
	else
		return S;
}

int FLT::printSystem(const char *file,System &S, char *inputs[], char *outputs[], int accuracy)
{
    size_t l,i,j,p;
    double value;
    FILE *Fs=fopen(file,"wt");
    if (!Fs)
        return 1;
    Rule *R;
    Membership *P;
    char type[MAX_CHAR_LONG], temp[MAX_CHAR_LONG];
    char **InName,**OutName; // Names
	
    size_t in = S.inputs();
    size_t out = S.outputs();
    size_t *N = new size_t[out];
    for (j=0;j<out;j++)
        N[j] = S.rules(j);
    InName = new char *[in];
    OutName = new char *[out];

    for (i=0;i<in;i++)
    {
        InName[i] = new char [MAX_CHAR_LONG];
        if (inputs==NULL)
            sprintf(InName[i],"%s%lu",U_FISInput,i+1);
        else sprintf(InName[i],"%s",inputs[i]);
    }
    for (i=0;i<out;i++)
    {
        OutName[i] = new char [MAX_CHAR_LONG];
        if (outputs==NULL)
            sprintf(OutName[i],"%s%lu",U_FISOutput,i+1);
        else sprintf(OutName[i],"%s",outputs[i]);
    }
    for (j=0;j<out;j++)
    {
        for (l=0;l<N[j];l++)
        {
			R = S.readRule(j,l);
            if (!R)
            {
                fclose(Fs);
                printf("%s\n",E_RulesExceded);
                for (i=0;i<in;i++)
                    delete [] InName[i];
                delete [] InName;
                for (i=0;i<out;i++)
                    delete [] OutName[i];
                delete [] OutName;
                delete [] N;
                return 2;
            }
			
            fprintf(Fs,"%s ",U_If);
            for (i=0;i<in;i++)
            {
                P = R->readFunction(i);
                fprintf(Fs,"%s %s %s(",InName[i],U_Is,MF_NAMES[(int)P->type()]);
                sprintf(temp,"%%.%dg; ",accuracy);
                for (p=0;p<P->num_params();p++)
                    fprintf(Fs,temp,P->read(p));
                fseek(Fs,-2,SEEK_CUR);
                fprintf(Fs,") %s ",U_And);
            }
            fseek(Fs,-5,1);
            fprintf(Fs,"\n   %s %s = ",U_Then,OutName[j]);
            value = R->readTSK(0);
            sprintf(temp,"%%.%dg",accuracy);
            if (value) fprintf(Fs,temp,value);
            for (i=0;i<in;i++)
            {
                value = R->readTSK(i+1);
                sprintf(temp," %%+ .%dg*%%s",accuracy);
                if (value!=0)
                    fprintf(Fs,temp,value,InName[i]);
            }
            fprintf(Fs,"\n");
        }
    }
    fclose(Fs);
    delete [] N;
    for (i=0;i<in;i++)
        delete [] InName[i];
    delete [] InName;
    for (i=0;i<out;i++)
        delete [] OutName[i];
    delete [] OutName;
    return 0;
}

Array1D<double> FLT::evaluate(const double * const point, System &S, System &C)
{
    const size_t n = S.outputs(); // Order
    const size_t m = C.outputs(); // Number of control inputs

	size_t  M_max = 0, N_max = 0;
    size_t *M = new size_t[n]; // Number of rules for each output of the plant
    for (size_t i=0;i<n;i++)
    {
        *(M+i) = S.rules(i);
        if (M[i]>M_max)
            M_max = M[i];
    }
    size_t *N = new size_t[m]; // Number of rules for each output of the controller
    for (size_t j=0;j<m;j++)
    {
        *(N+j) = C.rules(j);
        if (N[j]>N_max)
            N_max = N[j];
    }

    // Check if the plant and the controller are coherent
    if (S.outputs() != C.inputs())
    {
        delete [] M;
        delete [] N;
        return Array1D<double>();
    }
    if (n != S.inputs()-C.outputs())
    {
        delete [] M;
        delete [] N;
        return Array1D<double>();
    }
	
	if (S.checkInputLimits(point))
		WARNINGMSG(W_InputsLimitsExceeded)

    // State vector extension
    double *X = new double [n+1];
    X[0] = 1.0;
    for (size_t i=0;i<n;i++)
        *(X+i+1) = *(point+i);

    double temp, at, Sum_activ;
    double *Sum_W = new double[n];
    double **activ = new double *[N_max]; // Matching degree of the rules of the controller (omega)
    double **ct = new double *[n+1]; // Variable coefficient for the controller
    double *u = new double [m]; // Control signal
    double **W = new double *[M_max]; // Matching degree of the rules of the plant
    double *bt = new double [m]; // Variable coefficient 'b' for the plant
    Rule *R;
	Membership *P;

    for (size_t l=0;l<M_max;l++)
        W[l] = new double[n];
    for (size_t r=0;r<N_max;r++)
        activ[r] = new double[m];
    for (size_t i=0;i<=n;i++)
        ct[i] = new double[m];
		
    // Calculation of dX 
    for (size_t j=0;j<m;j++)
    {
        Sum_activ = 0.0;
        for (size_t r=0;r<N[j];r++)
        {
            R = C.readRule(j,r);
            activ[r][j] = 1.0;
            for (size_t k=0;k<n;k++)
			{
				P = R->readFunction(k);
				if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
					activ[r][j] *= P->eval(point[k]);
			}
            Sum_activ +=activ[r][j];
        } // activ[r][j]
		
        u[j] = 0.0;
        for (size_t k=0;k<=n;k++)
        {
            ct[k][j] = 0.0;
            for (size_t r=0;r<N[j];r++)
                ct[k][j] += activ[r][j]*C.readRule(j,r)->readTSK(k);
            ct[k][j] /= (Sum_activ+M_EPS); // ct[k][j]
            u[j] += X[k] * ct[k][j];
        } // u[j]
    }

    for (size_t i=0;i<n;i++)
    {
        Sum_W[i] = 0.0;
        for (size_t l=0;l<M[i];l++)
        {
            W[l][i] = 1.0;
            R = S.readRule(i,l);
            for (size_t k=0;k<n;k++)
			{
				P = R->readFunction(k);
				if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(point[k]);
			}
            for (size_t j=0;j<m;j++)
			{
				P = R->readFunction(j+n);
				if (P->type()!=ANYMF) // ANYMF doesn't affect to the system
					W[l][i] *= P->eval(u[j]);
			}
            // W[l][i]
            Sum_W[i] += W[l][i];
        }
    }
	
	Array1D<double> dX(n, 0.0);
    for (size_t i=0;i<n;i++)
    {
        for (size_t j=0;j<m;j++)
        {
            bt[j] = 0.0;
            for (size_t l=0;l<M[i];l++)
                bt[j] += W[l][i]*S.readRule(i,l)->readTSK(j+n+1);
//            bt[j]=bt[j]/Sum_W[i]; // bt[j][i] // The division is made at the end
        }
		
        for (size_t k=0;k<=n;k++)
        {
            at = 0.0;
            for (size_t l=0;l<M[i];l++)
                at += W[l][i]*S.readRule(i,l)->readTSK(k);
//            at=at/Sum_W[i]; // at[k][i] // The division is made at the end
            // Feedback loop
            temp = 0.0;
            for (size_t j=0;j<m;j++)
                temp += bt[j]*ct[k][j];
            dX[i] += (at+temp)*X[k];
        }
        dX[i] /= (Sum_W[i]+M_EPS);

    }
	
//	if (S.checkOutputLimits(&dX[0]))
//		WARNINGMSG(W_OutputsLimitsExceeded)
	
    // Free the dynamic memory
    delete [] M;
    delete [] N;
    delete [] bt;
    for (size_t i=0;i<=n;i++)
        delete [] ct[i];
    delete [] ct;
    for (size_t r=0;r<N_max;r++)
        delete [] activ[r];
    delete [] activ;
    for (size_t l=0;l<M_max;l++)
        delete [] W[l];
    delete [] W;
    delete [] u;
    delete [] Sum_W;
    delete [] X;

    return dX;
}

System FLT::initialController(System &Plant)
{
	size_t n = Plant.outputs(); // Number of inputs for the Controller
	size_t m = Plant.inputs() - n; // Number of outputs for the Controller
	if (m<=0)
		return System();

	size_t *num_rules_C = new size_t[m];
	num_rules_C[0] = 0;
	for (size_t i=0; i<n; i++)
		num_rules_C[0] += Plant.rules(i); // The controller will have all the systems rules in each output
	for (size_t j=1; j<m; j++)
		num_rules_C[j] = num_rules_C[0];
	
	System C(n, m, num_rules_C);
	
	delete []num_rules_C;
	
	Rule *pPrule; // Rule to save the rules of the Plant
	Rule Crule(n); // Rule to prepare the rules of the Controller
	Membership *mP;
	for (size_t i=0; i<=n; i++)
		Crule.changeTSK(0.0, i); // Initialize all consequents = 0

	// Output Limits
	double *CminOut = new double[m];
	double *CmaxOut = new double[m];
	
	for (size_t j=0; j<m; j++)
	{
		size_t r_idx = 0;
		for (size_t i=0; i<n; i++) // Only the 'n' first inputs are used in the Controller
		{
			for (size_t r=0; r<Plant.rules(i); r++)
			{
				pPrule = Plant.readRule(i, r);
				for (size_t q=0; q<n; q++) // Only the 'n' first functions are used in the Controller
				{
					mP = pPrule->readFunction(q);
					if(Crule.changeFunction(*mP, q))
					{
						delete []CminOut;
						delete []CmaxOut;
						#ifndef MATLAB_MEX_FILE
							cout<<E_MFDef<<endl;
						#endif
						#ifdef MATLAB_MEX_FILE
							mexPrintf("%s\n\n",E_MFDef);
						#endif
						return System();
					}
				}
				if(C.changeRule(Crule, r_idx, j))
				{
					delete []CminOut;
					delete []CmaxOut;
					#ifndef MATLAB_MEX_FILE
						cout<<E_MFDef<<endl;
					#endif
					#ifdef MATLAB_MEX_FILE
						mexPrintf("%s\n\n",E_MFDef);
					#endif
					return System();
				}
				r_idx++;
			}
		}
		// Copy the limits for the outputs of the Controller
		CminOut[j] = Plant.in_low(j + n);
		CmaxOut[j] = Plant.in_high(j + n);
	}
	
	// Input Limits
	double *CminIn = new double[n];
	double *CmaxIn = new double[n];

	for (size_t i=0; i<n; i++)
	{
		// Copy the limits for the inputs of the Controller
		CminIn[i] = Plant.in_low(i);
		CmaxIn[i] = Plant.in_high(i);
	}
	// Change the limits
	C.in_low(CminIn);
	C.in_high(CmaxIn);
	C.out_low(CminOut);
	C.out_high(CmaxOut);

	delete []CminIn;
	delete []CmaxIn;
	delete []CminOut;
	delete []CmaxOut;
	
	return C;
}

System FLT::subSystem(System &S, size_t nrules, size_t *outputs, size_t *rules)
{
	size_t idx;
    size_t n = S.inputs();
    size_t *selectOutputs = new size_t[S.outputs()];
	size_t *old2newOutputs = new size_t [S.outputs()];
    Rule *R;
    for (size_t i=0;i<S.outputs();i++)
    {
	    selectOutputs[i] = 0;
		old2newOutputs[i] = -1;
	}
	size_t m = 0;
	for (size_t i=0;i<nrules;i++)
	{
		if (selectOutputs[outputs[i]]==0)
		{
			selectOutputs[outputs[i]] = 1;
			m++;
		}
	}
	for (size_t i=0,idx=0;i<S.outputs();i++)
	{
		if (selectOutputs[i]!=0)
		{
			old2newOutputs[i] = idx;
			idx++;
		}
	}
    size_t *N = new size_t[m];
    for (size_t i=0;i<m;i++)
        N[i] = 0;
    for (size_t i=0;i<nrules;i++)
        N[(old2newOutputs[(outputs[i])])]++; // Count the number of rules for each output
	
    System Subsystem(n,m,N);
    for (size_t j=0;j<m;j++)
        N[j] = 0; // N[j] will be a counter in the next for
	
	for (size_t i=0;i<n;i++)
	{
		Array1D<double> limits = S.in_low();
		Subsystem.in_low(&limits[0]);
		limits = S.in_high();
		Subsystem.in_high(&limits[0]);
	}
	
	Array1D<double> Llimits(m);
	Array1D<double> Hlimits(m);
    for (size_t i=0;i<nrules;i++)
    {
		R = S.readRule(*outputs, *rules);
        if(!R)
        {
            delete []N;
			delete []old2newOutputs;
            delete []selectOutputs;
            return System();
        }		
        size_t out = old2newOutputs[(*outputs)];
		Llimits[out] = S.out_low(*outputs);
		Hlimits[out] = S.out_high(*outputs);
	
        if(Subsystem.changeRule(*R,N[out],out))
        {
            delete []N;
			delete []old2newOutputs;
            delete []selectOutputs;
            return System();
        }
        N[out]++; // Increment the rules counter
        outputs++;
        rules++;
    }
	
	Subsystem.out_low(&Llimits[0]);
	Subsystem.out_high(&Hlimits[0]);

	delete []N;
	delete []old2newOutputs;
    delete []selectOutputs;
	
    return Subsystem;
}

Array2D<double> FLT::extractPoints(System &S, unsigned int numMeshPoints, double precision, bool addMesh, bool onlyStatePoints)
{
	if (numMeshPoints<2)
	{
		cout<<"Warning, numMeshPoints in extractPoints must be >1, changed to 2.\n"<<endl;
		numMeshPoints = 2;
	}
	
	using namespace TNT;
	size_t n;
	if (onlyStatePoints)
		n = S.outputs();
	else
		n = S.inputs();
		
	size_t totalNumRules = 0;
	for (size_t j=0;j<S.outputs();j++) // Calculate the total number of rules of the system
		for (size_t r=0;r<S.rules(j);r++)
			totalNumRules += S.rules(j);
			
	Array2D<size_t> numPoints(totalNumRules,n,(size_t)0);
	
	for (size_t rt=0,j=0;j<S.outputs();j++)
	{// Calculate the number of point for each input
		for (size_t r=0;r<S.rules(j);r++,rt++)
		{
			Rule *R=S.readRule(j,r);
			for (size_t i=0;i<n;i++)
			{
				if (addMesh)
					numPoints[rt][i] += numMeshPoints;
					
				Membership *M = R->readFunction(i);
				switch(M->type())
				{
					case ANYMF:
					case CONSTMF:
					case SMF:
					case SIGMF:
					case ZMF:
					{
						numPoints[rt][i] += numMeshPoints;
						break;
					}
					case GAUSSMF:
					case GBELLMF:
					case SIG2MF:
					case TRAPMF:
					case TRIMF:
					case PIMF:
					case PSIGMF:
					{
						numPoints[rt][i] += 3;
						break;
					}
					case GAUSS2MF:
					{
						numPoints[rt][i] += 4;
						break;
					}
					default:
						Array2D <double> points;
						return points;
				}
			}
		}
	}
	
	size_t totalNumPoints = 0;
	for (size_t rt=0;rt<totalNumRules;rt++)
	{
		size_t temp = 1;
		for (size_t i=0;i<n;i++)
			temp *= numPoints[rt][i];
		totalNumPoints += temp;
	}
	
	Array2D<double>Points(n,totalNumPoints);
	size_t actualPoint = 0;
	for (size_t rt=0,j=0;j<S.outputs();j++)
	{
		for (size_t r=0;r<S.rules(j);r++,rt++)
		{// Extract the points of each rule
			Rule *R = S.readRule(j,r);
			size_t totalRepeat = 1;
			for (size_t i=0;i<n;i++)
					totalRepeat *= numPoints[rt][i];
			for (size_t i=0;i<n;i++)
			{
				Array1D<double> iPoints(numPoints[rt][i]); // The points of the input i
				Membership *M = R->readFunction(i);
				switch(M->type())
				{
					case ANYMF:
					case CONSTMF:
					{
						double d = S.in_low(i);
						double inc = (S.in_high(i) - d) / (numMeshPoints - 1.0);
						for (size_t p=0;p<numMeshPoints;p++,d+=inc)
							iPoints[p] = d;
						break;
					}
					case GAUSSMF:
					{
						iPoints[0] = M->read(0);
						iPoints[1] = M->read(0) - 0.83*M->read(1);
						iPoints[2] = M->read(0) + 0.83*M->read(1);
						break;
					}
					case GAUSS2MF:
					{
						//iPoints[0] = M->read(0)+(M->read(2)-M->read(0))/2; // If prefer only the middle of centers (3 points)
						iPoints[0] = M->read(0);
						iPoints[1] = M->read(2);
						iPoints[2] = M->read(0) - 0.83*M->read(1);
						iPoints[3] = M->read(2) + 0.83*M->read(3);
						break;
					}
					case GBELLMF:
					{
						iPoints[0] = M->read(2);
						iPoints[1] = M->read(2) - M->read(0);
						iPoints[2] = M->read(2) + M->read(0);
						break;
					}
					case PIMF:
					{
						iPoints[0] = M->read(0) + (M->read(1) - M->read(0)) / 2.0;
						iPoints[1] = M->read(2) + (M->read(3) - M->read(2)) / 2.0;
						iPoints[2] = iPoints[0] + (iPoints[1] - iPoints[0]) / 2.0;
						break;
					}
					case SMF:
					{
						double d = M->read(0) + (M->read(1) - M->read(0))/2.0;
						double inc = (S.in_high(i) - d) / (numMeshPoints - 1.0);
						for (size_t p=0;p<numMeshPoints;p++,d+=inc)
							iPoints[p] = d;
						break;
					}
					case SIGMF:
					{
						double d = M->read(1);
						double inc = (S.in_high(i) - d) / (numMeshPoints - 1.0);
						for (size_t p=0;p<numMeshPoints;p++,d+=inc)
							iPoints[p] = d;
						break;
					}
					case PSIGMF:
					case SIG2MF:
					{
						iPoints[0] = M->read(1);
						iPoints[1] = M->read(3);
						iPoints[2] = iPoints[0] + (iPoints[1] - iPoints[0]) / 2.0;
						break;
					}
					case TRAPMF:
					{
						iPoints[0] = M->read(1) + (M->read(2) - M->read(1)) / 2.0;
						iPoints[1] = iPoints[0] - (M->read(1) - M->read(0)) / 2.0;
						iPoints[2] = iPoints[0] + (M->read(3) - M->read(2)) / 2.0;
						break;
					}
					case TRIMF:
					{
						iPoints[0] = M->read(1);
						iPoints[1] = iPoints[0] - (M->read(1) - M->read(0)) / 2.0;
						iPoints[2] = iPoints[0] + (M->read(2) - M->read(1)) / 2.0;
						break;
					}
					case ZMF:
					{
						double d = M->read(1) - (M->read(1) - M->read(0)) / 2.0;
						double dec = (S.in_low(i) - d) / (numMeshPoints - 1.0);
						for (size_t p=0;p<numMeshPoints;p++,d-=dec)
							iPoints[p] = d;
						break;
					}
					default:
						Array2D <double> points;
						return points;
				}
				
				if (addMesh)
				{
					size_t dim = iPoints.dim();
					double d = S.in_low(i);
					double inc = (S.in_high(i) - d) / (numMeshPoints - 1.0);
					for (size_t p=dim-numMeshPoints;p<dim;p++,d+=inc)
						iPoints[p] = d;
				}
				
				// Fill Points with iPoints
				size_t repeat = 1;
				for (size_t i1=i+1;i1<n;i1++)
					repeat *= numPoints[rt][i1];
				size_t irepeat = 0;
				while (irepeat < totalRepeat)
				{
					for (size_t p=0;p<iPoints.dim();p++)
						for (size_t i1=0;i1<repeat;i1++,irepeat++)
							Points[i][irepeat+actualPoint] = iPoints[p];
				}
			}
			actualPoint += totalRepeat;
		}
	}
	
	// Locate repeated and out of limits points
	size_t numExcluded = 0;
	Array1D <bool> Excluded(totalNumPoints,false);
	for (size_t p=0;p<totalNumPoints;p++)
	{
		if (Excluded[p])
			continue;
		for (size_t pp = p + 1;pp<totalNumPoints;pp++)
		{
			if (Excluded[pp])
				continue;
			bool toExclude = false;
			double norm = 0.0;
			for (size_t i=0;(i<n) && !toExclude;i++)
			{ // Exclude out of limits points
				if ( (Points[i][pp]<S.in_low(i)) || (Points[i][pp]>S.in_high(i)) )
					toExclude = true;
				norm += fabs (Points[i][pp]);
			}
			if (norm < precision) // Exclude near the origin points
				toExclude = true;
			if (!toExclude)
			{
				toExclude=true;
				for (size_t i=0;(i<n) && toExclude;i++)
				{ // Exclude repeated points
					if (fabs(Points[i][p] - Points[i][pp]) > precision)
						toExclude = false;
				}
			}
			
			if (toExclude)
			{
				Excluded[pp] = true;
				numExcluded++;
			}
		}
	}

	// Remove excluded points
	Array2D<double> diferentPoints(n,totalNumPoints-numExcluded);
	for (size_t p=0,j=0;p<totalNumPoints;p++)
	{
		if (!Excluded[p])
		{
			for (size_t i=0;i<n;i++)
				diferentPoints[i][j] = Points[i][p];
			j++;
		}
	}

	return diferentPoints;
}
