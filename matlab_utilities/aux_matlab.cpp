#include "aux_matlab.hpp"

using namespace FLT;

// readModel ##################################################################
int FLT::readModel(const mxArray *model, System &S)
{
	if (mxIsChar(model))
	{// The model is a TXT or FIS file
		char name[MAX_CHAR_LONG], extension[4];
		bool end = false;
		
		if (mxGetString(model,name,MAX_CHAR_LONG))
			return 1;
		
		for (size_t i=0;(i<MAX_CHAR_LONG && end==0);i++)
		{
			if(name[i]=='\0')
			{
				i = i-3;
				for (size_t j=0;j<4;j++)
					extension[j] = name[i+j];
				end = true;
			}
		}
		if (!end)
			return 1;
			
		if(!strcasecmp(extension,"FIS"))
		{// The model is a FIS file
			mxArray *right_arg[1], *left_arg[1];
			
			right_arg[0] = mxCreateString(name);
			if(mexCallMATLAB(1,left_arg,1,right_arg,"readfis"))
			{
				mxDestroyArray(right_arg[0]);
				return 1;
			}
			if(FIS2System(left_arg[0],S))
			{
				mxDestroyArray(right_arg[0]);
				return 1;
			}
			mxDestroyArray(left_arg[0]);
			mxDestroyArray(right_arg[0]);
		}
		else if(!strcasecmp(extension,"TXT"))
		{// Is TXT file
			S = TXT2System(&name[0]);
			if (!S.outputs())
				return 2;
		}
		else
			return 1;
	}
	else
	{
		// Is a FIS variable
		if (!mxIsStruct(model) || (mxGetNumberOfFields(model)<10 && mxGetNumberOfFields(model)>13))
			return 3;
		if(FIS2System(model,S))
			return 4;
	}
	return 0;
}

// FIS2System #################################################################
int FLT::FIS2System(const mxArray *FIS, System &S)
{
	mxArray *temp, *inputs, *outputs, *rules;
	char text[MAX_CHAR_LONG];
	size_t n,m,M_in,M_out,M_R,j,i;

	// Initial checks
	i = (size_t)mxGetNumberOfFields(FIS);
	if (i<10 || i>13)
	{
		mexPrintf("%s\n",E_NoFIS);
		mexPrintf("%s\n\n",U_SeeManual);
		return 1;
	}

	// Type must be SUGENO
	temp = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"type"));
	if (mxGetString(temp,text,MAX_CHAR_LONG))
	{
		mexPrintf("%s\n",E_Model);
		mexPrintf("%s\n",E_NoFIS);
		mexPrintf("%s\n\n",U_SeeManual);
		return 1;
	}
	if (strcasecmp(text,"SUGENO"))
	{
		mexPrintf("%s\n\n",E_NoSugeno);
		return 1;
	}
	
	// Inference must be PROD
	temp = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"andMethod"));
	if (mxGetString(temp,text,MAX_CHAR_LONG))
	{
		mexPrintf("%s.\n",E_Model);
		mexPrintf("%s\n",E_NoFIS);
		mexPrintf("%s\n\n",U_SeeManual);
		return 1;
	}
	if (strcasecmp(text,"PROD"))
	{
		mexPrintf("%s\n\n",E_NoProd);
		return 1;
	}
	
	// Defuzzification method must be WTAVER
	temp = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"defuzzMethod"));
	if (mxGetString(temp,text,MAX_CHAR_LONG))
	{
		mexPrintf("%s\n",E_Model);
		mexPrintf("%s\n",E_NoFIS);
		mexPrintf("%s\n\n",U_SeeManual);
		return 1;
	}
	if (strcasecmp(text,"WTAVER"))
	{
		mexPrintf("%s\n\n",E_NoSumProm);
		return 1;
	}

	// Reads model parameters

	inputs = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"input")); // Inputs struct
	outputs = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"output")); // Outputs struct
	rules = mxGetFieldByNumber(FIS,0,mxGetFieldNumber(FIS,"rule")); // Rules struct

	M_in = (size_t)mxGetN(inputs); // Number of inputs
	M_out = (size_t)mxGetN(outputs); // Number of outputs
	M_R = (size_t)mxGetN(rules); // Total number of rules (the sum of the rules of all outputs)

	if (M_out==0)
	{
		mexPrintf("%s\n\n",E_NoOutputs);
		mxDestroyArray(inputs);
		mxDestroyArray(outputs);
		mxDestroyArray(rules);
		return 2;
	}
	
	if (M_in==0)
	{
		mexPrintf("%s\n\n",E_NoInputs);
		mxDestroyArray(inputs);
		mxDestroyArray(outputs);
		mxDestroyArray(rules);
		return 2;
	}
	
	if (M_R<M_out)
	{
		mexPrintf("%s\n\n",E_NoRules);
		mxDestroyArray(inputs);
		mxDestroyArray(outputs);
		mxDestroyArray(rules);
		return 2;
	}

	// Creates the model and reads data
	size_t l,p,v,r;
	int v_antec, v_conseq, check_conseq, error;
	size_t *num_rules = new size_t[M_out];
	for (i=0;i<M_out;i++)
		*(num_rules+i)=0;
		
	double *low_in = new double[M_in];
	double *high_in = new double[M_in];
	double *low_out = new double[M_out];
	double *high_out = new double[M_out];
	double d;
	bool constant, linear;
	char type[ MAX_SIZE_TYPE_NAME];
	char TSKtype[ MAX_SIZE_TYPE_NAME];
	mxArray *antecedent, *consequent, *mf_antec, *TSK;
	TYPE_MF T = ANYMF;
	Membership *P;
	Rule *R;
	System Sist;

	// Counts and checks rules
	for (l=0;l<M_R;l++)
	{
		temp = mxGetFieldByNumber(rules,l,mxGetFieldNumber(rules,"weight")); // Weight must be 1
		if (mxGetScalar(temp)!=1)
		{
			mexPrintf("%s %lu %s 1.\n\n",O_RuleWeigh,l+1,O_DifferentOf);
			mxDestroyArray(inputs);
			mxDestroyArray(outputs);
			mxDestroyArray(rules);
			delete [] num_rules;
			delete [] low_in;
			delete [] high_in;
			delete []low_out;
			delete [] high_out;
			return 3;
		}
		
		temp = mxGetFieldByNumber(rules,l,mxGetFieldNumber(rules,"connection")); // Connection must be 1
		if (mxGetScalar(temp)!=1)
		{
			mexPrintf("%s %lu %s 1, AND.\n\n",O_RuleConnection,l+1,O_DifferentOf);
			mxDestroyArray(inputs);
			mxDestroyArray(outputs);
			mxDestroyArray(rules);
			delete [] num_rules;
			delete [] low_in;
			delete [] high_in;
			delete []low_out;
			delete [] high_out;
			return 4;
		}
		consequent = mxGetFieldByNumber(rules,l,mxGetFieldNumber(rules,"consequent")); // Get the consequent (struct)
		if (mxGetN(consequent)!=M_out)
		{
			mexPrintf("%s %lu %s",O_ConsequentSize,l+1,O_IsNotCorrect);
			mxDestroyArray(inputs);
			mxDestroyArray(outputs);
			mxDestroyArray(rules);
			delete [] num_rules;
			delete [] low_in;
			delete [] high_in; 
			delete []low_out;
			delete [] high_out;
			return 5;
		}
		check_conseq = 0; // Check that only 1 consequent for rule are defined
		for (i=0;i<M_out;i++)
		{
			v_conseq = (int)*(mxGetPr(consequent) + i);
			v_conseq--; // Number of consequent for the rule
			if (v_conseq>=0)
			{
				num_rules[i]++;
				if (check_conseq==1)
				{
					mexPrintf("%s\n",E_No1Conseq);
					mxDestroyArray(inputs);
					mxDestroyArray(outputs);
					mxDestroyArray(rules);
					delete [] num_rules;
					delete [] low_in;
					delete [] high_in;
					delete []low_out;
					delete [] high_out;
					return 6;
				}
				check_conseq = 1;
			}
		}
	}
	
	// Create the System
	Sist.initialize(M_in,M_out,num_rules);

	// Use 'num_rules' to index the rules of each output
	mxArray * limits;
	for (j=0;j<M_out;j++)
	{
		*(num_rules+j) = 0;
		limits = mxGetFieldByNumber(outputs,j,mxGetFieldNumber(outputs,"range")); // Get j-output range
		low_out[j] = *mxGetPr(limits);
		high_out[j] = *(mxGetPr(limits) + 1);
}

	for (i=0;i<M_in;i++)
	{
		limits = mxGetFieldByNumber(inputs,i,mxGetFieldNumber(inputs,"range"));// Get i-input range
		low_in[i] = *mxGetPr(limits);
		high_in[i] = *(mxGetPr(limits) + 1);
	}

	Sist.in_low(low_in);
	Sist.in_high(high_in);
	Sist.out_low(low_out);
	Sist.out_high(high_out);

	delete [] low_in;
	delete [] high_in;
	delete []low_out;
	delete [] high_out;

	for (l=0;l<M_R;l++)
	{
		antecedent = mxGetFieldByNumber(rules,l,mxGetFieldNumber(rules,"antecedent")); // Get the antecedent (struct)
		if (mxGetN(antecedent)!=M_in)
		{
			mexPrintf("%s %lu %s",O_AntecedentSize,l+1,O_IsNotCorrect);
			mxDestroyArray(inputs);
			mxDestroyArray(outputs);
			mxDestroyArray(rules);
			mxDestroyArray(antecedent);
			delete [] num_rules;
			return 1;
		}
		consequent = mxGetFieldByNumber(rules,l,mxGetFieldNumber(rules,"consequent")); // Get the consequent (struct)
		for (i=0;i<M_out;i++)
		{
			v_conseq = (int)*(mxGetPr(consequent) + i);
			v_conseq--; // Number of consequent for the rule
			if (v_conseq>=0)
			{
				TSK = mxGetFieldByNumber(outputs,i,mxGetFieldNumber(outputs,"mf")); // Get the consequent of the output 'i' (struct) // Mf
				if (mxGetString(mxGetFieldByNumber(TSK,v_conseq,mxGetFieldNumber(TSK,"type")),TSKtype, MAX_SIZE_TYPE_NAME)) // Consecuent type (TSK or Singleton)
				{
					mexPrintf("%s\n",E_Model);
					mexPrintf("%s %lu.\n\n",E_Conseq,l+1);
					mxDestroyArray(inputs);

					mxDestroyArray(outputs);
					mxDestroyArray(rules);
					mxDestroyArray(antecedent);
					mxDestroyArray(consequent);
					delete [] num_rules;
					return 1;
				}
				constant = !strcasecmp(TSKtype,"CONSTANT");
				linear = !strcasecmp(TSKtype,"LINEAR");
				for (p=0;p<M_in;p++)
				{
					R = Sist.readRule(i,num_rules[i]);
					if (R==NULL)
					{
						mexPrintf("%s\n",E_BadModel);
						mxDestroyArray(inputs);
						mxDestroyArray(outputs);
						mxDestroyArray(rules);
						mxDestroyArray(antecedent);
						mxDestroyArray(consequent);
						delete [] num_rules;
						return 1;
					}
					v_antec = (int)*(mxGetPr(antecedent) + p);
					v_antec--; // Number of antecedent for the rule
					if (v_antec<0) // For use of 'none' in MATLAB to select ANYMF membership function
					{
						T = ANYMF;
						if (R->changeFunction(T,p))
						{
							mexPrintf("%s %s %lu, %s %lu.\n\n",E_MFTypeDef,l+1,p+1,O_Rule,O_Input);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 1;
						}
						P = createMF(T);
					}
					else
					{
						mf_antec = mxGetFieldByNumber(inputs,p,mxGetFieldNumber(inputs,"mf")); // Get the antecedent of the input 'p' (struct) // Mf
						if (mxGetString(mxGetFieldByNumber(mf_antec,v_antec,mxGetFieldNumber(mf_antec,"type")),type, MAX_SIZE_TYPE_NAME)) // Read the membership function type
						{
							mexPrintf("%s\n",E_Model);
							mexPrintf("%s %lu, %s %lu.\n\n",E_MFType,l+1,O_Input,p+1);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 1;
						}
						for (int iT=0;iT<=(int)ZMF;iT++) // Check if the MF type
						{
							T = (TYPE_MF)iT;
							if (!strcasecmp(MF_NAMES[T],type))
								break;
						}
						if (T>ZMF)
						{
							mexPrintf("%s %s %lu, %s %lu.\n\n",E_BadMF,O_Rule,l+1,O_Input,p+1);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 1;
						}
						if (MF_PARAM_NUMBER[(int)T]!=mxGetN(mxGetFieldByNumber(mf_antec,v_antec,mxGetFieldNumber(mf_antec,"params")))) 
						{
							mexPrintf("%s %s %lu. %s %lu.\n\n",E_NumberParam,O_MF,p+1,O_Rule,l+1);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 1;
						}
						if (R->changeFunction(T,p))
						{
							mexPrintf("%s %s %lu, %s %lu.\n\n",E_MFTypeDef,l+1,p+1,O_Rule,O_Input);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 7;
						}
						P = createMF(T);
						for (v=0;v<P->num_params();v++)
						{
							d=*(mxGetPr(mxGetFieldByNumber(mf_antec,v_antec,mxGetFieldNumber(mf_antec,"params")))+v);
							// Change the parameters for use with MATLAB (change the definition of Gauss functions)
							if (T==GAUSSMF || T==GAUSS2MF)
							{
								if (v==0 || v==2)
								{
									d = M_SQRT2*d;
									error = P->edit(v+1,d);
								}
								else if (v==1 || v==3)
									error = P->edit(v-1,d);
							}
							else
								error = P->edit(v,d);
						}
						if (error || P->test()>0)
						{
							mexPrintf("%s\n\n",E_ParamMF);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							delete P;
							P = NULL;
							return 8;
						}
					}
					if(R->changeFunction(*P,p))
					{
						mexPrintf("%s\n\n",E_MFDef);
						mxDestroyArray(inputs);
						mxDestroyArray(outputs);
						mxDestroyArray(rules);
						mxDestroyArray(antecedent);
						mxDestroyArray(consequent);
						delete [] num_rules;
						delete P;
						P = NULL;
						return 9;
					}
					delete P;
					P = NULL;
					for (v=0;v<=M_in;v++)
					{
						if(linear)
						{
							d = *(mxGetPr(mxGetFieldByNumber(TSK,v_conseq,mxGetFieldNumber(TSK,"params")))+v);
							if (v!=M_in)
							{
								if(R->changeTSK(d,v+1))
								{
									mexPrintf("%s\n\n",E_ParamConseq);
									mxDestroyArray(inputs);
									mxDestroyArray(outputs);
									mxDestroyArray(rules);
									mxDestroyArray(antecedent);
									mxDestroyArray(consequent);
									delete [] num_rules;
									return 11;
								}
							}
							else
							{
								if(R->changeTSK(d,0))
								{
									mexPrintf("%s\n\n",E_ParamConseq);
									mxDestroyArray(inputs);
									mxDestroyArray(outputs);
									mxDestroyArray(rules);
									mxDestroyArray(antecedent);
									mxDestroyArray(consequent);
									delete [] num_rules;
									return 11;
								}
							}
						}
						else // constant
						{
							if (v!=0)
							{
								if(R->changeTSK(0.0,v))
								{
									mexPrintf("%s\n\n",E_ParamConseq);
									mxDestroyArray(inputs);
									mxDestroyArray(outputs);
									mxDestroyArray(rules);
									mxDestroyArray(antecedent);
									mxDestroyArray(consequent);
									delete [] num_rules;
									return 11;
								}
							}
							else
							{
								d = *(mxGetPr(mxGetFieldByNumber(TSK,v_conseq,mxGetFieldNumber(TSK,"params")))+v);
								if(R->changeTSK(d,0))
								{
									mexPrintf("%s\n\n",E_ParamConseq);
									mxDestroyArray(inputs);
									mxDestroyArray(outputs);
									mxDestroyArray(rules);
									mxDestroyArray(antecedent);
									mxDestroyArray(consequent);
									delete [] num_rules;
									return 11;
								}
							}
						}
					}
				}
				num_rules[i]++;
			}
		}
	}
	delete [] num_rules;
	S = Sist;
	return 0;
}

// System2FIS #################################################################
int FLT::System2FIS(System &S, mxArray *FIS)
{
	size_t i,j,r,p,index,L=0;
	const size_t in = S.inputs();
	const size_t out = S.outputs();
	size_t *M = new size_t[out];
	
	Membership *P;
	Rule *R;
	
	for (i=0;i<out;i++)
	{
		M[i] = S.rules(i);
		L += M[i];
		if (M[i]<1)
		{
			mexPrintf("%s\n",E_NoRules);
			delete [] M;
			return 1;
		}
	}
	char type[MAX_CHAR_LONG];

	if (out<1)
	{
		mexPrintf("%s\n",E_NoOutputs);
		delete [] M;
		return 2;
	}
	if (in<1)
	{
		mexPrintf("%s\n",E_NoInputs);
		delete [] M;
		return 3;
	}

	// Create the FIS struct
	mxSetFieldByNumber(FIS,0,0,mxCreateString(U_FISName));
	mxSetFieldByNumber(FIS,0,1,mxCreateString("sugeno"));
	mxSetFieldByNumber(FIS,0,2,mxCreateString("prod"));
	mxSetFieldByNumber(FIS,0,3,mxCreateString("probor"));
	mxSetFieldByNumber(FIS,0,4,mxCreateString("wtaver"));
	mxSetFieldByNumber(FIS,0,5,mxCreateString("prod"));
	mxSetFieldByNumber(FIS,0,6,mxCreateString("max"));

	// Create the input struct
	mwSize dim[] = {1,in};
	const char* InOutNames[] = {"name","range","mf"};
	const char* mfNames[] = {"name","type","params"};
	mxArray *inputs = mxCreateStructArray(2,dim,3,InOutNames);
	char char_num[MAX_CHAR_LONG], inputName[MAX_CHAR_LONG], nom_mf[MAX_CHAR_LONG];
	mxArray *range, *mf, *param;
	mxArray *outputs;
	double d;
	TYPE_MF temp;

	for (i=0;i<in;i++)
	{
		sprintf(&inputName[0],"%s%lu",U_FISInput,i+1);
		mxSetFieldByNumber(inputs,i,0,mxCreateString(inputName));
		range = mxCreateDoubleMatrix(1,2,mxREAL);
		*(mxGetPr(range)) = S.in_low(i);
		*(mxGetPr(range)+1) = S.in_high(i);
		mxSetFieldByNumber(inputs,i,1,range);
		dim[1] = L;
		mf = mxCreateStructArray(2,dim,3,mfNames);
		
		index=0;
		for (j=0;j<out;j++)
		{
			for (r=0;r<M[j];r++,index++)
			{
				R = S.readRule(j,r);
				if (R==NULL)
				{
					delete [] M;

					return 4;
				}
				P = R->readFunction(i);
				if (P==NULL || P->test()>0)
				{
					delete [] M;
					return 5;
				}
				
				sprintf(nom_mf,"%s%lu%s%lu%s%lu",U_RuleAbbreviation,r+1,U_InputAbbreviation,i+1,U_OutputAbreviation,j+1);
				mxSetFieldByNumber(mf,index,0,mxCreateString(nom_mf));
				temp = P->type();
				for (p=0;p<=strlen(MF_NAMES[(int)temp]);p++)
					type[p] = tolower(MF_NAMES[(int)temp][p]);
					
				mxSetFieldByNumber(mf,index,1,mxCreateString(type));
				param = mxCreateDoubleMatrix(1,P->num_params(),mxREAL); // Do not use MF_PARAM_NUMBER[(int)temp], because P->num_params() are more general
				for (p=0;p<P->num_params();p++)
					*(mxGetPr(param)+p) = P->read(p);

				// Change the parameters for use with MATLAB (change the definition of Gauss functions)
				if (temp==GAUSSMF || temp==GAUSS2MF)
				{
					d = mxGetScalar(param);
					(*(mxGetPr(param))) = (*(mxGetPr(param)+1))/M_SQRT2;
					*(mxGetPr(param) + 1) = d;
				}
				if (temp==GAUSS2MF)
				{
					d = *(mxGetPr(param) + 2);
					(*(mxGetPr(param) + 2)) = (*(mxGetPr(param)+3))/M_SQRT2;
					*(mxGetPr(param) + 3) = d;
				}
				mxSetFieldByNumber(mf,index,2,param);
			}
		}
		mxSetFieldByNumber(inputs,i,2,mf);
	}
	mxSetFieldByNumber(FIS,0,7,inputs);

	// Create the output struct
	dim[1] = out;
	outputs = mxCreateStructArray(2,dim,3,InOutNames);
	char outputName[MAX_CHAR_LONG];
	for (j=0;j<out;j++)
	{
		sprintf(&outputName[0],"%s%lu",U_FISOutput,j+1);
		mxSetFieldByNumber(outputs,j,0,mxCreateString(outputName));
		range = mxCreateDoubleMatrix(1,2,mxREAL);
		*(mxGetPr(range)) = S.out_low(j);
		*(mxGetPr(range) + 1) = S.out_high(j);
		mxSetFieldByNumber(outputs,j,1,range);
		dim[1] = M[j];
		mf = mxCreateStructArray(2,dim,3,mfNames);
		for (r=0;r<M[j];r++)
		{
			R = S.readRule(j,r);
			sprintf(&nom_mf[0],"%s%lu%s%lu",U_RuleAbbreviation,r+1,U_OutputAbreviation,j + 1);
			mxSetFieldByNumber(mf,r,0,mxCreateString(nom_mf));
			mxSetFieldByNumber(mf,r,1,mxCreateString("linear"));
			param = mxCreateDoubleMatrix(1,in+1,mxREAL);
			for (i=0;i<in;i++)
				*(mxGetPr(param)+i) = R->readTSK(i+1);
			*(mxGetPr(param)+in) = R->readTSK(0);
			mxSetFieldByNumber(mf,r,2,param);
		}
		mxSetFieldByNumber(outputs,j,2,mf);
	}
	mxSetFieldByNumber(FIS,0,8,outputs);

	// Create the rule struct
	dim[1] = L;
	const char* ruleNames[] = {"antecedent","consequent","weight","connection"};
	mxArray *rules = mxCreateStructArray(2,dim,4,ruleNames);
	mxArray *antecedent, *consequent;
	size_t input = 0, output = 0, value;
	index = 0;
	for (j=0;j<out;j++)
	{
		for (r=0;r<M[j];r++,index++)
		{
			antecedent = mxCreateDoubleMatrix(1,in,mxREAL);
			mxSetFieldByNumber(rules,index,0,antecedent);
			for (i=0;i<in;i++)
				*(mxGetPr(antecedent)+i) = index + 1;
			consequent = mxCreateDoubleMatrix(1,out,mxREAL);
			*(mxGetPr(consequent)+j) = (r + 1);
			mxSetFieldByNumber(rules,index,1,consequent);
			mxSetFieldByNumber(rules,index,2,mxCreateDoubleScalar(1)); // Weight=1
			mxSetFieldByNumber(rules,index,3,mxCreateDoubleScalar(1)); // Connection=1
		}
	}
	mxSetFieldByNumber(FIS,0,9,rules);
	delete [] M;
	return 0;
}
