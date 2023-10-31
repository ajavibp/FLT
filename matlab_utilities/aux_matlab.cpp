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
		if (!mxIsClass(model, "sugfis"))
		{
			mexPrintf("%s of class ''sugfis''\n\n",E_NoFIS);
			return 3;
		}

		if(FIS2System(model,S))
		{
			mexPrintf("%s. FIS2System fails.\n\n",E_NoFIS);
			return 4;
		}

	}
	return 0;
}

// FIS2System #################################################################
int FLT::FIS2System(const mxArray *FIS, System &S)
{
	mxArray *temp, *inputs, *outputs, *rules, *char_array[1];
	char text[MAX_CHAR_LONG];
	size_t n,m,M_in,M_out,M_R,j,i;

	// --- Initial checks ---

	// Inference must be PROD
	temp = mxGetProperty(FIS,0,"AndMethod");
	
	mexCallMATLAB(1, char_array, 1, &temp, "char");	mxGetString(char_array[0],text,MAX_CHAR_LONG); // Dirty way to convert String to char array (MATLAB API can't manager Strings objetcs)
	if (strcasecmp(text,"PROD"))
	{
		mexPrintf("%s\n\n",E_NoProd);
		return 1;
	}

	// Defuzzification method must be WTAVER
	temp = mxGetProperty(FIS,0,"DefuzzificationMethod");

	mexCallMATLAB(1, char_array, 1, &temp, "char");	mxGetString(char_array[0],text,MAX_CHAR_LONG); // Dirty way to convert String to char array

	if (strcasecmp(text,"WTAVER"))
	{
		mexPrintf("%s\n\n",E_NoSumProm);
		return 1;
	}

	// Reads model parameters

	inputs = mxGetProperty(FIS,0,"Inputs");
	outputs = mxGetProperty(FIS,0,"Outputs");
	rules = mxGetProperty(FIS,0,"Rules");

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
		temp = mxGetProperty(rules,l,"Weight");
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

		temp = mxGetProperty(rules,l,"Connection"); // Connection must be 1
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

		consequent = mxGetProperty(rules,l,"Consequent"); // Get the consequent (struct)
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
		limits = mxGetProperty(outputs,j,"Range");
		low_out[j] = *mxGetPr(limits);
		high_out[j] = *(mxGetPr(limits) + 1);
	}

	for (i=0;i<M_in;i++)
	{
		limits = mxGetProperty(inputs,i,"Range");
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
		antecedent = mxGetProperty(rules,l,"Antecedent");
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

		consequent = mxGetProperty(rules,l,"Consequent"); // Get the consequent (struct)
		for (i=0;i<M_out;i++)
		{
			v_conseq = (int)*(mxGetPr(consequent) + i);
			v_conseq--; // Number of consequent for the rule
			if (v_conseq>=0)
			{
				TSK = mxGetProperty(outputs,i,"MembershipFunctions"); // Get the consequent of the output 'i' (struct) // MembershipFunctions
				
				temp = mxGetProperty(TSK,v_conseq,"type");
				mexCallMATLAB(1, char_array, 1, &temp, "char"); mxGetString(char_array[0],TSKtype,MAX_CHAR_LONG); // Dirty way to convert String to char array

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
							mexPrintf("%s %s %lu, %s %lu\n\n",E_MFTypeDef,l+1,p+1,O_Rule,O_Input);
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
						mf_antec = mxGetProperty(inputs,p,"MembershipFunctions"); // Get the antecedent of the input 'p' (struct) // Mf

						temp = mxGetProperty(mf_antec,v_antec,"Type");
						mexCallMATLAB(1, char_array, 1, &temp, "char"); mxGetString(char_array[0],type,MAX_CHAR_LONG); // Dirty way to convert String to char array

						for (int iT=0;iT<=(int)ZMF;iT++) // Check if the MF type
						{
							T = (TYPE_MF)iT;
							if (!strcasecmp(MF_NAMES[T],type))
								break;
						}
						if (T>ZMF)
						{
							mexPrintf("%s %s %lu, %s %lu\n\n",E_BadMF,O_Rule,l+1,O_Input,p+1);
							mxDestroyArray(inputs);
							mxDestroyArray(outputs);
							mxDestroyArray(rules);
							mxDestroyArray(antecedent);
							mxDestroyArray(consequent);
							delete [] num_rules;
							return 1;
						}

						temp = mxGetProperty(mf_antec,v_antec,"Parameters");

						if (MF_PARAM_NUMBER[(int)T]>mxGetN(temp)) 
						{
							mexPrintf("%s %s %lu. %s %lu\n\n",E_NumberParam,O_MF,p+1,O_Rule,l+1);
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
							mexPrintf("%s %s %lu, %s %lu\n\n",E_MFTypeDef,l+1,p+1,O_Rule,O_Input);
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
							d=*(mxGetPr(temp)+v); // temp is 'mf_antec' before
							
							
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
							temp = mxGetProperty(TSK,v_conseq,"Parameters");
							d = *(mxGetPr(temp)+v);
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
								d = *(mxGetPr(temp)+v); // temp is TSK before
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

	char type[MAX_CHAR_LONG];
	
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


	// Create the input struct -----
	char inputName[MAX_CHAR_LONG], nom_mf[MAX_CHAR_LONG];
	mxArray *range, *mf, *param;
	double d;
	TYPE_MF temp;

	mxArray *inputs = mxGetProperty(FIS,0,"Inputs"); // Get the input
	if(!inputs)
	{
		mexPrintf("Error reading ''Input'' fisvar.\n");
		delete [] M;
		return 1;
	}
	for (i=0;i<in;i++)
	{
		sprintf(&inputName[0],"%s%lu",U_FISInput,i+1);

		mxSetProperty(inputs, i, "Name", mxCreateString(inputName)); // Change the name

		range = mxGetProperty(inputs,i,"Range"); // Get the range vector
		if(!range)
		{
			mexPrintf("Error reading inputs ''Range'' values.\n");
			delete [] M;
			return 1;
		}
		*(mxGetPr(range)) = S.in_low(i);
		*(mxGetPr(range)+1) = S.in_high(i); // Change the range vector

		mxSetProperty(inputs,i,"Range",range); // Save the new range vector
		
		mf = mxGetProperty(inputs,i,"MembershipFunctions"); // Get the Membership function structure
		if(!mf)
		{
			mexPrintf("Error reading inputs ''MembershipFunctions'' values.\n");
			delete [] M;
			return 1;
		}

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
				mxSetProperty(mf,index,"Name",mxCreateString(nom_mf)); // Save the new name of the Membership Function

				temp = P->type();
				for (p=0;p<=strlen(MF_NAMES[(int)temp]);p++)
					type[p] = tolower(MF_NAMES[(int)temp][p]);
					
				mxSetProperty(mf,index,"Type",mxCreateString(type)); // Save the new type of the Membership Function

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
				
				mxSetProperty(mf,index,"Parameters",param); // Save the new type of the Membership Function
			}
		}
		mxSetProperty(inputs,i,"MembershipFunctions",mf);
	}
	mxSetProperty(FIS,0,"Inputs",inputs);
	

	// Create the output struct -----
	mxArray *outputs = mxGetProperty(FIS,0,"Outputs"); // Get the output
	if(!outputs)
	{
		mexPrintf("Error reading ''Input'' fisvar.\n");
		delete [] M;
		return 1;
	}

	char outputName[MAX_CHAR_LONG];
	for (j=0;j<out;j++)
	{
		sprintf(&outputName[0],"%s%lu",U_FISOutput,j+1);
		mxSetProperty(outputs, j, "Name", mxCreateString(outputName)); // Change the name
		range = mxGetProperty(outputs,j,"Range"); // Get the range vector
		if(!range)
		{
			mexPrintf("Error reading output ''Range'' values.\n");
			delete [] M;
			return 1;
		}
		*(mxGetPr(range)) = S.out_low(j);
		*(mxGetPr(range)+1) = S.out_high(j); // Change the range vector
		mxSetProperty(outputs,j,"Range",range); // Save the new range vector

		mf = mxGetProperty(outputs,j,"MembershipFunctions"); // Get the Membership function structure
		if(!mf)
		{
			mexPrintf("Error reading ''MembershipFunctions'' values.\n");
			delete [] M;
			return 1;
		}
		for (r=0;r<M[j];r++)
		{
			R = S.readRule(j,r);
			sprintf(&nom_mf[0],"%s%lu%s%lu",U_RuleAbbreviation,r+1,U_OutputAbreviation,j + 1);
			mxSetProperty(mf,r,"Name",mxCreateString(nom_mf)); // Save the new name of the Membership Function
			mxSetProperty(mf,r,"Type",mxCreateString("linear")); // Save the new type of the Membership Function
			
			param = mxCreateDoubleMatrix(1,in+1,mxREAL);
			for (i=0;i<in;i++)
				*(mxGetPr(param)+i) = R->readTSK(i+1);
			*(mxGetPr(param)+in) = R->readTSK(0);
			mxSetProperty(mf,r,"Parameters",param); // Save the new type of the Membership Function
		}
		mxSetProperty(outputs,j,"MembershipFunctions",mf);
	}
	mxSetProperty(FIS,0,"Outputs",outputs);

	// Create the rule struct -----
	index=0;
	for (j=0;j<out;j++)
	{
		for (r=0;r<M[j];r++,index++)
		{
			mxArray *newRule = mxCreateDoubleMatrix(1,in+out+2,mxREAL);
			mxDouble *pointerNewRule = mxGetPr(newRule);
			for (i=0;i<in;i++,pointerNewRule++)
				*pointerNewRule = index + 1; // Antecedents (+1 because MATLAB arrays start in 1)
		
			*(pointerNewRule + j) = r + 1; // Consequents (+1 because MATLAB arrays start in 1)
			pointerNewRule += out;

			*pointerNewRule = 1; // Weight
			pointerNewRule++;
			*pointerNewRule = 1; // Connection

			mxArray *newFIS;
			mxArray *parameters[2];
			parameters[0] = FIS;
			parameters[1] = newRule;
			mexCallMATLAB(1,&newFIS,2,parameters,"addRule"); // I don't like using mexCallMATLAB here, but I don't know how to do it any other way.
			mxSetProperty(FIS,0,"Rules",mxGetProperty(newFIS,0,"Rules"));
		}
	}
	
	delete [] M;
	return 0;
}
