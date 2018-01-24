/*
 * checkInput.cpp
 *
 *  Created on: Mar 18, 2012
 *
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#include <stdlib.h>
#include <iostream>

#include "checkInput.h"

using namespace std;

inputVCF::inputVCF(const vector<string> & FileNameVCF,string FileNamePed,string OutputName,string LocationFile,bool Multi,bool VarOnly,int Method,double MRate,bool PosOrder)
{
	fileNameVCF=FileNameVCF;
	fileNamePed=FileNamePed;
	outputName=OutputName;
	locationFile=LocationFile;
	//multi=Multi;
	varOnly=VarOnly;
	method=Method;
	mRate=MRate;
	posOrder=PosOrder;
	allLine=false;
}

inputVCF::inputVCF()
{
}

bool inputVCF::set_fileNameVCF(const vector<string> & FileNameVCF)
{
	fileNameVCF=FileNameVCF;
	return true;
}

bool inputVCF::set_fileNamePed(string FileNamePed)
{
	fileNamePed=FileNamePed;
	return true;
}

bool inputVCF::set_outputName(string OutputName)
{
	outputName=OutputName;
	return true;
}

bool inputVCF::set_locationFile(string LocationFile)
{
	locationFile=LocationFile;
	return true;
}

//bool inputVCF::set_multi(bool Multi)
//{
//	multi=Multi;
//	return true;
//}

bool inputVCF::set_varOnly(bool VarOnly)
{
	varOnly=VarOnly;
	return true;
}

bool inputVCF::set_method(int Method)
{
	method=Method;
	return true;
}

bool inputVCF::set_mRate(double MRate)
{
	mRate=MRate;
	return true;
}

bool inputVCF::set_genoProbN(const vector<double> & GenoProbN)
{
	genoProbN=GenoProbN;
	return true;
}

bool inputVCF::set_genoProbK(const vector<double> & GenoProbK)
{
	genoProbK=GenoProbK;
	return true;
}

bool inputVCF::set_genoProbXN(const vector<double> & GenoProbXN)
{
	genoProbXN=GenoProbXN;
	return true;
}

bool inputVCF::set_genoProbXK(const vector<double> & GenoProbXK)
{
	genoProbXK=GenoProbXK;
	return true;
}

bool inputVCF::set_posOrder(bool PosOrder)
{
	posOrder=PosOrder;
	return true;
}

bool inputVCF::set_allLine(bool AllLine)
{
	allLine=AllLine;
	return true;
}

bool inputVCF::set_diffOnly(bool DiffOnly)
{
	diffOnly = DiffOnly;
	if(diffOnly){
		allLine=false;
		varOnly=true;
	}
	return true;
}

bool inputVCF::set_numBurnIn(int NumBurnIn)
{
	numBurnIn=NumBurnIn;
	return true;
}

bool inputVCF::set_numRep(int NumRep)
{
	numRep=NumRep;
	return true;
}

bool inputVCF::set_LCriteria(double LC)
{
	LCriteria = LC;
	return true;
}

int checkInputVCF(int argc, char* argv[], inputVCF & vcfInfo)
{
	vector<string> fileNameVCF;
	string fileNamePed;
	string outputName;
	string locationFile;
	bool varOnly=false;
	int method=1;
	double mRate=1e-7;
	vector<double> genoProbN;
	vector<double> genoProbK;
	vector<double> genoProbXN;
	vector<double> genoProbXK;
	bool posOrder=false;
	bool allLine=false;
	bool diffOnly = false;
	int numBurnIn=-999;
	int numRep=-999;
	double lrc = 1;

	int returnVal=0;

	for(int i=2;i<argc;i++)
	{
		if(argv[i][0]=='-')
		{
			string strOption(argv[i]);
			string option=strOption.substr(1);

			//vcf
			if(option=="vcfFile")
			{
				bool flag=true;
				while(true)
				{
					i++;
					if(flag)
					{
						if(i==argc || argv[i][0]=='-')
						{
							cout<<"The vcf file hasn't been set. Please check input."<<endl;
							return -1;
						}
					}
					flag=false;
					if(i==argc || argv[i][0]=='-')
					{
						i--;
						break;
					}
					string vcfTmp(argv[i]);
					fileNameVCF.push_back(vcfTmp);
				}
			}

			//ped file name
			else if(option=="pedFile")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The ped file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					fileNamePed=argv[i];
				}
			}

			//location file
			else if(option=="l")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The location file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					locationFile=argv[i];
				}
			}

			//output file name
			else if(option=="output")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The output file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					outputName=argv[i];
				}
			}

			//varOnly
			else if(option=="v")
			{
				varOnly=true;
			}

			//method
			else if(option=="method")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Method hasn't been set. The default method (BN) will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					method=atoi(argv[i]);
				}
			}

			//mutation rate
			else if(option=="mRate")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Mutation rate hasn't been set. The default (1e-7) will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					mRate=atof(argv[i]);
				}
			}

			//genotype probability for New variant
			else if(option=="genoProbN")
			{
				vector<double> gpTmp(3);
				bool flag=true;
				for(int j=0;j<3;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.9985,0.001,0.0005) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						gpTmp[j]=atof(argv[i]);
					}
				}
				if(flag)
				{
					genoProbN=gpTmp;
				}
			}

			//genotype probability for Known variant
			else if(option=="genoProbK")
			{
				vector<double> gpTmp(3);
				bool flag=true;
				for(int j=0;j<3;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbK hasn't been set. The default (0.45,0.1,0.45) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						gpTmp[j]=atof(argv[i]);
					}
				}
				if(flag)
				{
					genoProbK=gpTmp;
				}
			}

			//genotype probability for New variant on X chromosome for male
			else if(option=="genoProbXN")
			{
				vector<double> gpTmp(3,0);
				bool flag=true;
				for(int j=0;j<2;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.999,0.001) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						if(j==0)
						{
							gpTmp[0]=atof(argv[i]);
						}
						else
						{
							gpTmp[2]=atof(argv[i]);
						}
					}
				}
				if(flag)
				{
					genoProbXN=gpTmp;
				}
			}

			//genotype probability for Known variant on X chromosome for male
			else if(option=="genoProbXK")
			{
				vector<double> gpTmp(3,0);
				bool flag=true;
				for(int j=0;j<2;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.5,0.5) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						if(j==0)
						{
							gpTmp[0]=atof(argv[i]);
						}
						else
						{
							gpTmp[2]=atof(argv[i]);
						}
					}
				}
				if(flag)
				{
					genoProbXK=gpTmp;
				}
			}

			else if(option=="o")
			{
				posOrder=true;
			}

			else if(option=="a")
			{
				allLine=true;
			}

			else if(option=="numBurnIn")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Number of burn in times hasn't been set. The default 1000 will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					numBurnIn=atoi(argv[i]);
				}
			}

			else if(option=="numRep")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Number of MCMC repeat times hasn't been set. The default 100000 will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					numRep=atoi(argv[i]);
				}
			}

			else if(option == "LRC")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cerr<< "Likelihood ratio criteria is not set. The default will be used."<<endl;
					i--;
					returnVal = 1;
					continue;
				}
				else
				{
					lrc=atof(argv[i]);
				}
			}

			else if(option == "d")
			{
				diffOnly = true;
			}

			else
			{
				cout<<"Cannot recognize option: \""<<option<<"\" in the command."<<endl;
				returnVal=1;
				continue;
			}
		}
		else
		{
			cout<<"Cannot recognize parameter: \""<<argv[i]<<"\" in the command."<<endl;
			returnVal=1;
		}
	}

	//check the input

	if(fileNameVCF.size()==0)
	{
		cout<<"The name of vcf file must be set. Please input the vcf file name."<<endl;
		return -1;
	}

	if(fileNamePed.size()==0)
	{
		cout<<"The name of ped file must be set. Please input the ped file name."<<endl;
		return -1;
	}

	if(outputName.size()==0)
	{
		cout<<"The name of output file must be set. Please input the output file name."<<endl;
		return -1;
	}

	if(method<1 || method>3)
	{
		cout<<"Method could only be 1 or 3. The default method (BN) will be used."<<endl;
		method=1;
		returnVal=1;
	}

	if(mRate<0 || mRate>0.5)
	{
		cout<<"Mutation rate is set out of range. The default (1e-7) will be used."<<endl;
		mRate=1e-7;
		returnVal=1;
	}

	if(varOnly==true && allLine==true)
	{
		cout<<"varOnly is setted, allLine is blocked."<<endl;
		allLine=false;
		returnVal=1;
	}

	if(numBurnIn<0)
	{
		if(numBurnIn != -999){
			cout<<"Number of burn in cannot be less than 0. The default 1000*n will be used."<<endl;
			numBurnIn=-1;
			returnVal=1;
		}
		else{
			numBurnIn = -1;
		}
	}

	if(numRep<=0)
	{
		if(numRep != -999){
			cout<<"Number of MCMC repeat times cannot be less than 1. The default 20000*n will be used."<<endl;
			numRep=-1;
			returnVal=1;
		}
		else{
			numRep = -1;
		}
	}

	if(lrc<0)
	{
		cerr<<"Likelihood ration criteria is not set correctly. The default will be used."<<endl;
		lrc=1;
		returnVal = 1;
	}

	vcfInfo.set_fileNameVCF(fileNameVCF);
	vcfInfo.set_fileNamePed(fileNamePed);
	vcfInfo.set_outputName(outputName);
	vcfInfo.set_locationFile(locationFile);
	vcfInfo.set_varOnly(varOnly);
	vcfInfo.set_method(method);
	vcfInfo.set_mRate(mRate);
	vcfInfo.set_genoProbN(genoProbN);
	vcfInfo.set_genoProbK(genoProbK);
	vcfInfo.set_genoProbXN(genoProbXN);
	vcfInfo.set_genoProbXK(genoProbXK);
	vcfInfo.set_posOrder(posOrder);
	vcfInfo.set_allLine(allLine);
	vcfInfo.set_diffOnly(diffOnly);
	vcfInfo.set_numBurnIn(numBurnIn);
	vcfInfo.set_numRep(numRep);
	vcfInfo.set_LCriteria(lrc);
	return returnVal;
}

inputLK::inputLK(const std::string & FileNameLK, const std::string & FileNamePed, const std::string & OutputName, int Method, double MRate)
{
	fileNameLK=FileNameLK;
	fileNamePed=FileNamePed;
	outputName=OutputName;
	method=Method;
	mRate=MRate;
}

inputLK::inputLK()
{
}

bool inputLK::set_fileNameLK(string FileNameLK)
{
	fileNameLK=FileNameLK;
	return true;
}

bool inputLK::set_fileNamePed(string FileNamePed)
{
	fileNamePed=FileNamePed;
	return true;
}

bool inputLK::set_genoProbK(const vector<double> & GenoProbK)
{
	genoProbK=GenoProbK;
	return true;
}

bool inputLK::set_genoProbN(const vector<double> & GenoProbN)
{
	genoProbN=GenoProbN;
	return true;
}

bool inputLK::set_genoProbXK(const vector<double> & GenoProbXK)
{
	genoProbXK=GenoProbXK;
	return true;
}

bool inputLK::set_genoProbXN(const vector<double> & GenoProbXN)
{
	genoProbXN=GenoProbXN;
	return true;
}

bool inputLK::set_mRate(double MRate)
{
	mRate=MRate;
	return true;
}

bool inputLK::set_method(int Method)
{
	method=Method;
	return true;
}

bool inputLK::set_outputName(string OutputName)
{
	outputName=OutputName;
	return true;
}

bool inputLK::set_lkType(int LKType)
{
	lkType=LKType;
	return true;
}

bool inputLK::set_numBurnIn(int NumBurnIn)
{
	numBurnIn=NumBurnIn;
	return true;
}

bool inputLK::set_numRep(int NumRep)
{
	numRep=NumRep;
	return true;
}

bool inputLK::set_LCriteria(double LC)
{
	LCriteria = LC;
	return true;
}

int checkInputLK(int argc, char* argv[], inputLK & lkInfo)
{
	string fileNameLK;
	string fileNamePed;
	string outputName;
	int method=1;
	int lkType=1;
	double mRate=1e-7;
	vector<double> genoProbN;
	vector<double> genoProbK;
	vector<double> genoProbXN;
	vector<double> genoProbXK;
	int numBurnIn=1000;
	int numRep=100000;
	double lrc = 1;

	int returnVal=0;

	for(int i=2;i<argc;i++)
	{
		if(argv[i][0]=='-')
		{
			string strOption(argv[i]);
			string option=strOption.substr(1);

			//lk
			if(option=="lkFile")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The likelihood file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					fileNameLK=argv[i];
				}
			}

			//ped file name
			else if(option=="pedFile")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The ped file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					fileNamePed=argv[i];
				}
			}

			//output file name
			else if(option=="output")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"The output file hasn't been set. Please check input."<<endl;
					return -1;
				}
				else
				{
					outputName=argv[i];
				}
			}

			//method
			else if(option=="method")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Method hasn't been set. The default method (BN) will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					method=atoi(argv[i]);
				}
			}

			//mutation rate
			else if(option=="mRate")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Mutation rate hasn't been set. The default (1e-7) will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					mRate=atof(argv[i]);
				}
			}

			//genotype probability for New variant
			else if(option=="genoProbN")
			{
				vector<double> gpTmp(3);
				bool flag=true;
				for(int j=0;j<3;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.9985,0.001,0.0005) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						gpTmp[j]=atof(argv[i]);
					}
				}
				if(flag)
				{
					genoProbN=gpTmp;
				}
			}

			//genotype probability for Known variant
			else if(option=="genoProbK")
			{
				vector<double> gpTmp(3);
				bool flag=true;
				for(int j=0;j<3;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbK hasn't been set. The default (0.45,0.1,0.45) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						gpTmp[j]=atof(argv[i]);
					}
				}
				if(flag)
				{
					genoProbK=gpTmp;
				}
			}

			//genotype probability for New variant on X chromosome for male
			else if(option=="genoProbXN")
			{
				vector<double> gpTmp(3,0);
				bool flag=true;
				for(int j=0;j<2;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.999,0.001) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						if(j==0)
						{
							gpTmp[0]=atof(argv[i]);
						}
						else
						{
							gpTmp[2]=atof(argv[i]);
						}
					}
				}
				if(flag)
				{
					genoProbXN=gpTmp;
				}
			}

			//genotype probability for Known variant on X chromosome for male
			else if(option=="genoProbXK")
			{
				vector<double> gpTmp(3,0);
				bool flag=true;
				for(int j=0;j<2;j++)
				{
					i++;
					if(i==argc || argv[i][0]=='-')
					{
						cout<<"genoProbN hasn't been set. The default (0.5,0.5) will be used."<<endl;
						i--;
						flag=false;
						returnVal=1;
						break;
					}
					else
					{
						if(j==0)
						{
							gpTmp[0]=atof(argv[i]);
						}
						else
						{
							gpTmp[2]=atof(argv[i]);
						}
					}
				}
				if(flag)
				{
					genoProbXK=gpTmp;
				}
			}

			else if(option=="lkType")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Likelihood type hasn't been set. The default normal(n) will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					string lkTypeStr(argv[i]);
					if(lkTypeStr=="n")
					{
						lkType=1;
					}
					else if(lkTypeStr=="log10")
					{
						lkType=2;
					}
					else if(lkTypeStr=="ln")
					{
						lkType=3;
					}
					else if(lkTypeStr=="PS")
					{
						lkType=4;
					}
					else
					{
						cout<<"Cannot recognize the likelihood type: "<<lkTypeStr<<". The default normal (n) will be used."<<endl;
						lkType=1;
						returnVal=1;
					}
				}
			}

			else if(option=="numBurnIn")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Number of burn in times hasn't been set. The default 1000 will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					numBurnIn=atoi(argv[i]);
				}
			}

			else if(option=="numRep")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cout<<"Number of MCMC repeat times hasn't been set. The default 100000 will be used."<<endl;
					i--;
					returnVal=1;
					continue;
				}
				else
				{
					numRep=atoi(argv[i]);
				}
			}

			else if(option == "LRC")
			{
				i++;
				if(i==argc || argv[i][0]=='-')
				{
					cerr<< "Likelihood ratio criteria is not set. The default will be used."<<endl;
					i--;
					returnVal = 1;
					continue;
				}
				else
				{
					lrc=atof(argv[i]);
				}
			}

			else
			{
				cout<<"Cannot recognize option: \""<<option<<"\" in the command."<<endl;
				returnVal=1;
				continue;
			}
		}
		else
		{
			cout<<"Cannot recognize parameter: \""<<argv[i]<<"\" in the command."<<endl;
			returnVal=1;
		}
	}

	//check the input

	if(fileNameLK.size()==0)
	{
		cout<<"The name of likelihood file must be set. Please input the likelihood file name."<<endl;
		return -1;
	}

	if(fileNamePed.size()==0)
	{
		cout<<"The name of ped file must be set. Please input the ped file name."<<endl;
		return -1;
	}

	if(outputName.size()==0)
	{
		cout<<"The name of output file must be set. Please input the output file name."<<endl;
		return -1;
	}

	if(method<1 || method>3)
	{
		cout<<"Method could only be 1 or 3. The default method (BN) will be used."<<endl;
		method=1;
		returnVal=1;
	}

	if(mRate<0 || mRate>0.5)
	{
		cout<<"Mutation rate is set out of range. The default (1e-7) will be used."<<endl;
		mRate=1e-7;
		returnVal=1;
	}

	if(numBurnIn<0)
	{
		cout<<"Number of burn in cannot be less than 0. The default 1000 will be used."<<endl;
		numBurnIn=1000;
		returnVal=1;
	}

	if(numRep<=0)
	{
		cout<<"Number of MCMC repeat times cannot be less than 1. The default 100000 will be used."<<endl;
		numRep=100000;
		returnVal=1;
	}

	if(lrc<0)
	{
		cerr<<"Likelihood ration criteria is not set correctly. The default will be used."<<endl;
		lrc=1;
		returnVal = 1;
	}

	lkInfo.set_fileNameLK(fileNameLK);
	lkInfo.set_fileNamePed(fileNamePed);
	lkInfo.set_genoProbK(genoProbK);
	lkInfo.set_genoProbN(genoProbN);
	lkInfo.set_genoProbXK(genoProbXK);
	lkInfo.set_genoProbXN(genoProbXN);
	lkInfo.set_mRate(mRate);
	lkInfo.set_method(method);
	lkInfo.set_lkType(lkType);
	lkInfo.set_outputName(outputName);
	lkInfo.set_numBurnIn(numBurnIn);
	lkInfo.set_numRep(numRep);
	lkInfo.set_LCriteria(lrc);

	return returnVal;
}
