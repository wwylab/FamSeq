/*
 * FamSeqPro.cpp
 * main function of FamSeq
 *
 *  Created on: Mar 7, 2012
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>

#include "family.h"
#include "normal.h"
#include "file.h"
#include "checkInput.h"

using namespace std;

int main(int argc, char* argv[])
{
	if(argc==1)
	{
		cout<<endl;
		cout<<"Program: FamSeq (Sequence calling using pedigree information)"<<endl;
		cout<<"Version: 1.0.2"<<endl<<endl;
		cout<<"Usage:\tFamSeq <input type> [options]"<<endl<<endl;
		cout<<"Input type: \tvcf\t\tinput vcf file"<<endl;
		cout<<"\t\tLK\t\tinput likelihood file"<<endl;
		cout<<endl;
		cout<<"Type FamSeq -h for help." <<endl;
		cout<<endl;

		return -1;
	}

	//vcf
	if((strcmp(argv[1],"vcf"))==0)
	{
		//show help message for vcf file input
		if(argc==2)
		{
			cout<<endl;
			cout<<"Usage:\tFamSeq vcf [options]"<<endl<<endl;
			cout<<"Call variants when the input data is in a vcf file."<<endl;
			cout<<endl;
			cout<<"Type FamSeq -h for help." <<endl;
			return -1;
		}

		inputVCF vcfInfo;
		int checkReturn=checkInputVCF(argc,argv,vcfInfo);
		if(checkReturn<0)
		{
			return -1;
		}

		if(checkReturn>0)
		{
			cout<<"There are some improper parameters in the command line. Some parameters are set to default."<<endl;
		}

		family fam;
		if(!setFam(vcfInfo,fam))
		{
			cout<<"Cannot set family."<<endl;
			return -1;
		}

		if(vcfInfo.get_fileNameVCF().size()==1)
		{
			//callGenoMVCF((vcfInfo.get_fileNameVCF())[0],vcfInfo.get_outputName(),fam,vcfInfo.get_varOnly(),vcfInfo.get_method());
			callGenoMVCF(vcfInfo,fam);
		}
		else
		{
			//
		}

	}
	//likelihood
	else if((strcmp(argv[1],"LK"))==0)
	{
		if(argc==2)
		{
			cout<<endl;
			cout<<"Usage:\tFamSeq LK [options]"<<endl<<endl;
			cout<<"Call variants when the input data is in a vcf file."<<endl;
			cout<<endl;
			cout<<"Type FamSeq -h for help." <<endl;
			return -1;
		}

		inputLK lkInfo;
		int checkReturn=checkInputLK(argc,argv,lkInfo);
		if(checkReturn<0)
		{
			return -1;
		}

		if(checkReturn>0)
		{
			cout<<"There are some improper parameters in the command line. Some parameters are set to default."<<endl;
		}

		family fam;
		if(!setFam(lkInfo,fam))
		{
			cout<<"Cannot set Family."<<endl;
			return  -1;
		}

		callGenoLK(lkInfo,fam);
	}
	else if((strcmp(argv[1],"-h"))==0){
		cout<<"FamSeq: Version: 1.0.2"<<endl;
		cout<<"Usage:\tFamSeq <input type> [options]"<<endl<<endl;
		cout<<"FamSeq accepts two kinds of input files: vcf file and likelihood only format file. ";
		cout<<"If the input is vcf file, type 'FamSeq vcf [options]' in the command line. ";
		cout<<"Type 'FamSeq LK [options]' if the input is likelihood only format. The user can only use only of them.";
		cout<<endl<<endl;
		cout<<"Options:"<<endl<<endl;
		cout<<"-vcfFile\tThe name of input vcf file."<<endl<<endl;
		cout<<"-lkFile\t\tThe name of input likelihood only format file."<<endl<<endl;
		cout<<"-lkType\tThe likelihood type stored in the likelihood only format file. n:normal(default); log10: log10 scaled; ln: ln scaled; PS: phred scaled."<<endl<<endl;
		cout<<"-pedFile\tThe name of the file storing the pedigree information."<<endl<<endl;
		cout<<"-output\t\tThe name of output file"<<endl<<endl;
		cout<<"-method\t\tChoose the method used in variant calling. 1(default): Bayesian network; 2: Elston-Stewart algorithm; 3: MCMC."<<endl<<endl;
		cout<<"-mRate\t\tMutation rate. The default value is 1e-7"<<endl<<endl;
		cout<<"-v\t\tOnly record the position at which the genotype is not RR in the output file. (R: reference allele, A: alternative allele)."<<endl<<endl;
		cout<<"-a\t\tRecord all the position in the output file."<<endl<<endl;
		cout<<"-genoProbN\tGenotype probability of three kinds of genotype for autosome in population (Pr(G)) when the variant is not in dbSNP. The default value is:  0.9985, 0.001 and 0.0005. The dbSNP position should be provided in column ID in input vcf file. "<<endl<<endl;
		cout<<"-genoProbK\tGenotype probability of three kinds of genotype for autosome in population (Pr(G)) when the variant is in dbSNP. The default value is: 0.45, 0.1 and 0.45."<<endl<<endl;
		cout<<"-genoProbXN\tGenotype probability of two kinds of genotype for chromosome X for male in population (Pr(G)) when the variant is not in dbSNP. The default value is: 0.999 and 0.001."<<endl<<endl;
		cout<<"-genoProbXK\tGenotype probability of two kinds of genotype for chromosome X for male in population (Pr(G)) when the variant is in dbSNP. The default value is: 0.5 and 0.5."<<endl<<endl;
		cout<<"-numBurnIn\tNumber of burn in when the user chooses the MCMC method. The default value is 1,000n, where n is the number of individuals in the pedigree."<<endl<<endl;
		cout<<"-numRep\t\tNumber of iteration times when the user chooses MCMC method. The default value is 20,000n. "<<endl<<endl;
	}
	else
	{
		cout<<"Cannot recognize the input type: \""<<argv[1]<<"\"."<<endl;
		cout<<"The input type can only be vcf or LK"<<endl;
		cout<<endl;
		cout<<"Type FamSeq -h for help." <<endl;
		return -1;
	}
	return 0;
}
