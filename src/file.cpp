/*
 * file.cpp
 *
 *  Created on: Mar 16, 2012
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <climits>
#include <algorithm>

#include "file.h"
#include "normal.h"

using namespace std;

bool readPed(string fileName, vector<individual> & mem)
{
	ifstream fin;
	fin.open(fileName.c_str());
	if(!fin.is_open())
	{
		cout<<"Cannot open "<<fileName<<endl;
		return false;
	}

	string headings;
	getline(fin,headings);

	mem.clear();
	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);
		if(fline.size()<2)
		{
			break;
		}
		int id, mid, fid, gender;
		string sName;
		istringstream instr(fline);
		instr>>id;
		instr>>mid;
		instr>>fid;
		instr>>gender;
		instr>>sName;


		individual indTmp(id,gender,0,0,mid,fid,sName);
		mem.push_back(indTmp);
	}
	fin.close();
	fin.clear();
	return true;
}

bool readVCF(std::string fileName, std::string & header, std::vector<std::vector<int> > & pos, std::vector<std::vector<string> > & id, std::vector<std::vector<char> > & ref, std::vector<std::vector<char> > & alt, std::vector<std::vector<int> > & dp, std::vector<std::vector<std::vector<string> > > & others, std::vector<std::vector<std::string> > & genotype, std::vector<std::vector<int> > & pl1, std::vector<std::vector<int> > & pl2, std::vector<std::vector<int> > & pl3)
{
	ifstream fin(fileName.c_str());
	if(!fin.is_open())
	{
		cout<<"Cannot open "<<fileName<<endl;
		return false;
	}

	pos.clear();
	id.clear();
	ref.clear();
	alt.clear();
	dp.clear();
	others.clear();
	genotype.clear();
	pl1.clear();
	pl2.clear();
	pl3.clear();

	pos.resize(24);
	id.resize(24);
	ref.resize(24);
	alt.resize(24);
	dp.resize(24);
	others.resize(24);
	for(int i=0;i<24;i++)
	{
		others[i].resize(3);
	}
	genotype.resize(24);
	pl1.resize(24);
	pl2.resize(24);
	pl3.resize(24);

	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);
	}
	return true;
}


bool callGenoMVCF(const inputVCF & vcfInfo, family fam)
//bool callGenoMVCF(std::string vcfFileName, std::string outFileName, family fam, bool VarOnly, int method)
{
	//int NumRep = 2000000;
	//int NumRep = 10000;

	string vcfFileName=(vcfInfo.get_fileNameVCF())[0];
	string outFileName=vcfInfo.get_outputName();
	string locationFile=vcfInfo.get_locationFile();
	bool VarOnly=vcfInfo.get_varOnly();
	int method=vcfInfo.get_method();
	bool allLine=vcfInfo.get_allLine();
	vector<double> gProbN=fam.get_genoProbN();
	vector<double> gProbK=fam.get_genoProbK();
	vector<double> gProbXN=fam.get_genoProbXN();
	vector<double> gProbXK=fam.get_genoProbXK();
	bool diffOnly = vcfInfo.get_diffOnly();
	if(diffOnly){
		allLine = false;
		VarOnly = true;
	}

	ifstream fin(vcfFileName.c_str());
	if(!fin.is_open())
	{
		cout<<"Cannot open "<<vcfFileName<<endl;
		return false;
	}

	string title="";

	ofstream fout(outFileName.c_str());

	bool flagFSPInfo=true;
	bool flagFSPTag=true;
	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);

		if(fline[0]=='#')
		{
			if(title.size()>2)
			{
				fout<<title<<endl;
				if(title.substr(2,6)=="FORMAT" && fline.substr(2,4)=="INFO" && flagFSPTag)
				{
					fout<<"##FORMAT=<ID=GPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calculated by individual-based Method\">"<<endl;
					fout<<"##FORMAT=<ID=FPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calculated by FamSeqPro\">"<<endl;
					fout<<"##FORMAT=<ID=FGT,Number=1,Type=String,Description=\"Genotype called by FamSeqPro\">"<<endl;
					flagFSPTag=false;
				}

				//Show some information in FSP
				if(fline.substr(2,6)=="contig" && flagFSPInfo)
				{
					fout<<"##FS mutation rate="<<fam.get_mRate()<<" "<<endl;
					fout<<"##FS genotype frequency in pupulation (Rare): "<<gProbN[0]<<":"<<gProbN[1]<<":"<<gProbN[2]<<endl;
					fout<<"##FS genotype frequency in population (Common): "<<gProbK[0]<<":"<<gProbK[1]<<":"<<gProbK[2]<<endl;
					fout<<"##FS genotype frequency for chromosome X of male in population (Rare): "<<gProbXN[0]<<":"<<gProbXN[1]<<":"<<gProbXN[2]<<endl;
					fout<<"##FS genotype frequency for chromosome X of male in population (Common): "<<gProbXK[0]<<":"<<gProbXK[1]<<":"<<gProbXK[2]<<endl;
					flagFSPInfo=false;
				}
			}
			title=fline;
			continue;
		}
		else
		{
			if(flagFSPTag)
			{
				fout<<"##FORMAT=<ID=GPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calbulated by Single Method\">"<<endl;
				fout<<"##FORMAT=<ID=FPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calculated by FamSeqPro\">"<<endl;
				fout<<"##FORMAT=<ID=FGT,Number=1,Type=String,Description=\"Genotype called by FamSeqPro\">"<<endl;
				flagFSPTag=false;
			}

			if(flagFSPInfo)
			{
				fout<<"##FS mutation rate="<<fam.get_mRate()<<" "<<endl;
				fout<<"##FS genotype frequency in pupulation (Rare): "<<gProbN[0]<<":"<<gProbN[1]<<":"<<gProbN[2]<<endl;
				fout<<"##FS genotype frequency in population (Common): "<<gProbK[0]<<":"<<gProbK[1]<<":"<<gProbK[2]<<endl;
				fout<<"##FS genotype frequency for chromosome X of male in population (Rare): "<<gProbXN[0]<<":"<<gProbXN[1]<<":"<<gProbXN[2]<<endl;
				fout<<"##FS genotype frequency for chromosome X of male in population (Common): "<<gProbXK[0]<<":"<<gProbXK[1]<<":"<<gProbXK[2]<<endl;
				flagFSPInfo=false;
			}
			break;
		}
	}
	fin.close();
	fin.clear();

	vector<string> vsTitle=split(title,"\t");
	for(int i=0;i<9;i++)
	{
		fout<<vsTitle[i]<<'\t';
	}
	//mapV2P:map index from vcf file to ped file
	//mapP2V:map index from ped file to vcf file
	vector<int> mapV2P(vsTitle.size()-9,-1), mapP2V(fam.get_numInd(),-1);
	vector<individual> mem=fam.get_member();
	for(size_t i=0;i<mapV2P.size();i++)
	{
		for(size_t j=0;j<mapP2V.size();j++)
		{
			if(vsTitle[9+i]==mem[j].get_sName())
			{
				mapV2P[i]=j;
				mapP2V[j]=i;
				break;
			}
		}
	}

	fam.set_mapP2V(mapP2V);
	fam.set_mapV2P(mapV2P);

	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			fout<<vsTitle[9+i]<<'\t';
		}
	}
	fout<<endl;


	bool flagLocation=false;
	vector<vector<int> > location;
	//1-22,X,Y,MT
	if(locationFile.size()>0)
	{
		flagLocation=true;
		ifstream finL(locationFile.c_str());
		if(!finL.is_open())
		{
			cout<<"Cannot open "<<locationFile<<endl;
			return false;
		}
		for(int i=0;i<25;i++)
		{
			vector<int> tmp;
			location.push_back(tmp);
		}

		while(!finL.eof())
		{
			string fline;
			getline(finL,fline);
			if(fline.size()<2)
			{
				break;
			}
			vector<string> vsLine=split(fline,"\t");
			int chrTmp=0;
			if(vsLine[0]=="X")
			{
				chrTmp=23;
			}
			else if(vsLine[0]=="Y")
			{
				chrTmp=24;
			}
			else if(vsLine[0]=="MT")
			{
				chrTmp=25;
			}
			else
			{
				chrTmp=atoi(vsLine[0].c_str());
			}
			if(chrTmp==0)
			{
				continue;
			}
			int posTmp=0;
			posTmp=atoi(vsLine[1].c_str());
			if(posTmp==0)
			{
				continue;
			}
			location[chrTmp-1].push_back(posTmp);
		}

		for(int i=0;i<25;i++)
		{
			sort(location[i].begin(),location[i].end());
		}
		finL.close();
		finL.clear();
	}

	fin.open(vcfFileName.c_str());
	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);
		if(fline.size()<2)
		{
			break;
		}

		if(fline[0]=='#')
		{
			continue;
		}

		vector<string> vsLine=split(fline,"\t");

		//selected the position first
		if(flagLocation)
		{
			int chrTmp=0;
			if(vsLine[0]=="X" || vsLine[0]=="chrX")
			{
				chrTmp=23;
			}
			else if(vsLine[0]=="Y" || vsLine[0]=="chrY")
			{
				chrTmp=24;
			}
			else if(vsLine[0]=="MT")
			{
				chrTmp=25;
			}
			else
			{
				if(vsLine[0].substr(0,3)=="chr")
				{
					chrTmp=atoi(vsLine[0].substr(3).c_str());
				}
				else
				{
					chrTmp=atoi(vsLine[0].c_str());
				}
			}
			if(chrTmp==0)
			{
				continue;
			}
			int posTmp=0;
			posTmp=atoi(vsLine[1].c_str());
			if(posTmp==0)
			{
				continue;
			}

			int indexFind=binSearch(location[chrTmp-1],posTmp);
			if(indexFind<0)
			{
				continue;
			}
		}

		if(vsLine[3]=="." || vsLine[3]=="-")
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		if(vsLine[3].size()!=1 || vsLine[4].size()!=1)
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		if(VarOnly && (vsLine[4]=="." || vsLine[4]=="-"))
		{
			continue;
		}

		if(vsLine[0]=="Y" || vsLine[0]=="chrY")
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		if(vsLine[0]=="MT")
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		int chrTmp=0;
		if(vsLine[0].substr(0,3)=="chr"){
			chrTmp=atoi(vsLine[0].substr(3).c_str());
		}
		else{
			chrTmp=atoi(vsLine[0].c_str());
		}

		if(!((0<chrTmp && chrTmp<23) || vsLine[0]=="X" || vsLine[0]=="chrX" || vsLine[0]=="CHRX"))
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		//
		bool KnowSNP=false;
		if(vsLine[2]!=".")
		{
			KnowSNP=true;
		}

		int chrType=0;
		if(vsLine[0]=="X" || vsLine[0]=="chrX" || vsLine[0]=="CHRX")
		{
			chrType=1;
		}

		vector<string> vsInfo=split(vsLine[8],":");
		int numMiss=0;
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				if(vsLine[9+i].size()<5)
				{
					numMiss++;
				}
			}
		}
		if(numMiss==int(fam.get_realNumInd()))
		{
			if(allLine)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			continue;
		}

		int indDP=-1,indPL=-1,indGT=-1;
		for(size_t i=0;i<vsInfo.size();i++)
		{
			if(vsInfo[i]=="DP")
			{
				indDP=i;
			}

			if(vsInfo[i]=="PL" || vsInfo[i] == "GL")
			{
				indPL=i;
			}

			if(vsInfo[i]=="GT")
			{
				indGT=i;
			}
		}

		if(numMiss==0)
		{
			if(indPL<0)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			else
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<vsLine[8];
				fout<<":GPP:FPP:FGT\t";

				dMatrix<double> LKTmp(fam.get_numInd(),3,1);
				//vector<string> PLTmp;
				//vector<string> GTTmp;

				//int diff = 0;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						vector<string> vsGT=split(vsLine[9+i],":");
						if(vsGT.size()!=vsInfo.size())
						{
							//cout<<vsInfo.size()<<'\t'<<vsGT.size()<<'\t';
							//diff = vsInfo.size() - vsGT.size();
							//cout<<diff<<endl;
							continue;
						}
						//fout<<vsGT[indDP]<<'\t';
						//GTTmp.push_back(vsGT[indGT]);
						vector<string> vsPL=split(vsGT[indPL],",");
						for(int j=0;j<3;j++)
						{
							//PLTmp.push_back(vsPL[j]);
							double lkTmp=atof(vsPL[j].c_str());
							lkTmp=pow(10.0,-fabs(lkTmp)/10.0);
							LKTmp(mapV2P[i],j)=lkTmp;
						}
					}
				}

				fam.set_LK(LKTmp);

				//cout<<fline<<endl;
				//cout<<LKTmp<<endl;

				//different method
				//BN
				if(method==1)
				{
//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbBN(KnowSNP,chrType);
//					}
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else if(method==2)
				{
//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbPeeling(KnowSNP,chrType);
//					}
					if(!fam.calPostProbPeeling(KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else if(method==3)
				{
					int numBurnIn, numRep;
					if(vcfInfo.get_numBurinIn()<0){
						numBurnIn = 1000*int(fam.get_realNumInd());
					}
					else{
						numBurnIn = vcfInfo.get_numBurinIn();
					}
					if(vcfInfo.get_numRep()<0){
						numRep = 20000*int(fam.get_realNumInd());
					}
					else{
						numRep = vcfInfo.get_numRep();
					}

//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbMCMC(numBurnIn,numRep,KnowSNP,chrType);
//					}
					if(!fam.calPostProbMCMC(numBurnIn,numRep,KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else
				{
				}

				dMatrix<double> postProb=fam.get_postProb();
				vector<int> postRlt=fam.get_postRlt();
				dMatrix<double> postProbSingle=fam.get_postProbSingle();

				int indTmp=0;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<":";
						//cout<<'\t'<<diff<<endl;
						/*
						if(diff>0)
						{
							cout<<vsLine[1]<<endl;
							for(int j=0;j<diff;j++)
							{
								fout<<"NA:";
							}
						}
						*/

						if(-10*log10(postProbSingle(indTmp,0))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,0)))<<',';
						}
						if(-10*log10(postProbSingle(indTmp,1))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,1)))<<',';
						}
						if(-10*log10(postProbSingle(indTmp,2))==numeric_limits<double>::infinity())
						{
							fout<<99999<<':';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,2)))<<':';
						}
						if(-10*log10(postProb(indTmp,0))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,0)))<<',';
						}
						if(-10*log10(postProb(indTmp,1))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,1)))<<',';
						}
						if(-10*log10(postProb(indTmp,2))==numeric_limits<double>::infinity())
						{
							fout<<99999<<':';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,2)))<<':';
						}
						if(postRlt[indTmp]==0)
						{
							fout<<"0/0\t";
						}
						else if(postRlt[indTmp]==1)
						{
							fout<<"0/1\t";
						}
						else
						{
							fout<<"1/1\t";
						}
						indTmp++;
					}
				}
				fout<<endl;
			}
		}
		else
		{
			if(indPL<0)
			{
				for(int i=0;i<9;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[9+i]<<'\t';
					}
				}
				fout<<endl;
			}
			else
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<vsLine[8];
				fout<<":GPP:FPP:FGT\t";

				dMatrix<double> LKTmp(fam.get_numInd(),3,1);
				//vector<string> PLTmp;
				//vector<string> GTTmp;
				//int numFormat=0;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						if(vsLine[9+i].size()<5)
						{
							//fout<<"0\t";
							//GTTmp.push_back("NA");
							for(int j=0;j<3;j++)
							{
								//PLTmp.push_back("1");
								LKTmp(mapV2P[i],j)=1;
							}
						}
						else
						{
							vector<string> vsGT=split(vsLine[9+i],":");
							//numFormat=int(vsGT.size());
							if(vsGT.size()!=vsInfo.size()){
								continue;
							}
							//fout<<vsGT[indDP]<<'\t';
							//GTTmp.push_back(vsGT[indGT]);
							vector<string> vsPL=split(vsGT[indPL],",");
							for(int j=0;j<3;j++)
							{
								//PLTmp.push_back(vsPL[j]);
								double lkTmp=atof(vsPL[j].c_str());
								lkTmp=pow(10.0,-fabs(lkTmp)/10.0);
								LKTmp(mapV2P[i],j)=lkTmp;
							}
						}
					}
				}

				fam.set_LK(LKTmp);

				//cout<<fline<<endl;
				//cout<<LKTmp<<endl;

				//different method
				//BN
				if(method==1)
				{
//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbBN(KnowSNP,chrType);
//					}
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else if(method==2)
				{
//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbPeeling(KnowSNP,chrType);
//					}
					if(!fam.calPostProbPeeling(KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else if(method==3)
				{
					int numBurnIn, numRep;
					if(vcfInfo.get_numBurinIn()<0){
						numBurnIn = 1000*int(fam.get_realNumInd());
					}
					else{
						numBurnIn = vcfInfo.get_numBurinIn();
					}
					if(vcfInfo.get_numRep()<0){
						numRep = 20000*int(fam.get_realNumInd());
					}
					else{
						numRep = vcfInfo.get_numRep();
					}

//					for(int ii = 0;ii<NumRep;ii++){
//						fam.calPostProbMCMC(numBurnIn,numRep,KnowSNP,chrType);
//					}
					if(!fam.calPostProbMCMC(numBurnIn,numRep,KnowSNP,chrType))
					{
						cout<<"Warning: this variant hasn't been calculated: "<<endl;
						cout<<fline<<endl;
						for(size_t i=0;i<mapV2P.size();i++)
						{
							if(mapV2P[i]>=0)
							{
								fout<<vsLine[9+i]<<":NA:NA:NA\t";
							}
						}
						fout<<endl;
						continue;
					}
				}
				else
				{
				}

				dMatrix<double> postProb=fam.get_postProb();
				vector<int> postRlt=fam.get_postRlt();
				dMatrix<double> postProbSingle=fam.get_postProbSingle();

				int indTmp=0;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						if(vsLine[9+i].size()<5)
						{
							for(size_t j=0;j<vsInfo.size();j++)
							{
								fout<<"NA:";
							}
						}
						else
						{
							fout<<vsLine[9+i]<<":";
						}
						if(-10*log10(postProbSingle(indTmp,0))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,0)))<<',';
						}
						if(-10*log10(postProbSingle(indTmp,1))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,1)))<<',';
						}
						if(-10*log10(postProbSingle(indTmp,2))==numeric_limits<double>::infinity())
						{
							fout<<99999<<':';
						}
						else
						{
							fout<<fabs(-10*log10(postProbSingle(indTmp,2)))<<':';
						}
						if(-10*log10(postProb(indTmp,0))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,0)))<<',';
						}
						if(-10*log10(postProb(indTmp,1))==numeric_limits<double>::infinity())
						{
							fout<<99999<<',';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,1)))<<',';
						}
						if(-10*log10(postProb(indTmp,2))==numeric_limits<double>::infinity())
						{
							fout<<99999<<':';
						}
						else
						{
							fout<<fabs(-10*log10(postProb(indTmp,2)))<<':';
						}
						if(postRlt[indTmp]==0)
						{
							fout<<"0/0\t";
						}
						else if(postRlt[indTmp]==1)
						{
							fout<<"0/1\t";
						}
						else
						{
							fout<<"1/1\t";
						}
						indTmp++;
					}
				}
				fout<<endl;
			}
		}
	}


	/*
	string vcfFileName=(vcfInfo.get_fileNameVCF())[0];
	string outFileName=vcfInfo.get_outputName();
	bool VarOnly=vcfInfo.get_varOnly();
	int method=vcfInfo.get_method();
	bool allLine=vcfInfo.get_allLine();

	ifstream fin(vcfFileName.c_str());
	if(!fin.is_open())
	{
		cout<<"Cannot open "<<vcfFileName<<endl;
		return false;
	}

	string title;
	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);

		if(fline[0]=='#')
		{
			title=fline;
			continue;
		}
		else
		{
			break;
		}
	}
	fin.close();
	fin.clear();

	ofstream fout(outFileName.c_str());
	vector<string> vsTitle=split(title,"\t");
	for(int i=0;i<8;i++)
	{
		fout<<vsTitle[i]<<'\t';
	}

	vector<int> mapV2P(vsTitle.size()-9,-1), mapP2V(fam.get_numInd(),-1);
	vector<individual> mem=fam.get_member();
	for(size_t i=0;i<mapV2P.size();i++)
	{
		for(size_t j=0;j<mapP2V.size();j++)
		{
			if(vsTitle[9+i]==mem[j].get_sName())
			{
				mapV2P[i]=j;
				mapP2V[j]=i;
				break;
			}
		}
	}

	fam.set_mapP2V(mapP2V);
	fam.set_mapV2P(mapV2P);

	//
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			fout<<"DP_"<<vsTitle[9+i]<<'\t';
		}
	}

	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			for(int j=0;j<3;j++)
			{
				fout<<"PL_"<<vsTitle[9+i]<<'\t';
			}
		}
	}

	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			for(int j=0;j<3;j++)
			{
				fout<<"PP_"<<vsTitle[9+i]<<'\t';
			}
		}
	}

	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			fout<<"GT_S_"<<vsTitle[9+i]<<'\t';
		}
	}

	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			fout<<"GT_F_"<<vsTitle[9+i]<<'\t';
		}
	}
	fout<<endl;

	fin.open(vcfFileName.c_str());
	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);
		if(fline.size()<2)
		{
			break;
		}

		if(fline[0]=='#')
		{
			continue;
		}

		vector<string> vsLine=split(fline,"\t");

		if(vsLine[3].size()!=1 || vsLine[4].size()!=1)
		{
			if(allLine)
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<endl;
			}
			continue;
		}

		if(VarOnly && vsLine[4]==".")
		{
			continue;
		}

		if(vsLine[0]=="Y")
		{
			if(allLine)
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<endl;
			}
			continue;
		}

		if(vsLine[0].size()>2)
		{
			if(allLine)
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<endl;
			}
			continue;
		}

		bool KnowSNP=false;
		if(vsLine[2]!=".")
		{
			KnowSNP=true;
		}

		int chrType=0;
		if(vsLine[0]=="X")
		{
			chrType=1;
		}

		vector<string> vsInfo=split(vsLine[8],":");

		int numMiss=0;
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				if(vsLine[9+i].size()<5)
				{
					numMiss++;
				}
			}
		}
		if(numMiss==int(fam.get_realNumInd()))
		{
			if(allLine)
			{
				for(int i=0;i<8;i++)
				{
					fout<<vsLine[i]<<'\t';
				}
				fout<<endl;
			}
			continue;
		}

		int indDP=-1,indPL=-1,indGT=-1;
		for(size_t i=0;i<vsInfo.size();i++)
		{
			if(vsInfo[i]=="DP")
			{
				indDP=i;
			}

			if(vsInfo[i]=="PL")
			{
				indPL=i;
			}

			if(vsInfo[i]=="GT")
			{
				indGT=i;
			}
		}

		if(numMiss==0)
		{
			for(int i=0;i<8;i++)
			{
				fout<<vsLine[i]<<'\t';
			}

			if(indPL<0)
			{
				//DP
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						vector<string> vsGT=split(vsLine[9+i],":");
						if(indDP<int(vsGT.size()))
						{
							fout<<vsGT[indDP]<<'\t';
						}
						else
						{
							fout<<"NA\t";
						}
					}
				}
				//PL
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						for(int j=0;j<3;j++)
						{
							fout<<"NA\t";
						}
					}
				}
				//Posterior probability
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						for(int j=0;j<3;j++)
						{
							fout<<"NA\t";
						}
					}
				}
				//GT
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t'<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
				}
				fout<<endl;
			}
			else
			{
				dMatrix<double> LKTmp(fam.get_numInd(),3,1);
				vector<string> PLTmp;
				vector<string> GTTmp;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						vector<string> vsGT=split(vsLine[9+i],":");
						fout<<vsGT[indDP]<<'\t';
						GTTmp.push_back(vsGT[indGT]);
						vector<string> vsPL=split(vsGT[indPL],",");
						for(int j=0;j<3;j++)
						{
							PLTmp.push_back(vsPL[j]);
							double lkTmp=atoi(vsPL[j].c_str());
							lkTmp=pow(10.0,-lkTmp/10.0);
							LKTmp(mapV2P[i],j)=lkTmp;
						}
					}
				}

				for(size_t i=0;i<PLTmp.size();i++)
				{
					fout<<PLTmp[i]<<'\t';
				}

				fam.set_LK(LKTmp);

				//different method
				//BN
				if(method==1)
				{
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						//cout<<"This variant hasn't bee calculated: "<<endl;
						//cout<<fline<<endl;
						for(unsigned int i=0;i<fam.get_realNumInd();i++)
						{
							for(int j=0;j<3;j++)
							{
								fout<<"NA\t";
							}
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<GTTmp[i]<<'\t';
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<"NA\t";
						}
						fout<<endl;
						continue;
					}
				}

				//peeling
				if(method==2)
				{
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						cout<<"This variant hasn't bee calculated: "<<endl;
						cout<<fline<<endl;
						for(unsigned int i=0;i<fam.get_realNumInd();i++)
						{
							for(int j=0;j<3;j++)
							{
								fout<<"NA\t";
							}
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<GTTmp[i]<<'\t';
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<"NA\t";
						}
						fout<<endl;
						continue;
					}
				}
				dMatrix<double> postProb=fam.get_postProb();
				for(int i=0;i<postProb.get_row();i++)
				{
					for(int j=0;j<postProb.get_column();j++)
					{
						fout<<postProb(i,j)<<'\t';
					}
				}

				for(size_t i=0;i<GTTmp.size();i++)
				{
					if(GTTmp[i]=="0/0")
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
					else if(GTTmp[i]=="0/1")
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else if(GTTmp[i]=="1/1")
					{
						fout<<vsLine[4][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else
					{
						fout<<"NA\t";
					}
				}

				vector<int> postRlt=fam.get_postRlt();
				for(size_t i=0;i<postRlt.size();i++)
				{
					if(postRlt[i]==0)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
					else if(postRlt[i]==1)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else if(postRlt[i]==2)
					{
						fout<<vsLine[4][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else
					{
						fout<<"NA\t";
					}
				}
				fout<<endl;
			}
		}
		else
		{
			for(int i=0;i<8;i++)
			{
				fout<<vsLine[i]<<'\t';
			}

			if(indPL<0)
			{
				//DP
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						vector<string> vsGT=split(vsLine[9+i],":");
						if(indDP<int(vsGT.size()))
						{
							fout<<vsGT[indDP]<<'\t';
						}
						else
						{
							fout<<"NA\t";
						}
					}
				}
				//PL
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						for(int j=0;j<3;j++)
						{
							fout<<"NA\t";
						}
					}
				}
				//Posterior probability
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						for(int j=0;j<3;j++)
						{
							fout<<"NA\t";
						}
					}
				}
				//GT
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t'<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
				}
				fout<<endl;
			}
			else
			{
				dMatrix<double> LKTmp(fam.get_numInd(),3,1);
				vector<string> PLTmp;
				vector<string> GTTmp;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						if(vsLine[9+i].size()<5)
						{
							fout<<"0\t";
							GTTmp.push_back("NA");
							for(int j=0;j<3;j++)
							{
								PLTmp.push_back("1");
								LKTmp(mapV2P[i],j)=1;
							}
						}
						else
						{
							vector<string> vsGT=split(vsLine[9+i],":");
							fout<<vsGT[indDP]<<'\t';
							GTTmp.push_back(vsGT[indGT]);
							vector<string> vsPL=split(vsGT[indPL],",");
							for(int j=0;j<3;j++)
							{
								PLTmp.push_back(vsPL[j]);
								double lkTmp=atoi(vsPL[j].c_str());
								lkTmp=pow(10.0,-lkTmp/10.0);
								LKTmp(mapV2P[i],j)=lkTmp;
							}
						}
					}
				}
				for(size_t i=0;i<PLTmp.size();i++)
				{
					fout<<PLTmp[i]<<'\t';
				}

				fam.set_LK(LKTmp);
				//BN
				if(method==1)
				{
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						cout<<"This variant hasn't bee calculated: "<<endl;
						cout<<fline<<endl;
						//cout<<LKTmp<<endl;
						for(unsigned int i=0;i<fam.get_realNumInd();i++)
						{
							for(int j=0;j<3;j++)
							{
								fout<<"NA\t";
							}
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<GTTmp[i]<<'\t';
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<"NA\t";
						}
						fout<<endl;
						continue;
					}
				}

				//peeling
				if(method==2)
				{
					if(!fam.calPostProbBN(KnowSNP,chrType))
					{
						cout<<"This variant hasn't bee calculated: "<<endl;
						cout<<fline<<endl;
						//cout<<LKTmp<<endl;
						for(unsigned int i=0;i<fam.get_realNumInd();i++)
						{
							for(int j=0;j<3;j++)
							{
								fout<<"NA\t";
							}
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<GTTmp[i]<<'\t';
						}
						for(size_t i=0;i<GTTmp.size();i++)
						{
							fout<<"NA\t";
						}
						fout<<endl;
						continue;
					}
				}

				dMatrix<double> postProb=fam.get_postProb();
				for(int i=0;i<postProb.get_row();i++)
				{
					for(int j=0;j<postProb.get_column();j++)
					{
						fout<<postProb(i,j)<<'\t';
					}
				}

				for(size_t i=0;i<GTTmp.size();i++)
				{
					if(GTTmp[i]=="0/0")
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
					else if(GTTmp[i]=="0/1")
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else if(GTTmp[i]=="1/1")
					{
						fout<<vsLine[4][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else
					{
						fout<<"NA\t";
					}
				}

				vector<int> postRlt=fam.get_postRlt();
				for(size_t i=0;i<postRlt.size();i++)
				{
					if(postRlt[i]==0)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[3][0]<<'\t';
					}
					else if(postRlt[i]==1)
					{
						fout<<vsLine[3][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else if(postRlt[i]==2)
					{
						fout<<vsLine[4][0]<<"/"<<vsLine[4][0]<<'\t';
					}
					else
					{
						fout<<"NA\t";
					}
				}
				fout<<endl;
			}
		}
	}
	*/

	fin.close();
	fin.clear();
	fout.close();
	fout.clear();
	return true;
}

bool callGenoLK(const inputLK & lkInfo, family fam)
{
	string lkFileName=lkInfo.get_fileNameLK();
	string outFileName=lkInfo.get_outputName();
	int method=lkInfo.get_method();
	int lkType=lkInfo.get_lkType();
	vector<double> gProbN=fam.get_genoProbN();
	vector<double> gProbK=fam.get_genoProbK();
	vector<double> gProbXN=fam.get_genoProbXN();
	vector<double> gProbXK=fam.get_genoProbXK();

	ifstream fin(lkFileName.c_str());
	if(!fin.is_open())
	{
		cout<<"Cannot open "<<lkFileName<<endl;
		return false;
	}

	ofstream fout(outFileName.c_str());

	string title;

	getline(fin,title);

	fout<<"##FORMAT=<ID=GPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calculated by individual-base Method\">"<<endl;
	fout<<"##FORMAT=<ID=FPP,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled for posterior probability calculated by FamSeqPro\">"<<endl;
	fout<<"##FORMAT=<ID=FGT,Number=1,Type=String,Description=\"Genotype called by FamSeqPro\">"<<endl;

	fout<<"##FS mutation rate="<<fam.get_mRate()<<" "<<endl;
	fout<<"##FS genotype frequency in pupulation: "<<gProbN[0]<<":"<<gProbN[1]<<":"<<gProbN[2]<<endl;

	vector<string> vsTitle=split(title,"\t");
	vector<int> mapV2P(vsTitle.size(),-1), mapP2V(fam.get_numInd(),-1);
	vector<individual> mem=fam.get_member();

	for(size_t i=0;i<mapV2P.size();i++)
	{
		for(size_t j=0;j<mapP2V.size();j++)
		{
			if(vsTitle[i]==mem[j].get_sName())
			{
				mapV2P[i]=j;
				mapP2V[j]=i;
				break;
			}
		}
	}

	fam.set_mapP2V(mapP2V);
	fam.set_mapV2P(mapV2P);

	fout<<"#FORMAT\t";
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			fout<<vsTitle[i]<<'\t';
		}
	}
	fout<<endl;

	while(!fin.eof())
	{
		string fline;
		getline(fin,fline);
		if(fline.size()<2)
		{
			break;
		}

		vector<string> vsLine=split(fline,"\t");
		dMatrix<double> LKTmp(fam.get_numInd(),3,1);
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				vector<string> vsPL=split(vsLine[i],",");
				for(int j=0;j<3;j++)
				{
					double lkTmp=atof(vsPL[j].c_str());
					if(lkType==1)
					{
						LKTmp(mapV2P[i],j)=lkTmp;
					}
					else if(lkType==2)
					{
						lkTmp=pow(10.0,lkTmp);
						LKTmp(mapV2P[i],j)=lkTmp;
					}
					else if(lkType==3)
					{
						lkTmp=exp(lkTmp);
						LKTmp(mapV2P[i],j)=lkTmp;
					}
					else if(lkType==4)
					{
						lkTmp=pow(10.0,-lkTmp/10.0);
						LKTmp(mapV2P[i],j)=lkTmp;
					}
				}
			}
		}

		fam.set_LK(LKTmp);

		//different method
		//BN

		fout<<"LK:GPP:FPP:FGT\t";
		if(method==1)
		{
			if(!fam.calPostProbBN())
			{
				cout<<"Warning: this variant hasn't been calculated: "<<endl;
				cout<<fline<<endl;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[i]<<":NA:NA:NA\t";
					}
				}
				fout<<endl;
				continue;
			}
		}
		else if(method==2)
		{
			if(!fam.calPostProbPeeling())
			{
				cout<<"Warning: this variant hasn't been calculated: "<<endl;
				cout<<fline<<endl;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[i]<<":NA:NA:NA\t";
					}
				}
				fout<<endl;
				continue;
			}
		}
		else if(method==3)
		{
			if(!fam.calPostProbMCMC(lkInfo.get_numBurinIn(),lkInfo.get_numRep()))
			{
				cout<<"Warning: this variant hasn't been calculated: "<<endl;
				cout<<fline<<endl;
				for(size_t i=0;i<mapV2P.size();i++)
				{
					if(mapV2P[i]>=0)
					{
						fout<<vsLine[i]<<":NA:NA:NA\t";
					}
				}
				fout<<endl;
				continue;
			}
		}
		else
		{
		}

		dMatrix<double> postProb=fam.get_postProb();
		vector<int> postRlt=fam.get_postRlt();
		dMatrix<double> postProbSingle=fam.get_postProbSingle();

		int indTmp=0;
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				fout<<vsLine[i]<<":";
				if(-10*log10(postProbSingle(indTmp,0))==numeric_limits<double>::infinity())
				{
					fout<<99999<<',';
				}
				else
				{
					fout<<fabs(-10*log10(postProbSingle(indTmp,0)))<<',';
				}
				if(-10*log10(postProbSingle(indTmp,1))==numeric_limits<double>::infinity())
				{
					fout<<99999<<',';
				}
				else
				{
					fout<<fabs(-10*log10(postProbSingle(indTmp,1)))<<',';
				}
				if(-10*log10(postProbSingle(indTmp,2))==numeric_limits<double>::infinity())
				{
					fout<<99999<<':';
				}
				else
				{
					fout<<fabs(-10*log10(postProbSingle(indTmp,2)))<<':';
				}
				if(-10*log10(postProb(indTmp,0))==numeric_limits<double>::infinity())
				{
					fout<<99999<<',';
				}
				else
				{
					fout<<fabs(-10*log10(postProb(indTmp,0)))<<',';
				}
				if(-10*log10(postProb(indTmp,1))==numeric_limits<double>::infinity())
				{
					fout<<99999<<',';
				}
				else
				{
					fout<<fabs(-10*log10(postProb(indTmp,1)))<<',';
				}
				if(-10*log10(postProb(indTmp,2))==numeric_limits<double>::infinity())
				{
					fout<<99999<<':';
				}
				else
				{
					fout<<fabs(-10*log10(postProb(indTmp,2)))<<':';
				}
				if(postRlt[indTmp]==0)
				{
					fout<<"0/0\t";
				}
				else if(postRlt[indTmp]==1)
				{
					fout<<"0/1\t";
				}
				else
				{
					fout<<"1/1\t";
				}
				indTmp++;
			}
		}
		fout<<endl;
	}

	fin.close();
	fin.clear();
	fout.close();
	fout.clear();

	return true;
}

bool setFam(const inputVCF & vcfInfo, family &fam)
{
	vector<individual> mem;
	if(!readPed(vcfInfo.get_fileNamePed(),mem))
	{
		cout<<"Cannot read Ped file: "<<vcfInfo.get_fileNamePed()<<"."<<endl;
		return false;
	}

	fam=family(mem,vcfInfo.get_mRate());

	if(vcfInfo.get_genoProbN().size()!=0)
	{
		fam.set_genoProbN(vcfInfo.get_genoProbN());
	}

	if(vcfInfo.get_genoProbK().size()!=0)
	{
		fam.set_genoProbK(vcfInfo.get_genoProbK());
	}

	if(vcfInfo.get_genoProbXN().size()!=0)
	{
		fam.set_genoProbXN(vcfInfo.get_genoProbXN());
	}

	if(vcfInfo.get_genoProbXK().size()!=0)
	{
		fam.set_genoProbXK(vcfInfo.get_genoProbXK());
	}

	fam.set_lc(vcfInfo.get_LCriteria());

	if(!fam.init())
	{
		cout<<"Cannot initiate family. Please check ped file."<<endl;
		return false;
	}
	return true;
}

bool setFam(const inputLK & lkInfo, family & fam)
{
	vector<individual> mem;
	if(!readPed(lkInfo.get_fileNamePed(),mem))
	{
		cout<<"Cannot read Ped file: "<<lkInfo.get_fileNamePed()<<"."<<endl;
		return false;
	}

	fam=family(mem,lkInfo.get_mRate());

	if(lkInfo.get_genoProbN().size()!=0)
	{
		fam.set_genoProbN(lkInfo.get_genoProbN());
	}

	if(lkInfo.get_genoProbK().size()!=0)
	{
		fam.set_genoProbK(lkInfo.get_genoProbK());
	}

	if(lkInfo.get_genoProbXN().size()!=0)
	{
		fam.set_genoProbXN(lkInfo.get_genoProbXN());
	}

	if(lkInfo.get_genoProbXK().size()!=0)
	{
		fam.set_genoProbXK(lkInfo.get_genoProbXK());
	}

	fam.set_lc(lkInfo.get_LCriteria());

	if(!fam.init())
	{
		cout<<"Cannot initiate family. Please check ped file."<<endl;
		return false;
	}
	return true;
}
