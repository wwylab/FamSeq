/*
 * family.cu
 * GPU Version
 *
 *  Created on: Mar 7, 2012
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#include <cmath>
#include <algorithm>
#include <numeric>
#include <stdlib.h>
#include <time.h>

#include "family.h"

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cuda.h"

using namespace std;

#define NUM_THREAD 512
#define NUM_BLOCK 8

individual::individual(int id, int gender, int age, int state, int mId, int fId, string sName)
{
	m_id=id;
	m_gender=gender;
	m_age=age;
	m_state=state;
	m_mId=mId;
	m_fId=fId;
	m_sName=sName;
}

bool individual::set_id(int id)
{
	m_id=id;
	return true;
}

bool individual::set_age(int age)
{
	m_age=age;
	return true;
}

bool individual::set_gender(int gender)
{
	m_gender=gender;
	return true;
}

bool individual::set_sate(int state)
{
	m_state=state;
	return true;
}

bool individual::set_mId(int mId)
{
	m_mId=mId;
	return true;
}

bool individual::set_fId(int fId)
{
	m_fId=fId;
	return true;
}

bool individual::set_sName(string sName)
{
	m_sName=sName;
	return true;
}

family::family()
{
}

family::family(const vector<individual> & mem, double mutationRate)
{
	member=mem;
	numInd=member.size();

	//fulfill the whole family
	//realNumInd=numInd;

	//fulfillFamily();

	mRate=mutationRate;

	//set Pr(genotype)
	genoProbN.clear();
	/*
	genoProbN.push_back(0.999*0.999);
	genoProbN.push_back(2*0.999*0.001);
	genoProbN.push_back(0.000001);
	*/
	genoProbN.push_back(0.9985);
	genoProbN.push_back(0.001);
	genoProbN.push_back(0.0005);

	genoProbK.clear();
	genoProbK.push_back(0.45);
	genoProbK.push_back(0.1);
	genoProbK.push_back(0.45);

	genoProbXN.clear();
	genoProbXN.push_back(0.999);
	genoProbXN.push_back(0);
	genoProbXN.push_back(0.001);

	genoProbXK.clear();
	genoProbXK.push_back(0.5);
	genoProbXK.push_back(0);
	genoProbXK.push_back(0.5);

	postProb=dMatrix<double> (numInd,3);
	postProbSingle=dMatrix<double> (numInd,3);

	likelihood=dMatrix<double> (numInd,3);

	//set indicators
	flagInit=false;
	flagLK=false;
	flagPB=false;
	flagPBS=false;
	flagGender=false;
}

family & family::operator =(const family & fm)
{
	if(this==&fm)
	{
		return *this;
	}

	member=fm.get_member();
	numInd=fm.get_numInd();
	realNumInd=fm.get_realNumInd();
	child=fm.get_child();
	parent=fm.get_parent();
	spouse=fm.get_spouse();
	mRate=fm.get_mRate();
	flagInit=fm.get_flagInit();
	flagLK=fm.get_flagLK();
	flagPB=fm.get_flagPB();
	flagPBS=fm.get_flagPBS();
	flagGender=fm.get_flagGender();
	genoProbN=fm.get_genoProbN();
	genoProbK=fm.get_genoProbK();
	genoProbXN=fm.get_genoProbXN();
	genoProbXK=fm.get_genoProbXK();
	pcp2=fm.get_pcp2();
	pcp2Xf=fm.get_pcp2Xf();
	pcp2Xm=fm.get_pcp2Xm();
	postProb=fm.get_postProb(false);
	postProbSingle=fm.get_postProbSingle(false);
	likelihood=fm.get_LK();
	mapP2V=fm.get_mapP2V();
	mapV2P=fm.get_mapV2P();

	return *this;
}

bool family::fulfillFamily()
{
	for(unsigned int i=0;i<realNumInd;i++)
	{
		int fid=member[i].get_fId();
		int mid=member[i].get_mId();

		int indM=-1;
		int indF=-1;
		for(unsigned int j=0;j<realNumInd;j++)
		{
			if(member[j].get_id()==mid)
			{
				indM=j;
			}
			if(member[j].get_id()==fid)
			{
				indF=j;
			}
		}

		if(indM>=0 && indF<0)
		{
			numInd=numInd+1;
			member[i].set_fId(-(int)numInd);
			individual indTmp(-(int)numInd,1,0,0,0,0,"NA");
			member.push_back(indTmp);
		}

		if(indF>=0 && indM<0)
		{
			numInd=numInd+1;
			member[i].set_mId(-((int)numInd));
			individual indTmp(-((int)numInd),2,0,0,0,0,"NA");
			member.push_back(indTmp);
		}
	}
	return true;
}

bool family::checkPed(){
	for(size_t i=0;i<parent.size();i++){
		if(parent[i].size()==2){
			if(member[parent[i][0]].get_gender() != 2){
				cerr<<"Sample "<<member[parent[i][0]].get_id()<<"'s a mother while she is not a female."<<endl;
				return false;
			}

			if(member[parent[i][1]].get_gender() != 1){
				cerr<<"Sample "<<member[parent[i][0]].get_id()<<"'s a father while she is not a male."<<endl;
				return false;
			}
		}
	}
	return true;
}

bool family::init()
{
	if(!setPCP())
	{
		return false;
	}

	if(!setRelation())
	{
		return false;
	}

	if(!checkPed()){
		return false;
	}

	return true;
}

bool family::set_mRate(double rate)
{
	if(mRate!=rate)
	{
		mRate=rate;
		if(!setPCP())
		{
			return false;
		}
	}
	return true;
}

bool family::set_lc(double lc)
{
	m_lc = lc;
	return true;
}

bool family::setPCP()
{
	//pcp2
	if(!calPCP2())
	{
		return false;
	}

	if(!calPCP2Xf())
	{
		return false;
	}

	if(!calPCP2Xm())
	{
		return false;
	}

	/*
	for(size_t i=0;i<pcp2Xf.size();i++)
	{
		cout<<pcp2Xf[i]<<endl;
	}

	for(size_t i=0;i<pcp2Xm.size();i++)
	{
		cout<<pcp2Xm[i]<<endl;
	}
	*/
	return true;
}

bool family::setRelation()
{
	child.clear();
	parent.clear();
	spouse.clear();

	child.resize(numInd);
	parent.resize(numInd);
	spouse.resize(numInd);

	for(unsigned int i=0;i<numInd;i++)
	{
		int mId=member[i].get_mId();
		int fId=member[i].get_fId();
		int indM=-1;
		int indF=-1;
		for(unsigned int j=0;j<numInd;j++)
		{
			if(mId==member[j].get_id())
			{
				indM=j;
			}
			if(fId==member[j].get_id())
			{
				indF=j;
			}
		}

		if((indM<0 && indF>=0) || (indM>=0 && indF<0))
		{
			cout<<"This is not a fulfill family. Please check the ped file."<<endl;
			return false;
		}

		if(indM>=0 && indF>=0)
		{
			//first one is mother
			//second one is father
			parent[i].push_back(indM);
			parent[i].push_back(indF);
			child[indM].push_back(i);
			child[indF].push_back(i);
			bool notFind=true;
			for(vector<int>::size_type j=0;j<spouse[indM].size();j++)
			{
				if(indF==spouse[indM][j])
				{
					notFind=false;
					break;
				}
			}
			if(notFind)
			{
				spouse[indM].push_back(indF);
				spouse[indF].push_back(indM);
			}
		}
	}
	return true;
}

bool family::set_fg(bool fg)
{
	flagGender=fg;
	return true;
}

bool family::set_mapP2V(vector<int> mP2V)
{
	mapP2V=mP2V;
	return true;
}

bool family::set_mapV2P(vector<int> mV2P)
{
	mapV2P=mV2P;
	realNumInd=0;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			realNumInd++;
		}
	}
	return true;
}

bool family::calPCP2()
{
	return calPCP2S(2,pcp2);
}

bool family::calPCP2Xf()
{
	pcp2Xf.clear();
	for(int i=0;i<3;i++)
	{
		dMatrix<double> dMTmp(3,3,0);
		pcp2Xf.push_back(dMTmp);
	}

	pcp2Xf[0](0,0)=(1.0-mRate)*(1.0-mRate);
	pcp2Xf[1](0,0)=2*mRate*(1.0-mRate);
	pcp2Xf[2](0,0)=mRate*mRate;

	pcp2Xf[0](0,2)=(1.0-mRate)*mRate;
	pcp2Xf[1](0,2)=(1.0-mRate)*(1.0-mRate)+mRate*mRate;
	pcp2Xf[2](0,2)=(1.0-mRate)*mRate;

	pcp2Xf[0](1,0)=(1.0-mRate)*(1.0-mRate)/2+mRate*(1.0-mRate)/2;
	pcp2Xf[1](1,0)=mRate*(1.0-mRate)+(1.0-mRate)*(1.0-mRate)/2+mRate*mRate/2;
	pcp2Xf[2](1,0)=mRate*mRate/2+mRate*(1.0-mRate)/2;

	pcp2Xf[0](1,2)=mRate*mRate/2+mRate*(1.0-mRate)/2;
	pcp2Xf[1](1,2)=mRate*(1.0-mRate)+(1.0-mRate)*(1.0-mRate)/2+mRate*mRate/2;
	pcp2Xf[2](1,2)=(1.0-mRate)*(1.0-mRate)/2+mRate*(1.0-mRate)/2;

	pcp2Xf[0](2,0)=(1.0-mRate)*mRate;
	pcp2Xf[1](2,0)=(1.0-mRate)*(1.0-mRate)+mRate*mRate;
	pcp2Xf[2](2,0)=(1.0-mRate)*mRate;

	pcp2Xf[0](2,2)=mRate*mRate;
	pcp2Xf[1](2,2)=2*mRate*(1.0-mRate);
	pcp2Xf[2](2,2)=(1.0-mRate)*(1.0-mRate);
	return true;
}

bool family::calPCP2Xm()
{
	pcp2Xm.clear();
	for(int i=0;i<3;i++)
	{
		dMatrix<double> dMTmp(3,3,0);
		pcp2Xm.push_back(dMTmp);
	}

	pcp2Xm[0](0,0)=1-mRate;
	pcp2Xm[2](0,0)=mRate;

	pcp2Xm[0](0,2)=1-mRate;
	pcp2Xm[2](0,2)=mRate;

	pcp2Xm[0](1,0)=0.5;
	pcp2Xm[2](1,0)=0.5;

	pcp2Xm[0](1,2)=0.5;
	pcp2Xm[2](1,2)=0.5;

	pcp2Xm[0](2,0)=mRate;
	pcp2Xm[2](2,0)=1-mRate;

	pcp2Xm[0](2,2)=mRate;
	pcp2Xm[2](2,2)=1-mRate;
	return true;
}

bool family::calPCP2S(int nAllele, std::vector<dMatrix<double> > & PCP2S)
{
	int nGeno=(nAllele+1)*nAllele/2;
	PCP2S.clear();
	for(int i=0;i<nGeno;i++)
	{
		dMatrix<double> dMTmp(nGeno,nGeno,0);
		PCP2S.push_back(dMTmp);
	}

	dMatrix<int> codeTable(nAllele,nAllele,0);
	dMatrix<int> decodeTable(nGeno,2,0);
	int count=0;
	for(int i=0;i<nAllele;i++)
	{
		for(int j=i;j<nAllele;j++)
		{
			codeTable(i,j)=count;
			codeTable(j,i)=count;
			decodeTable(count,0)=i;
			decodeTable(count,1)=j;
			count++;
		}
	}

	if(mRate==0)
	{
		for(int i=0;i<nGeno;i++)
		{
			for(int j=0;j<nGeno;j++)
			{
				//haplotype
				int i0,i1,j0,j1;
				i0=decodeTable(i,0);
				i1=decodeTable(i,1);
				j0=decodeTable(j,0);
				j1=decodeTable(j,1);
				PCP2S[codeTable(i0,j0)](i,j)=PCP2S[codeTable(i0,j0)](i,j)+0.25;
				PCP2S[codeTable(i0,j1)](i,j)=PCP2S[codeTable(i0,j1)](i,j)+0.25;
				PCP2S[codeTable(i1,j0)](i,j)=PCP2S[codeTable(i1,j0)](i,j)+0.25;
				PCP2S[codeTable(i1,j1)](i,j)=PCP2S[codeTable(i1,j1)](i,j)+0.25;
			}
		}
	}
	else
	{
		for(int i=0;i<nGeno;i++)
		{
			for(int j=0;j<nGeno;j++)
			{
				vector<double> p1(nAllele,mRate/(2*(nAllele-1)));
				vector<double> p2(nAllele,mRate/(2*(nAllele-1)));

				//haplotype 0 and 0
				p1[decodeTable(i,0)]=(1-mRate)/2;
				p2[decodeTable(j,0)]=(1-mRate)/2;
				for(int k=0;k<nAllele;k++)
				{
					for(int l=0;l<nAllele;l++)
					{
						PCP2S[codeTable(k,l)](i,j)=PCP2S[codeTable(k,l)](i,j)+p1[k]*p2[l];
					}
				}

				//haplotype 0 and 1
				p2[decodeTable(j,0)]=mRate/(2*(nAllele-1));
				p2[decodeTable(j,1)]=(1-mRate)/2;
				for(int k=0;k<nAllele;k++)
				{
					for(int l=0;l<nAllele;l++)
					{
						PCP2S[codeTable(k,l)](i,j)=PCP2S[codeTable(k,l)](i,j)+p1[k]*p2[l];
					}
				}

				//haplotype 1 and 0
				p1[decodeTable(i,0)]=mRate/(2*(nAllele-1));
				p2[decodeTable(j,1)]=mRate/(2*(nAllele-1));
				p1[decodeTable(i,1)]=(1-mRate)/2;
				p2[decodeTable(j,0)]=(1-mRate)/2;
				for(int k=0;k<nAllele;k++)
				{
					for(int l=0;l<nAllele;l++)
					{
						PCP2S[codeTable(k,l)](i,j)=PCP2S[codeTable(k,l)](i,j)+p1[k]*p2[l];
					}
				}

				//haplotype 1 and 1
				p2[decodeTable(j,0)]=mRate/(2*(nAllele-1));
				p2[decodeTable(j,1)]=(1-mRate)/2;
				for(int k=0;k<nAllele;k++)
				{
					for(int l=0;l<nAllele;l++)
					{
						PCP2S[codeTable(k,l)](i,j)=PCP2S[codeTable(k,l)](i,j)+p1[k]*p2[l];
					}
				}
			}
		}
	}

	return true;
}

bool family::set_genoProbK(const vector<double> & prob)
{
	genoProbK=prob;
	return true;
}

bool family::set_genoProbN(const vector<double> & prob)
{
	genoProbN=prob;
	return true;
}

bool family::set_genoProbXK(const vector<double> & prob)
{
	genoProbXK=prob;
	return true;
}

bool family::set_genoProbXN(const vector<double> & prob)
{
	genoProbXN=prob;
	return true;
}

dMatrix<double> family::get_postProb(bool flag) const
{
	if(flag)
	{
		if(!flagPB)
		{
			cout<<"The posterior probability (FamSeqPro) hasn't been calculated yet!"<<endl;
			return dMatrix<double> (1,1,-1);
		}
		dMatrix<double> rlt(realNumInd,3);
		int ind=0;
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				for(int j=0;j<3;j++)
				{
					rlt(ind,j)=postProb(mapV2P[i],j);
				}
				ind++;
			}
		}
		return rlt;
	}
	else
	{
		return postProb;
	}
}

dMatrix<double> family::get_postProbSingle(bool flag) const
{
	if(flag)
	{
		if(!flagPBS)
		{
			cout<<"The posterior probability (Single) hasn't been calculated yet!"<<endl;
			return dMatrix<double> (1,1,-1);
		}
		dMatrix<double> rlt(realNumInd,3);
		int ind=0;
		for(size_t i=0;i<mapV2P.size();i++)
		{
			if(mapV2P[i]>=0)
			{
				for(int j=0;j<3;j++)
				{
					rlt(ind,j)=postProbSingle(mapV2P[i],j);
				}
				ind++;
			}
		}
		return rlt;
	}
	else
	{
		return postProbSingle;
	}
}

vector<int> family::get_postRlt() const
{
	if(!flagPB)
	{
		cout<<"The posterior probability (FamSeqPro) hasn't been calculated yet!"<<endl;
		vector<int> tmp;
		return tmp;
	}

	vector<int> rlt;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double bigTmp=-1;
			int ind=-1;
			for(int j=0;j<3;j++)
			{
				if(bigTmp<postProb(mapV2P[i],j))
				{
					bigTmp=postProb(mapV2P[i],j);
					ind=j;
				}
			}
			rlt.push_back(ind);
		}
	}

	return rlt;
}

vector<int> family::get_postRltSingle() const
{
	if(!flagPBS)
	{
		cout<<"The posterior probability (Single) hasn't been calculated yet!"<<endl;
		vector<int> tmp;
		return tmp;
	}

	vector<int> rlt;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double bigTmp=-1;
			int ind=-1;
			for(int j=0;j<3;j++)
			{
				if(bigTmp<postProbSingle(mapV2P[i],j))
				{
					bigTmp=postProbSingle(mapV2P[i],j);
					ind=j;
				}
			}
			rlt.push_back(ind);
		}
	}

	return rlt;
}

bool family::set_LK(const dMatrix<double> & lk)
{
	/*
	if(lk.get_row()!=int(numInd))
	{
		for(int i=0;i<lk.get_row();i++)
		{
			for(int j=0;j<3;j++)
			{
				likelihood(i,j)=lk(i,j);
			}
		}
		for(unsigned int i=lk.get_row();i<numInd;i++)
		{
			for(int j=0;j<3;j++)
			{
				likelihood(i,j)=1;
			}
		}
	}
	else
	{
		for(int i=0;i<lk.get_row();i++)
		{
			for(int j=0;j<3;j++)
			{
				likelihood(i,j)=lk(i,j);
			}
		}
	}
	*/

	if(int(numInd)!=lk.get_row() || lk.get_column()!=3)
	{
		cout<<"The dimention of likelihood matrix is wrong. Cannot set likelihood."<<endl;
		return false;
	}

	for(int i=0;i<lk.get_row();i++)
	{
		for(int j=0;j<3;j++)
		{
			likelihood(i,j)=lk(i,j);
		}
	}

	flagLK=true;
	flagPB=false;
	flagPBS=false;
	return true;
}

bool decodeGeno(int genoCode, int numSample, vector<int> & geno){
	geno.resize(numSample);
	for(int i=0;i<numSample;i++){
		int rm = genoCode%3;
		geno[i] = rm;
		genoCode = (genoCode-rm)/3;
	}
	return true;
}


extern "C" __global__ void calPostProb(int* father, int* mother, double* lk, double* pcp2, double* postProb, double * genoFry, int numSample,int numGeno){
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
	double pTmp[60]={0};
	for(int i=bid*NUM_THREAD+tid;i<numGeno;i+=NUM_BLOCK*NUM_THREAD){
		int geno[20];
		int genoCode=i;
		for(int j=0;j<numSample;j++){
			int rm=genoCode%3;
			geno[j]=rm;
			genoCode=(genoCode-rm)/3;
		}

		double postTmp = 10000000;
		for(int j=0;j<numSample;j++){
			if(father[j]<0){
				postTmp = postTmp*genoFry[geno[j]]*lk[j*3+geno[j]];
			}
			else{
				int indexpcp = geno[j]*9+geno[mother[j]]*3+geno[father[j]];
				postTmp = postTmp*pcp2[indexpcp]*lk[j*3+geno[j]];
			}
		}
		for(int j=0;j<numSample;j++){
			pTmp[j*3+geno[j]] += postTmp;
		}
	}

	for(int i=0;i<numSample*3;i++){
		postProb[(bid*NUM_THREAD+tid)*numSample*3+i]=pTmp[i];
	}
}



/*
extern "C" __global__ void calPostProb(int* father, int* mother, double* lk, double* pcp2, double* postProb, double * genoFry, int numSample,int numGeno){
	for(int i=blockIdx.x;i<numGeno;i+=NUM_BLOCK){
		int geno[20];
		int genoCode=i;
		for(int j=0;j<numSample;j++){
			int rm=genoCode%3;
			geno[j]=rm;
			genoCode=(genoCode-rm)/3;
		}

		double postTmp = 10000000;
		for(int j=0;j<numSample;j++){
			if(father[j]<0){
				postTmp = postTmp*genoFry[geno[j]]*lk[j*3+geno[j]];
			}
			else{
				int indexpcp = geno[j]*9+geno[mother[j]]*3+geno[father[j]];
				postTmp = postTmp*pcp2[indexpcp]*lk[j*3+geno[j]];
			}
		}

		for(int j=0;j<numSample;j++){
			postProb[blockIdx.x*numSample*3+j*3+geno[j]]=postTmp+postProb[blockIdx.x*numSample*3+j*3+geno[j]];
		}
	}
}
*/


/*
extern "C" __global__ void calPostProb(int* father, int* mother, double * lk, double * pcp2, double * postProb, double * genoFry, int numSample, int numGeno){
	
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	if(index < numGeno){
		int geno[20];
		int genoCode = index;
		for(int i=0;i<numSample;i++){
			int rm = genoCode%3;
			geno[i] = rm;
			genoCode = (genoCode-rm)/3;
		}
		double postTmp = 10000000;
		for(int i=0;i<numSample;i++){
			if(father[i]<0){
				postTmp = postTmp*genoFry[geno[i]]*lk[i*3+geno[i]];
			}
			else{
				int indexpcp = geno[i]*9+geno[mother[i]]*3+geno[father[i]];
				postTmp = postTmp*pcp2[indexpcp]*lk[i*3+geno[i]];
			}
		}

		//for(int i=0;i<numSample;i++){
			//postProb[i*3+geno[i]]=postProb[i*3+geno[i]]+1;
			//atomicAdd((postProb+(i*3+geno[i])),float(postTmp));
			//atomicAdd((int*) postProb,int(postTmp));

		//}

		postProb[index] = postTmp;
		//postProb[index] = postTmp;
	}
}
*/

/*
extern "C" __global__ void sumPostProb(int numGeno, int numSample, double * postProb, double * postProbS){
	for(int i=0;i<numGeno;i++){
		int geno[100];
		int genoCode = i;
		for(int j=0;j<numSample;j++){
			int rm = genoCode%3;
			geno[j] = rm;
			genoCode = (genoCode-rm)/3;
		}

		for(int j=0;j<numSample;j++){
			postProbS[j*3+geno[j]]=postProbS[j*3+geno[j]]+postProb[i];
		}
	}
}
*/


extern "C" __global__ void sumPostProb(int numGeno, int numSample, double * postProb, double * postProbS){
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;

	double pTmp[60];
	for(int i=bid*NUM_THREAD+tid;i<numGeno;i+=NUM_BLOCK*NUM_THREAD){
		int geno[20];
		int genoCode=i;
		for(int j=0;j<numSample;j++){
			int rm=genoCode%3;
			geno[j]=rm;
			genoCode=(genoCode-rm)/3;
		}
		for(int j=0;j<numSample;j++){
			pTmp[j*3+geno[j]]+=postProb[i];
		}
	}
	for(int i=0;i<numSample;i++){
		for(int j=0;j<3;j++){
			postProbS[(bid*NUM_THREAD+tid)*numSample*3+i*3+j] = pTmp[i*3+j];
		}
	}

	/*
	for(int i=blockIdx.x;i<numGeno;i+=NUM_BLOCK){
		int geno[20];
		int genoCode=i;
		for(int j=0;j<numSample;j++){
			int rm=genoCode%3;
			geno[j]=rm;
			genoCode=(genoCode-rm)/3;
		}
		for(int j=0;j<numSample;j++){
			postProbS[blockIdx.x*numSample*3+j*3+geno[j]] = postProbS[blockIdx.x*numSample*3+j*3+geno[j]]+postProb[i];
		}
	}
	*/
}


extern "C" __global__ void calPostProbX(int * gender, int* father, int* mother, double * lk, double * pcp2M, double* pcp2F, double * postProb, double * genoFryM, double* genoFryF, int numSample, int numGeno){
	const int tid = threadIdx.x;
	const int bid = blockIdx.x;
	double pTmp[60]= {0};
	for(int i=bid*NUM_THREAD+tid;i<numGeno;i+=NUM_BLOCK*NUM_THREAD){
		int geno[20];
		int genoCode=i;
		for(int j=0;j<numSample;j++){
			int rm=genoCode%3;
			geno[j]=rm;
			genoCode=(genoCode-rm)/3;
		}

		double postTmp = 10000000;
		for(int j=0;j<numSample;j++){
			if(father[j]<0){
				if(gender[j]==1){
					postTmp = postTmp*genoFryM[geno[j]]*lk[j*3+geno[j]];
				}
				else{
					postTmp = postTmp*genoFryF[geno[j]]*lk[j*3+geno[j]];
				}
			}
			else{
				int indexpcp = geno[j]*9+geno[mother[j]]*3+geno[father[j]];
				if(gender[j]==1){
					postTmp = postTmp*pcp2M[indexpcp]*lk[j*3+geno[j]];
				}
				else{
					postTmp = postTmp*pcp2F[indexpcp]*lk[j*3+geno[j]];
				}	
			}
		}
		
		for(int j=0;j<numSample;j++){
			pTmp[j*3+geno[j]] += postTmp;
		}
	}
		
	for(int i=0;i<numSample*3;i++){
		postProb[(bid*NUM_THREAD+tid)*numSample*3+i]=pTmp[i];
	}
}

#if 1
bool family::calPostProbBN(bool Known, int chrType)
{
	if(!flagLK)
	{
		cout<<"Likelihood has not been set. Please set likelihood first."<<endl;
		return false;
	}

	if(!calPostProbSingle(Known,chrType))
	{
		//cout<<Known<<'\t'<<chrType<<endl;
		//cout<<likelihood<<endl;
		return false;
	}

	//cout<<likelihood<<endl;

	bool flagSingle=true;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double big=0;
			double sum=0;
			for(int j=0;j<3;j++)
			{
				if(big<likelihood(mapV2P[i],j))
				{
					big=likelihood(mapV2P[i],j);
				}
				sum=sum+likelihood(mapV2P[i],j);
			}
			big=big/sum;
			if(big<=m_lc)
			{
				flagSingle=false;
				break;
			}
		}
	}

	postProb.setZero();

	if(flagSingle)
	{
		if(chrType==0)
		{
			vector<double> genoProb;
			if(Known)
			{
				genoProb=genoProbK;
			}
			else
			{
				genoProb=genoProbN;
			}
			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					postProb(i,j)=likelihood(i,j)*genoProb[j];
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		else if(chrType==1)
		{
			vector<double> genoProbM;
			vector<double> genoProbF;
			if(Known)
			{
				genoProbM=genoProbXK;
			}
			else
			{
				genoProbM=genoProbXN;
			}
			if(Known)
			{
				genoProbF=genoProbK;
			}
			else
			{
				genoProbF=genoProbN;
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					if(member[i].get_gender()==1)
					{
						postProb(i,j)=likelihood(i,j)*genoProbM[j];
					}
					else
					{
						postProb(i,j)=likelihood(i,j)*genoProbF[j];
					}
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		flagPB=true;
		return true;
	}



	if(chrType==0)
	{
		/*
		set up data that will be used in BN
		*/
		vector<int> host_father(numInd,-1);
		vector<int> host_mother(numInd,-1);
		for(size_t i=0;i<numInd;i++){
			if(parent[i].size()!=0){
				host_mother[i]=parent[i][0];
				host_father[i]=parent[i][1];
			}
		}

		vector<double> host_lk(numInd*3,0);
		for(size_t i=0;i<numInd;i++){
			host_lk[i*3] = likelihood(i,0);
			host_lk[i*3+1] = likelihood(i,1);
			host_lk[i*3+2] = likelihood(i,2);
		}

		vector<double> host_pcp2(27,0);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					host_pcp2[i*9+j*3+k] = pcp2[i](j,k);
				}
			}
		}

		int numGeno = (int)pow(3.0,(double)numInd);

		cudaError_t cudaStatus;

		int * dev_father = 0;
		int * dev_mother = 0;
		double * dev_lk = 0;
		double * dev_pcp2 =0;
		//double * dev_postProb = 0;
		double * dev_genoFry = 0;
		double * dev_postProb = 0;
		//double * dev_postProbS = 0;


		
		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		// Allocate GPU buffers
		cudaStatus = cudaMalloc((void**)&dev_father, host_father.size() * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_mother, host_mother.size() * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_lk, host_lk.size() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);
			
			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_pcp2, host_pcp2.size() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		/*
		cudaStatus = cudaMalloc((void**)&dev_postProb, (int)pow(3.0,(double)numInd) * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			delete host_postProb;
			
			flagPB=false;
			return false;
		}
		*/

		cudaStatus = cudaMalloc((void**)&dev_genoFry, 3*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		
		cudaStatus = cudaMalloc((void**)&dev_postProb, (NUM_THREAD*NUM_BLOCK)*3*numInd*sizeof(double));
		//cudaStatus = cudaMalloc((void**)&dev_postProb, numGeno*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);
			
			flagPB=false;
			return false;
		}
		

		/*
		cudaStatus = cudaMalloc((void**)&dev_postProbS, (NUM_THREAD*NUM_BLOCK)*3*numInd*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			delete host_postProb;
			
			flagPB=false;
			return false;
		}
		*/
		
		

		/*
		copy data
		*/
		cudaStatus = cudaMemcpy(dev_father, &host_father[0], host_father.size() * sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed father!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);
			
			flagPB=false;
			return false;
		}
		
		/*
		int *test_father = new int[host_father.size()];
		cudaMemcpy(test_father, dev_father, host_father.size()*sizeof(int), cudaMemcpyDeviceToHost);
		for(size_t i =0; i<host_father.size();i++){
			cout<<host_father[i]<<'\t'<<test_father[i]<<endl;
		}
		delete test_father;
		*/
		

		cudaStatus = cudaMemcpy(dev_mother, &host_mother[0], host_mother.size() * sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed mother!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		/*
		int *test_mother = new int[host_mother.size()];
		cudaMemcpy(test_mother, dev_mother, host_mother.size()*sizeof(int), cudaMemcpyDeviceToHost);
		for(size_t i =0; i<host_mother.size();i++){
			cout<<host_mother[i]<<'\t'<<test_mother[i]<<endl;
		}
		delete test_mother;
		*/

		cudaStatus = cudaMemcpy(dev_lk, &host_lk[0], host_lk.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed lk!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		/*
		double* test_lk = new double[host_lk.size()];
		cudaMemcpy(test_lk, dev_lk, host_lk.size()*sizeof(double), cudaMemcpyDeviceToHost);
		for(size_t i=0;i<host_lk.size();i++){
			cout<<test_lk[i]<<'\t'<<host_lk[i]<<endl;
		}
		delete test_lk;
		*/
		
		
		
		cudaStatus = cudaMemcpy(dev_pcp2, &host_pcp2[0], host_pcp2.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed pcp2!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}
		

		/*
		double* test_pcp2 = new double[host_pcp2.size()];
		cudaMemcpy(test_pcp2, dev_pcp2, host_pcp2.size()*sizeof(double),cudaMemcpyDeviceToHost);
		for(size_t i=0;i<host_pcp2.size();i++)
		{
			cout<<test_pcp2[i]<<'\t'<<host_pcp2[i]<<endl;
		}
		delete test_pcp2;
		*/

		if(Known){
			cudaStatus = cudaMemcpy(dev_genoFry, &genoProbK[0], 3*sizeof(double), cudaMemcpyHostToDevice);

			/*
			double* test_genoFre = new double[3];
			cudaMemcpy(test_genoFre, dev_genoFry, 3*sizeof(double), cudaMemcpyDeviceToHost);
			for(size_t i=0;i<3;i++){
				cout<<test_genoFre[i]<<'\t'<<genoProbK[i]<<endl;
			}
			delete test_genoFre;
			*/
		}
		else{
			cudaStatus = cudaMemcpy(dev_genoFry, &genoProbN[0], 3*sizeof(double), cudaMemcpyHostToDevice);

			/*
			double* test_genoFre = new double[3];
			cudaMemcpy(test_genoFre, dev_genoFry, 3*sizeof(double), cudaMemcpyDeviceToHost);
			for(size_t i=0;i<3;i++){
				cout<<test_genoFre[i]<<'\t'<<genoProbN[i]<<endl;
			}
			delete test_genoFre;
			*/
		}
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed PenoProb!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			
			flagPB=false;
			return false;
		}

		
		/*
		for(int i=0;i<3*numInd*NUM_BLOCK;i++){
			host_postProb[i]=0;
		}

		//for(int i=0;i<3*numInd;i++)
		//{
		//	cout<<host_postProb[i]<<endl;
		//}

		cudaStatus = cudaMemcpy(dev_postProbS, host_postProb, 3*numInd*NUM_BLOCK*sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed host_postProb!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			delete host_postProb;
			
			flagPB=false;
			return false;
		}
		*/

		//cudaStatus = cudaMemset(dev_postProb, 0, (NUM_THREAD*NUM_BLOCK)*3*numInd*sizeof(double));
		

		//for(int i=0;i<3*numInd*NUM_THREAD;i++)
		//{
		//	cout<<host_postProb[i]<<endl;
		//}
		//for(int i=0;i<3*numInd;i++)
		//{
		//	host_posteriorProb[i]=0;
		//}


		//int numThread = 512;
		//calPostProb<<<(numGeno + numThread-1)/numThread, numThread>>>(dev_father, dev_mother, dev_lk, dev_pcp2, dev_postProb, dev_genoFry, numInd, numGeno);

		//calPostProb<<<(numGeno + NUM_THREAD-1)/NUM_THREAD, NUM_THREAD>>>(dev_father, dev_mother, dev_lk, dev_pcp2, dev_postProb, dev_genoFry, numInd, numGeno);

/*
		// Allocate CUDA events that we'll use for timing
		cudaError_t error;
		cudaEvent_t start;
		error = cudaEventCreate(&start);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		cudaEvent_t stop;
		error = cudaEventCreate(&stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(start, NULL);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		calPostProb<<<NUM_BLOCK, NUM_THREAD>>>(dev_father, dev_mother, dev_lk, dev_pcp2, dev_postProb, dev_genoFry, numInd, numGeno);
		
		// Record the stop event
		error = cudaEventRecord(stop, NULL);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Wait for the stop event to complete
		error = cudaEventSynchronize(stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		float msecTotal = 0.0f;
		error = cudaEventElapsedTime(&msecTotal, start, stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Compute and print the performance
		printf("Time= %.3f msec\n",msecTotal);

*/
		//sumPostProb<<<1,1>>>(dev_postProb, dev_posteriorProb, numGeno, numInd);

		//sumPostProb<<<NUM_BLOCK,NUM_THREAD>>>(numGeno, numInd, dev_postProb, dev_postProbS);

		calPostProb<<<NUM_BLOCK, NUM_THREAD>>>(dev_father, dev_mother, dev_lk, dev_pcp2, dev_postProb, dev_genoFry, numInd, numGeno);
		
		double * host_postProb = new double[NUM_THREAD*NUM_BLOCK*3*numInd];
		cudaStatus = cudaMemcpy(host_postProb, dev_postProb, NUM_THREAD*NUM_BLOCK*3*numInd*sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed p to p!");
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2);
			//cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFry);
			cudaFree((void*)dev_postProb);
			//cudaFree((void*)dev_postProbS);

			delete host_postProb;
			
			flagPB=false;
			return false;
		}
		

		postProb.setZero();

		for(int i=0;i<NUM_BLOCK*NUM_THREAD;i++){
			for(int j=0;j<numInd;j++){
				for(int k=0;k<3;k++){
					postProb(j,k)=postProb(j,k)+host_postProb[i*numInd*3+j*3+k];
				}
			}
		}
		
		/*
		for(int i=0;i<NUM_BLOCK;i++){
			for(int j=0;j<numInd;j++){
				for(int k=0;k<3;k++){
					//postProb(j,k)=postProb(j,k)+host_postProb[i*3*numInd+j*3+k];
					postProb(j,k)=1;
				}
			}
		}
		*/
		
		

		/*
		int geno[20]={0};
		double* pTmp = new double[numInd*3];
		for(int i=0;i<numInd*3;i++){
			pTmp[i]=0;
		}
		for(int i=0;i<numGeno;i++){
			
			for(int j=0;j<numInd;j++){
				geno[j]++;
				if(geno[j]==3){
					geno[j]=0;
				}
				else{
					break;
				}
			}
			
			for(int j=0;j<numInd;j++){
				//postProb(j,geno[j])=postProb(j,geno[j])+host_postProb[i];
				pTmp[j*3+geno[j]] = pTmp[j*3+geno[j]]+host_postProb[i];
			}
		}
		
		for(int i=0;i<numInd;i++){
			for(int j=0;j<3;j++){
				postProb(i,j)=pTmp[i*3+j];
			}
		}
		delete[] pTmp;
		*/
		
		

		//for(int i=0;i<3*numInd*NUM_THREAD;i++)
		//{
		//	cout<<host_postProb[i]<<endl;
		//}
		/*
		for(unsigned int i=0;i<numInd;i++){
			for(int j=0;j<3;j++){
				//postProb(i,j) = host_postProb[i*3+j];
				postProb(i,j) = 1;
			}
		}
		*/

		

		for(unsigned int i=0;i<numInd;i++){
			double sTmp = 0;
			for(int j=0;j<3;j++){
				sTmp = sTmp + postProb(i,j);
				//cout<<postProb(i,j)<<'\t';
			}
			//cout<<endl;
			if(sTmp<=0)
			{
				return false;
			}
			else
			{
				for(int j=0;j<3;j++){
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		
		cudaFree((void*)dev_father);
		cudaFree((void*)dev_mother);
		cudaFree((void*)dev_lk);
		cudaFree((void*)dev_pcp2);
		//cudaFree((void*)dev_postProb);
		cudaFree((void*)dev_genoFry);
		cudaFree((void*)dev_postProb);
		//cudaFree((void*)dev_postProbS);

		delete host_postProb;

	}
	else if(chrType==1)
	{
		vector<double> genoProbM;
		vector<double> genoProbF;
		if(Known)
		{
			genoProbM=genoProbXK;
		}
		else
		{
			genoProbM=genoProbXN;
		}
		if(Known)
		{
			genoProbF=genoProbK;
		}
		else
		{
			genoProbF=genoProbN;
		}

		/*
		set up data that will be used in BN
		*/
		vector<int> host_gender(numInd);
		for(size_t i=0;i<numInd;i++){
			host_gender[i] = member[i].get_gender();
		}

		vector<int> host_father(numInd,-1);
		vector<int> host_mother(numInd,-1);
		for(size_t i=0;i<numInd;i++){
			if(parent[i].size()!=0){
				host_mother[i]=parent[i][0];
				host_father[i]=parent[i][1];
			}
		}

		vector<double> host_lk(numInd*3,0);
		for(size_t i=0;i<numInd;i++){
			host_lk[i*3] = likelihood(i,0);
			host_lk[i*3+1] = likelihood(i,1);
			host_lk[i*3+2] = likelihood(i,2);
		}

		vector<double> host_pcp2M(27,0);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					host_pcp2M[i*9+j*3+k] = pcp2Xm[i](j,k);
				}
			}
		}

		vector<double> host_pcp2F(27,0);
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				for(int k=0;k<3;k++){
					host_pcp2F[i*9+j*3+k] = pcp2Xf[i](j,k);
				}
			}
		}

		int numGeno = (int)pow(3.0,(double)numInd);

		cudaError_t cudaStatus;

		int * dev_gender = 0;
		int * dev_father = 0;
		int * dev_mother = 0;
		double * dev_lk = 0;
		//male
		double * dev_pcp2M =0;
		//female
		double * dev_pcp2F = 0;
		//double * dev_posteriorProb = 0;
		//male
		double * dev_genoFryM = 0;
		//female
		double * dev_genoFryF = 0;
		double * dev_postProb = 0;

		
		// Choose which GPU to run on, change this on a multi-GPU system.
		cudaStatus = cudaSetDevice(0);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		// Allocate GPU buffers
		cudaStatus = cudaMalloc((void**)&dev_gender, host_gender.size() * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_father, host_father.size() * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_mother, host_mother.size() * sizeof(int));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_lk, host_lk.size() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_pcp2M, host_pcp2M.size() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);


			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_pcp2F, host_pcp2F.size() * sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_postProb, (NUM_THREAD*NUM_BLOCK)*3*numInd*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_genoFryF, 3*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMalloc((void**)&dev_genoFryM, 3*sizeof(double));
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMalloc failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}


		/*
		copy data
		*/
		cudaStatus = cudaMemcpy(dev_gender, &host_gender[0], host_gender.size() * sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMemcpy(dev_father, &host_father[0], host_father.size() * sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}
		
		/*
		int *test_father = new int[host_father.size()];
		cudaMemcpy(test_father, dev_father, host_father.size()*sizeof(int), cudaMemcpyDeviceToHost);
		for(size_t i =0; i<host_father.size();i++){
			cout<<host_father[i]<<'\t'<<test_father[i]<<endl;
		}
		delete test_father;
		*/
		
		

		cudaStatus = cudaMemcpy(dev_mother, &host_mother[0], host_mother.size() * sizeof(int), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		/*
		int *test_mother = new int[host_mother.size()];
		cudaMemcpy(test_mother, dev_mother, host_mother.size()*sizeof(int), cudaMemcpyDeviceToHost);
		for(size_t i =0; i<host_mother.size();i++){
			cout<<host_mother[i]<<'\t'<<test_mother[i]<<endl;
		}
		delete test_mother;
		*/
		

		cudaStatus = cudaMemcpy(dev_lk, &host_lk[0], host_lk.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		/*
		double* test_lk = new double[host_lk.size()];
		cudaMemcpy(test_lk, dev_lk, host_lk.size()*sizeof(double), cudaMemcpyDeviceToHost);
		for(size_t i=0;i<host_lk.size();i++){
			cout<<test_lk[i]<<'\t'<<host_lk[i]<<endl;
		}
		delete test_lk;
		*/

		cudaStatus = cudaMemcpy(dev_pcp2M, &host_pcp2M[0], host_pcp2M.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMemcpy(dev_pcp2F, &host_pcp2F[0], host_pcp2F.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		/*
		double* test_pcp2M = new double[host_pcp2M.size()];
		cudaMemcpy(test_pcp2M, dev_pcp2M, host_pcp2M.size()*sizeof(double),cudaMemcpyDeviceToHost);
		for(size_t i=0;i<host_pcp2M.size();i++)
		{
			cout<<test_pcp2M[i]<<'\t'<<host_pcp2M[i]<<endl;
		}
		delete test_pcp2M;

		double* test_pcp2F = new double[host_pcp2F.size()];
		cudaMemcpy(test_pcp2F, dev_pcp2F, host_pcp2F.size()*sizeof(double),cudaMemcpyDeviceToHost);
		for(size_t i=0;i<host_pcp2F.size();i++)
		{
			cout<<test_pcp2F[i]<<'\t'<<host_pcp2F[i]<<endl;
		}
		delete test_pcp2F;
		*/
		

		cudaStatus = cudaMemcpy(dev_genoFryF, &genoProbF[0], genoProbF.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}

		cudaStatus = cudaMemcpy(dev_genoFryM, &genoProbM[0], genoProbM.size() * sizeof(double), cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			flagPB=false;
			return false;
		}


/*
		// Allocate CUDA events that we'll use for timing
		cudaError_t error;
		cudaEvent_t start;
		error = cudaEventCreate(&start);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		cudaEvent_t stop;
		error = cudaEventCreate(&stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(start, NULL);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		calPostProbX<<<NUM_BLOCK, NUM_THREAD>>>(dev_gender, dev_father, dev_mother, dev_lk, dev_pcp2M, dev_pcp2F, dev_postProb, dev_genoFryM, dev_genoFryF, numInd, numGeno);
		
				
		// Record the stop event
		error = cudaEventRecord(stop, NULL);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Wait for the stop event to complete
		error = cudaEventSynchronize(stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		float msecTotal = 0.0f;
		error = cudaEventElapsedTime(&msecTotal, start, stop);

		if (error != cudaSuccess)
		{
			fprintf(stderr, "Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}

		// Compute and print the performance
		printf("Time= %.3f msec\n",msecTotal);
*/		
		//sumPostProb<<<1,1>>>(dev_postProb, dev_posteriorProb, numGeno, numInd);
		
		calPostProbX<<<NUM_BLOCK, NUM_THREAD>>>(dev_gender, dev_father, dev_mother, dev_lk, dev_pcp2M, dev_pcp2F, dev_postProb, dev_genoFryM, dev_genoFryF, numInd, numGeno);

		double * host_postProb = new double[NUM_THREAD*NUM_BLOCK*3*numInd];
		cudaStatus = cudaMemcpy(host_postProb, dev_postProb, NUM_THREAD*NUM_BLOCK*3*numInd*sizeof(double), cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) {
			fprintf(stderr, "cudaMemcpy failed!");
			cudaFree((void*)dev_gender);
			cudaFree((void*)dev_father);
			cudaFree((void*)dev_mother);
			cudaFree((void*)dev_lk);
			cudaFree((void*)dev_pcp2M);
			cudaFree((void*)dev_pcp2F);
			cudaFree((void*)dev_postProb);
			cudaFree((void*)dev_genoFryM);
			cudaFree((void*)dev_genoFryF);

			delete host_postProb;

			flagPB=false;
			return false;
		}
		
		postProb.setZero();
		for(int i=0;i<NUM_BLOCK*NUM_THREAD;i++){
			for(int j=0;j<numInd;j++){
				for(int k=0;k<3;k++){
					postProb(j,k)=postProb(j,k)+host_postProb[i*numInd*3+j*3+k];
				}
			}
		}
		
		for(unsigned int i=0;i<numInd;i++){
			double sTmp = 0;
			for(int j=0;j<3;j++){
				sTmp = sTmp + postProb(i,j);
				//cout<<postProb(i,j)<<'\t';
			}
			//cout<<endl;
			if(sTmp<=0)
			{
				return false;
			}
			else
			{
				for(int j=0;j<3;j++){
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}

		cudaFree((void*)dev_gender);
		cudaFree((void*)dev_father);
		cudaFree((void*)dev_mother);
		cudaFree((void*)dev_lk);
		cudaFree((void*)dev_pcp2M);
		cudaFree((void*)dev_pcp2F);
		cudaFree((void*)dev_postProb);
		cudaFree((void*)dev_genoFryM);
		cudaFree((void*)dev_genoFryF);

		delete host_postProb;
	}

	flagPB=true;
	return true;
}
#endif

#if 0
bool family::calPostProbBN(bool Known, int chrType)
{
	if(!flagLK)
	{
		cout<<"Likelihood has not been set. Please set likelihood first."<<endl;
		return false;
	}

	if(!calPostProbSingle(Known,chrType))
	{
		//cout<<Known<<'\t'<<chrType<<endl;
		//cout<<likelihood<<endl;
		return false;
	}

	//cout<<likelihood<<endl;

	bool flagSingle=true;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double big=0;
			double sum=0;
			for(int j=0;j<3;j++)
			{
				if(big<likelihood(mapV2P[i],j))
				{
					big=likelihood(mapV2P[i],j);
				}
				sum=sum+likelihood(mapV2P[i],j);
			}
			big=big/sum;
			if(big<m_lc)
			{
				flagSingle=false;
				break;
			}
		}
	}

	postProb.setZero();

	if(flagSingle)
	{
		if(chrType==0)
		{
			vector<double> genoProb;
			if(Known)
			{
				genoProb=genoProbK;
			}
			else
			{
				genoProb=genoProbN;
			}
			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					postProb(i,j)=likelihood(i,j)*genoProb[j];
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		else if(chrType==1)
		{
			vector<double> genoProbM;
			vector<double> genoProbF;
			if(Known)
			{
				genoProbM=genoProbXK;
			}
			else
			{
				genoProbM=genoProbXN;
			}
			if(Known)
			{
				genoProbF=genoProbK;
			}
			else
			{
				genoProbF=genoProbN;
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					if(member[i].get_gender()==1)
					{
						postProb(i,j)=likelihood(i,j)*genoProbM[j];
					}
					else
					{
						postProb(i,j)=likelihood(i,j)*genoProbF[j];
					}
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		flagPB=true;
		return true;
	}



	if(chrType==0)
	{
		vector<double> genoProb;
		if(Known)
		{
			genoProb=genoProbK;
		}
		else
		{
			genoProb=genoProbN;
		}

		vector<int> genotype(numInd,0);
		vector<double> indProb(numInd,0);

		while(true)
		{
			for(unsigned int i=0;i<numInd;i++)
			{
				if(parent[i].size()==0)
				{
					indProb[i]=genoProb[genotype[i]]*likelihood(i,genotype[i]);
				}
				else
				{
					indProb[i]=pcp2[genotype[i]](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,genotype[i]);
				}
			}

			double probAll=10000000;
			for(unsigned int i=0;i<numInd;i++)
			{
				//cout<<indProb[i]<<endl;
				probAll=probAll*indProb[i];
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				postProb(i,genotype[i])=postProb(i,genotype[i])+probAll;
			}

			unsigned int index=0;
			while(index<numInd)
			{
				genotype[index]=genotype[index]+1;
				if(genotype[index]==3)
				{
					genotype[index]=0;
					index++;
				}
				else
				{
					break;
				}
			}
			if(index==numInd)
			{
				break;
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			double sTmp=postProb.sum_row(i);
			if(sTmp<=0)
			{
				return false;
			}
			for(unsigned int j=0;j<3;j++)
			{
				postProb(i,j)=postProb(i,j)/sTmp;
			}
		}

		///////////////
//		cout<<postProb<<endl;
//		cout<<postProbSingle<<endl;
//		for(size_t ii=0;ii<child.size();ii++)
//		{
//			for(size_t jj=0;jj<child[ii].size();jj++)
//			{
//				cout<<child[ii][jj]<<'\t';
//			}
//			cout<<endl;
//		}
//		cout<<endl;
//
//		for(size_t ii=0;ii<parent.size();ii++)
//		{
//			for(size_t jj=0;jj<parent[ii].size();jj++)
//			{
//				cout<<parent[ii][jj]<<'\t';
//			}
//			cout<<endl;
//		}
//		cout<<endl;
//
//		for(size_t ii=0;ii<spouse.size();ii++)
//		{
//			for(size_t jj=0;jj<spouse[ii].size();jj++)
//			{
//				cout<<spouse[ii][jj]<<'\t';
//			}
//			cout<<endl;
//		}
//		cout<<endl;

	}
	else if(chrType==1)
	{
		vector<double> genoProbM;
		vector<double> genoProbF;
		if(Known)
		{
			genoProbM=genoProbXK;
		}
		else
		{
			genoProbM=genoProbXN;
		}
		if(Known)
		{
			genoProbF=genoProbK;
		}
		else
		{
			genoProbF=genoProbN;
		}

		vector<int> genotype(numInd,0);
		vector<double> indProb(numInd,0);

		while(true)
		{
			/*
			bool flagConti=false;
			for(unsigned int i=0;i<numInd;i++)
			{
				if(member[i].get_gender()==1 && genotype[i]==1)
				{
					flagConti=true;
					break;
				}
			}
			if(flagConti)
			{
				unsigned int index=0;
				while(index<numInd)
				{
					genotype[index]=genotype[index]+1;
					if(genotype[index]==3)
					{
						genotype[index]=0;
						index++;
					}
					else
					{
						break;
					}
				}
				if(index==numInd)
				{
					break;
				}
				continue;
			}
			*/

			for(unsigned int i=0;i<numInd;i++)
			{
				if(parent[i].size()==0)
				{
					if(member[i].get_gender()==1)
					{
						indProb[i]=genoProbM[genotype[i]]*likelihood(i,genotype[i]);
					}
					else
					{
						indProb[i]=genoProbF[genotype[i]]*likelihood(i,genotype[i]);
					}
				}
				else
				{
					if(member[i].get_gender()==1)
					{
						indProb[i]=pcp2Xm[genotype[i]](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,genotype[i]);
					}
					else
					{
						indProb[i]=pcp2Xf[genotype[i]](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,genotype[i]);
					}
				}
			}

			double probAll=10000000;
			for(unsigned int i=0;i<numInd;i++)
			{
				//cout<<indProb[i]<<endl;
				probAll=probAll*indProb[i];
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				postProb(i,genotype[i])=postProb(i,genotype[i])+probAll;
			}

			unsigned int index=0;
			while(index<numInd)
			{
				genotype[index]=genotype[index]+1;
				if(genotype[index]==3)
				{
					genotype[index]=0;
					index++;
				}
				else
				{
					break;
				}
			}
			if(index==numInd)
			{
				break;
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			double sTmp=postProb.sum_row(i);
			if(sTmp<=0)
			{
				return false;
			}
			for(unsigned int j=0;j<3;j++)
			{
				postProb(i,j)=postProb(i,j)/sTmp;
			}
		}
	}

	flagPB=true;
	return true;
}

#endif

bool family::calPostProbPeeling(bool Known, int chrType)
{
	if(!flagLK)
	{
		cout<<"Likelihood has not been set. Please set likelihood first."<<endl;
		return false;
	}

	if(!calPostProbSingle(Known,chrType))
	{
		return false;
	}

	//calPostProbSingle();
	bool flagSingle=true;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double big=0;
			double sum=0;
			for(int j=0;j<3;j++)
			{
				if(big<likelihood(mapV2P[i],j))
				{
					big=likelihood(mapV2P[i],j);
				}
				sum=sum+likelihood(mapV2P[i],j);
			}
			big=big/sum;
			if(big<m_lc)
			{
				flagSingle=false;
				break;
			}
		}
	}

	postProb.setZero();

	//flagSingle=false;

	if(flagSingle)
	{
		if(chrType==0)
		{
			vector<double> genoProb;
			if(Known)
			{
				genoProb=genoProbK;
			}
			else
			{
				genoProb=genoProbN;
			}
			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					postProb(i,j)=likelihood(i,j)*genoProb[j];
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		else if(chrType==1)
		{
			vector<double> genoProbM;
			vector<double> genoProbF;
			if(Known)
			{
				genoProbM=genoProbXK;
			}
			else
			{
				genoProbM=genoProbXN;
			}
			if(Known)
			{
				genoProbF=genoProbK;
			}
			else
			{
				genoProbF=genoProbN;
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					if(member[i].get_gender()==1)
					{
						postProb(i,j)=likelihood(i,j)*genoProbM[j];
					}
					else
					{
						postProb(i,j)=likelihood(i,j)*genoProbF[j];
					}
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		flagPB=true;
		return true;
	}

	if(chrType==0)
	{
		vector<double> genoProb;
		if(Known)
		{
			genoProb=genoProbK;
		}
		else
		{
			genoProb=genoProbN;
		}

		dMatrix<double> antProb(numInd,3,-1);
		//posProb[iInd](jInd,iGeno)
		//iInd: index of i
		//jInd: index of j (i's spouse)
		//iGeno: i's genotype
		vector<dMatrix<double> > posProb;
		for(unsigned int i=0;i<numInd;i++)
		{
			dMatrix<double> dMTmp(numInd,3,-1);
			posProb.push_back(dMTmp);
		}

		//prepare anterior probability
		//set founder's anterior probability
		for(unsigned int i=0;i<numInd;i++)
		{
			if(parent[i].size()==0)
			{
				for(unsigned int j=0;j<3;j++)
				{
					antProb(i,j)=genoProb[j];
				}
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			for(unsigned int j=0;j<3;j++)
			{
				double postTmp=1;
				for(unsigned int k=0;k<spouse[i].size();k++)
				{
					postTmp=postTmp*calPosProb(i,j,spouse[i][k],antProb,posProb);
				}
				postTmp=postTmp*likelihood(i,j)*calAntProb(i,j,antProb,posProb);
				postProb(i,j)=postTmp;
			}

			double sTmp=postProb.sum_row(i);
			if(sTmp==0)
			{
				//cout<<antProb<<endl;
				return false;
			}
			for(unsigned int j=0;j<3;j++)
			{
				postProb(i,j)=postProb(i,j)/sTmp;
			}
		}
	}
	else if(chrType==1)
	{
		vector<double> genoProbM;
		vector<double> genoProbF;
		if(Known)
		{
			genoProbM=genoProbXK;
		}
		else
		{
			genoProbM=genoProbXN;
		}
		if(Known)
		{
			genoProbF=genoProbK;
		}
		else
		{
			genoProbF=genoProbN;
		}

		dMatrix<double> antProb(numInd,3,-1);
		//posProb[iInd](jInd,iGeno)
		//iInd: index of i
		//jInd: index of j (i's spouse)
		//iGeno: i's genotype
		vector<dMatrix<double> > posProb;
		for(unsigned int i=0;i<numInd;i++)
		{
			dMatrix<double> dMTmp(numInd,3,-1);
			posProb.push_back(dMTmp);
		}

		//prepare anterior probability
		//set founder's anterior probability
		for(unsigned int i=0;i<numInd;i++)
		{
			if(parent[i].size()==0)
			{
				if(member[i].get_gender()==1)
				{
					for(unsigned int j=0;j<3;j++)
					{
						antProb(i,j)=genoProbM[j];
					}
				}
				else
				{
					for(unsigned int j=0;j<3;j++)
					{
						antProb(i,j)=genoProbF[j];
					}
				}
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			for(int j=0;j<3;j++)
			{
				double postTmp=1;
				for(unsigned int k=0;k<spouse[i].size();k++)
				{
					postTmp=postTmp*calPosProbX(i,j,spouse[i][k],antProb,posProb);
				}
				postTmp=postTmp*likelihood(i,j)*calAntProbX(i,j,antProb,posProb);
				postProb(i,j)=postTmp;
			}

			double sTmp=postProb.sum_row(i);
			if(sTmp==0)
			{
				//cout<<antProb<<endl;
				return false;
			}
			for(int j=0;j<3;j++)
			{
				postProb(i,j)=postProb(i,j)/sTmp;
			}
		}
	}
	else
	{
	}
	flagPB=true;
	return true;
}

bool family::calPostProbSingle(bool Known, int chrType)
{
	if(!flagLK)
	{
		cout<<"Likelihood has not been set. Please set likelihood first."<<endl;
		return false;
	}

	postProbSingle.setZero();

	if(chrType==0)
	{
		vector<double> genoProb;
		if(Known)
		{
			genoProb=genoProbK;
		}
		else
		{
			genoProb=genoProbN;
		}
		for(unsigned int i=0;i<numInd;i++)
		{
			for(int j=0;j<3;j++)
			{
				postProbSingle(i,j)=likelihood(i,j)*genoProb[j];
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			double sTmp=postProbSingle.sum_row(i);
			if(sTmp<=0)
			{
				return false;
			}
			for(unsigned int j=0;j<3;j++)
			{
				postProbSingle(i,j)=postProbSingle(i,j)/sTmp;
			}
		}
	}
	else if(chrType==1)
	{
		vector<double> genoProbM;
		vector<double> genoProbF;
		if(Known)
		{
			genoProbM=genoProbXK;
		}
		else
		{
			genoProbM=genoProbXN;
		}
		if(Known)
		{
			genoProbF=genoProbK;
		}
		else
		{
			genoProbF=genoProbN;
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			for(int j=0;j<3;j++)
			{
				if(member[i].get_gender()==1)
				{
					postProbSingle(i,j)=likelihood(i,j)*genoProbM[j];
				}
				else
				{
					postProbSingle(i,j)=likelihood(i,j)*genoProbF[j];
				}
			}
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			double sTmp=postProbSingle.sum_row(i);
			if(sTmp<=0)
			{
				return false;
			}
			for(unsigned int j=0;j<3;j++)
			{
				postProbSingle(i,j)=postProbSingle(i,j)/sTmp;
			}
		}
	}

	flagPBS=true;
	return true;
}

double family::calAntProb(int iInd, int iGeno, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb)
{
	//double returnVal=0;

	if(antProb(iInd,iGeno)>=0)
	{
		return antProb(iInd,iGeno);
	}

	//unsigned int numGeno=priorProb.get_column();

	/*
	if(parent[iInd].size()==1)
	{
		int p1=parent[iInd][0];

		//other child beside iInd
		vector<int> otherChild;
		for(unsigned int i=0;i<child[p1].size();i++)
		{
			if(child[p1][i]!=iInd && parent[child[p1][i]].size()==1)
			{
				otherChild.push_back(child[p1][i]);
			}
		}

		//other spouse of p1
		vector<int> otherSpouse1=spouse[p1];

		double sm=0;
		for(unsigned int i=0;i<3;i++)
		{
			double sf=0;
			for(unsigned int j=0;j<numGeno;j++)
			{
				double mc=1;
				for(unsigned int k=0;k<otherChild.size();k++)
				{
					double sc=0;
					for(unsigned int l=0;l<numGeno;l++)
					{
						double mcs=1;
						for(unsigned int m=0;m<spouse[otherChild[k]].size();m++)
						{
							mcs=mcs*calPosProb(otherChild[k],l,spouse[otherChild[k]][m],antProb,posProb);
						}
						sc=sc+mcs*likelihood(otherChild[k],l)*pcp2[l](i,j);
					}
					mc=mc*sc;
				}

				//no other spouse for p2
				double mf=1;


				sf=sf+genoProb[j]*mf*pcp2[iGeno](i,j)*mc/numGeno;
			}
			double mm=1;
			for(unsigned int j=0;j<otherSpouse1.size();j++)
			{
				mm=mm*calPosProb(p1,i,otherSpouse1[j],antProb,posProb);
			}
			sm=sm+calAntProb(p1,i,antProb,posProb)*likelihood(p1,i)*mm*sf;
		}

		antProb(iInd,iGeno)=sm;
		return sm;
	}
	else
	*/
	{
		int p1=parent[iInd][0];
		int p2=parent[iInd][1];

		//other child beside iInd
		vector<int> otherChild;
		for(unsigned int i=0;i<child[p1].size();i++)
		{
			for(unsigned int j=0;j<child[p2].size();j++)
			{
				if(child[p1][i]==child[p2][j] && child[p1][i]!=iInd)
				{
					otherChild.push_back(child[p1][i]);
					break;
				}
			}
		}

		//p1's other spouse beside p2
		vector<int> otherSpouse1;
		for(unsigned int i=0;i<spouse[p1].size();i++)
		{
			if(spouse[p1][i]!=p2)
			{
				otherSpouse1.push_back(spouse[p1][i]);
			}
		}

		//p2's other spouse beside p1
		vector<int> otherSpouse2;
		for(unsigned int i=0;i<spouse[p2].size();i++)
		{
			if(spouse[p2][i]!=p1)
			{
				otherSpouse2.push_back(spouse[p2][i]);
			}
		}

		double sm=0;
		for(unsigned int i=0;i<3;i++)
		{
			double sf=0;
			for(unsigned int j=0;j<3;j++)
			{
				double mc=1;
				for(unsigned int k=0;k<otherChild.size();k++)
				{
					double sc=0;
					for(unsigned int l=0;l<3;l++)
					{
						double mcs=1;
						for(unsigned int m=0;m<spouse[otherChild[k]].size();m++)
						{
							mcs=mcs*calPosProb(otherChild[k],l,spouse[otherChild[k]][m],antProb,posProb);
						}
						sc=sc+mcs*likelihood(otherChild[k],l)*pcp2[l](i,j);
					}
					mc=mc*sc;
				}

				double mf=1;
				for(unsigned int k=0;k<otherSpouse2.size();k++)
				{
					mf=mf*calPosProb(p2,j,otherSpouse2[k],antProb,posProb);
				}
				sf=sf+calAntProb(p2,j,antProb,posProb)*likelihood(p2,j)*mf*pcp2[iGeno](i,j)*mc;
			}
			double mm=1;
			for(unsigned int j=0;j<otherSpouse1.size();j++)
			{
				mm=mm*calPosProb(p1,i,otherSpouse1[j],antProb,posProb);
			}
			sm=sm+calAntProb(p1,i,antProb,posProb)*likelihood(p1,i)*mm*sf;
		}

		antProb(iInd,iGeno)=sm;
		return sm;
	}
}

double family::calAntProbX(int iInd, int iGeno, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb)
{
	//double returnVal=0;

	if(antProb(iInd,iGeno)>=0)
	{
		return antProb(iInd,iGeno);
	}

	//unsigned int numGeno=priorProb.get_column();
	{
		int p1=parent[iInd][0];
		int p2=parent[iInd][1];

		//other child beside iInd
		vector<int> otherChild;
		for(unsigned int i=0;i<child[p1].size();i++)
		{
			for(unsigned int j=0;j<child[p2].size();j++)
			{
				if(child[p1][i]==child[p2][j] && child[p1][i]!=iInd)
				{
					otherChild.push_back(child[p1][i]);
					break;
				}
			}
		}

		//p1's other spouse beside p2
		vector<int> otherSpouse1;
		for(unsigned int i=0;i<spouse[p1].size();i++)
		{
			if(spouse[p1][i]!=p2)
			{
				otherSpouse1.push_back(spouse[p1][i]);
			}
		}

		//p2's other spouse beside p1
		vector<int> otherSpouse2;
		for(unsigned int i=0;i<spouse[p2].size();i++)
		{
			if(spouse[p2][i]!=p1)
			{
				otherSpouse2.push_back(spouse[p2][i]);
			}
		}

		double sm=0;
		for(unsigned int i=0;i<3;i++)
		{
			double sf=0;
			for(unsigned int j=0;j<3;j++)
			{
				double mc=1;
				for(unsigned int k=0;k<otherChild.size();k++)
				{
					double sc=0;
					for(unsigned int l=0;l<3;l++)
					{
						double mcs=1;
						for(unsigned int m=0;m<spouse[otherChild[k]].size();m++)
						{
							mcs=mcs*calPosProbX(otherChild[k],l,spouse[otherChild[k]][m],antProb,posProb);
						}
						if(member[otherChild[k]].get_gender()==1)
						{
							if(member[p1].get_gender()==1)
							{
								sc=sc+mcs*likelihood(otherChild[k],l)*pcp2Xm[l](j,i);
							}
							else
							{
								sc=sc+mcs*likelihood(otherChild[k],l)*pcp2Xm[l](i,j);
							}
						}
						else
						{
							if(member[p1].get_gender()==1)
							{
								sc=sc+mcs*likelihood(otherChild[k],l)*pcp2Xf[l](j,i);
							}
							else
							{
								sc=sc+mcs*likelihood(otherChild[k],l)*pcp2Xf[l](i,j);
							}
						}
					}
					mc=mc*sc;
				}

				double mf=1;
				for(unsigned int k=0;k<otherSpouse2.size();k++)
				{
					mf=mf*calPosProbX(p2,j,otherSpouse2[k],antProb,posProb);
				}
				if(member[iInd].get_gender()==1)
				{
					if(member[p1].get_gender()==1)
					{
						sf=sf+calAntProbX(p2,j,antProb,posProb)*likelihood(p2,j)*mf*pcp2Xm[iGeno](j,i)*mc;
					}
					else
					{
						sf=sf+calAntProbX(p2,j,antProb,posProb)*likelihood(p2,j)*mf*pcp2Xm[iGeno](i,j)*mc;
					}
				}
				else
				{
					if(member[p1].get_gender()==1)
					{
						sf=sf+calAntProbX(p2,j,antProb,posProb)*likelihood(p2,j)*mf*pcp2Xf[iGeno](j,i)*mc;
					}
					else
					{
						sf=sf+calAntProbX(p2,j,antProb,posProb)*likelihood(p2,j)*mf*pcp2Xf[iGeno](i,j)*mc;
					}
				}
			}
			double mm=1;
			for(unsigned int j=0;j<otherSpouse1.size();j++)
			{
				mm=mm*calPosProbX(p1,i,otherSpouse1[j],antProb,posProb);
			}
			sm=sm+calAntProbX(p1,i,antProb,posProb)*likelihood(p1,i)*mm*sf;
		}

		antProb(iInd,iGeno)=sm;
		return sm;
	}
}

double family::calPosProb(int iInd, int iGeno, int jInd, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb)
{
	//double returnVal=0;

	if(posProb[iInd](jInd,iGeno)>=0)
	{
		return posProb[iInd](jInd,iGeno);
	}

	//unsigned int numGeno=priorProb.get_column();

	//jInd's other spouse beside iInd
	vector<int> otherSpouse;
	for(unsigned int i=0;i<spouse[jInd].size();i++)
	{
		if(spouse[jInd][i]!=iInd)
		{
			otherSpouse.push_back(spouse[jInd][i]);
		}
	}

	//iInd and jInd's child
	vector<int> allChild;
	for(unsigned int i=0;i<child[iInd].size();i++)
	{
		for(unsigned int j=0;j<child[jInd].size();j++)
		{
			if(child[iInd][i]==child[jInd][j])
			{
				allChild.push_back(child[iInd][i]);
			}
		}
	}

	double sj=0;
	for(unsigned int i=0;i<3;i++)
	{
		double ms=1;
		for(unsigned int j=0;j<otherSpouse.size();j++)
		{
			ms=ms*calPosProb(jInd,i,otherSpouse[j],antProb,posProb);
		}
		double mc=1;
		for(unsigned int j=0;j<allChild.size();j++)
		{
			double sc=0;
			for(unsigned int k=0;k<3;k++)
			{
				double mcs=1;
				for(unsigned int l=0;l<spouse[allChild[j]].size();l++)
				{
					mcs=mcs*calPosProb(allChild[j],k,spouse[allChild[j]][l],antProb,posProb);
				}
				sc=sc+pcp2[k](iGeno,i)*likelihood(allChild[j],k)*mcs;
			}
			mc=mc*sc;
		}
		sj=sj+calAntProb(jInd,i,antProb,posProb)*likelihood(jInd,i)*ms*mc;
	}

	posProb[iInd](jInd,iGeno)=sj;
	return sj;
}

double family::calPosProbX(int iInd, int iGeno, int jInd, dMatrix<double> & antProb, vector<dMatrix<double> > & posProb)
{
	//double returnVal=0;

	if(posProb[iInd](jInd,iGeno)>=0)
	{
		return posProb[iInd](jInd,iGeno);
	}

	//unsigned int numGeno=priorProb.get_column();

	//jInd's other spouse beside iInd
	vector<int> otherSpouse;
	for(unsigned int i=0;i<spouse[jInd].size();i++)
	{
		if(spouse[jInd][i]!=iInd)
		{
			otherSpouse.push_back(spouse[jInd][i]);
		}
	}

	//iInd and jInd's child
	vector<int> allChild;
	for(unsigned int i=0;i<child[iInd].size();i++)
	{
		for(unsigned int j=0;j<child[jInd].size();j++)
		{
			if(child[iInd][i]==child[jInd][j])
			{
				allChild.push_back(child[iInd][i]);
			}
		}
	}

	double sj=0;
	for(unsigned int i=0;i<3;i++)
	{
		double ms=1;
		for(unsigned int j=0;j<otherSpouse.size();j++)
		{
			ms=ms*calPosProbX(jInd,i,otherSpouse[j],antProb,posProb);
		}
		double mc=1;
		for(unsigned int j=0;j<allChild.size();j++)
		{
			double sc=0;
			for(unsigned int k=0;k<3;k++)
			{
				double mcs=1;
				for(unsigned int l=0;l<spouse[allChild[j]].size();l++)
				{
					mcs=mcs*calPosProbX(allChild[j],k,spouse[allChild[j]][l],antProb,posProb);
				}
				if(member[allChild[j]].get_gender()==1)
				{
					if(member[iInd].get_gender()==1)
					{
						sc=sc+pcp2Xm[k](i,iGeno)*likelihood(allChild[j],k)*mcs;
					}
					else
					{
						sc=sc+pcp2Xm[k](iGeno,i)*likelihood(allChild[j],k)*mcs;
					}
				}
				else
				{
					if(member[iInd].get_gender()==1)
					{
						sc=sc+pcp2Xf[k](i,iGeno)*likelihood(allChild[j],k)*mcs;
					}
					else
					{
						sc=sc+pcp2Xf[k](iGeno,i)*likelihood(allChild[j],k)*mcs;
					}
				}
			}
			mc=mc*sc;
		}
		sj=sj+calAntProbX(jInd,i,antProb,posProb)*likelihood(jInd,i)*ms*mc;
	}

	posProb[iInd](jInd,iGeno)=sj;
	return sj;
}

bool family::calPostProbMCMC(int numBurnIn, int numRep,bool Known,int chrType)
{
	if(!flagLK)
	{
		cout<<"Likelihood has not been set. Please set likelihood first."<<endl;
		return false;
	}

	if(!calPostProbSingle(Known,chrType))
	{
		//cout<<Known<<'\t'<<chrType<<endl;
		//cout<<likelihood<<endl;
		return false;
	}

	//cout<<likelihood<<endl;

	bool flagSingle=true;
	for(size_t i=0;i<mapV2P.size();i++)
	{
		if(mapV2P[i]>=0)
		{
			double big=0;
			double sum=0;
			for(int j=0;j<3;j++)
			{
				if(big<likelihood(mapV2P[i],j))
				{
					big=likelihood(mapV2P[i],j);
				}
				sum=sum+likelihood(mapV2P[i],j);
			}
			big=big/sum;
			if(big<m_lc)
			{
				flagSingle=false;
				break;
			}
		}
	}

	postProb.setZero();

	if(flagSingle)
	{
		if(chrType==0)
		{
			vector<double> genoProb;
			if(Known)
			{
				genoProb=genoProbK;
			}
			else
			{
				genoProb=genoProbN;
			}
			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					postProb(i,j)=likelihood(i,j)*genoProb[j];
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		else if(chrType==1)
		{
			vector<double> genoProbM;
			vector<double> genoProbF;
			if(Known)
			{
				genoProbM=genoProbXK;
			}
			else
			{
				genoProbM=genoProbXN;
			}
			if(Known)
			{
				genoProbF=genoProbK;
			}
			else
			{
				genoProbF=genoProbN;
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				for(int j=0;j<3;j++)
				{
					if(member[i].get_gender()==1)
					{
						postProb(i,j)=likelihood(i,j)*genoProbM[j];
					}
					else
					{
						postProb(i,j)=likelihood(i,j)*genoProbF[j];
					}
				}
			}

			for(unsigned int i=0;i<numInd;i++)
			{
				double sTmp=postProb.sum_row(i);
				if(sTmp<=0)
				{
					return false;
				}
				for(unsigned int j=0;j<3;j++)
				{
					postProb(i,j)=postProb(i,j)/sTmp;
				}
			}
		}
		flagPB=true;
		return true;
	}


	vector<int> genotype(numInd,0);
	for(unsigned int i=0;i<numInd;i++)
	{
		genotype[i]=rand()%3;
	}

	dMatrix<double> genoFry(numInd,3,0);
	for(int i=0;i<numBurnIn;i++)
	{
		estGenoProb(genotype,genoFry,Known,chrType);
	}

	genoFry.setZero();

	for(int i=0;i<numRep;i++)
	{
		estGenoProb(genotype,genoFry,Known,chrType);
	}

	for(unsigned int i=0;i<numInd;i++)
	{
		for(int j=0;j<3;j++)
		{
			postProb(i,j)=genoFry(i,j)/numRep;
		}
		if(postProb.sum_row(i)<=0)
		{
			return false;
		}
	}

	flagPB=true;
	return true;
}

bool family::estGenoProb(std::vector<int> & genotype, dMatrix<double> & genoFre,bool Known,int chrType)
{
	//normal chromosome
	if(chrType==0)
	{
		vector<double> genoProb;
		if(Known)
		{
			genoProb=genoProbK;
		}
		else
		{
			genoProb=genoProbN;
		}
		//cout<<likelihood<<endl;
		for(unsigned int i=0;i<numInd;i++)
		{
			vector<double> genoFreTmp(3,1000000);
			for(int j=0;j<3;j++)
			{
				//cout<<genoFreTmp[j]<<endl;
				if(parent[i].size()==0)
				{
					genoFreTmp[j]=genoFreTmp[j]*genoProb[j]*likelihood(i,j);
				}
				else
				{
					genoFreTmp[j]=genoFreTmp[j]*pcp2[j](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,j);
				}
				//cout<<genoFreTmp[j]<<endl;
				for(size_t k=0;k<child[i].size();k++)
				{
					int indChild=child[i][k];
					if(member[i].get_gender()==1)
					{
						genoFreTmp[j]=genoFreTmp[j]*pcp2[genotype[indChild]](genotype[parent[indChild][0]],j);
					}
					else
					{
						genoFreTmp[j]=genoFreTmp[j]*pcp2[genotype[indChild]](j,genotype[parent[indChild][1]]);
					}
					//cout<<genoFreTmp[j]<<endl;
				}
			}
			double sFre=0;
			for(int j=0;j<3;j++)
			{
				sFre=sFre+genoFreTmp[j];
			}
			if(sFre<=0)
			{
				for(int j=0;j<3;j++)
				{
					genoFreTmp[j]=0;
				}
			}
			else
			{
				for(int j=0;j<3;j++)
				{
					genoFreTmp[j]=genoFreTmp[j]/sFre;
				}
			}
			double rd=double(rand())/double(RAND_MAX);
			if(rd<genoFreTmp[0])
			{
				genotype[i]=0;
			}
			else if(rd>(1.0-genoFreTmp[2]))
			{
				genotype[i]=2;
			}
			else
			{
				genotype[i]=1;
			}

			for(int j=0;j<3;j++)
			{
				genoFre(i,j)=genoFre(i,j)+genoFreTmp[j];
			}
		}
	}

	//chromosome X
	else if(chrType==1)
	{
		vector<double> genoProbM;
		vector<double> genoProbF;
		if(Known)
		{
			genoProbM=genoProbXK;
			genoProbF=genoProbK;
		}
		else
		{
			genoProbM=genoProbXN;
			genoProbF=genoProbN;
		}

		for(unsigned int i=0;i<numInd;i++)
		{
			vector<double> genoFreTmp(3,1000000);
			if(member[i].get_gender()==1)
			{
				for(int j=0;j<3;j++)
				{
					if(parent[i].size()==0)
					{
						genoFreTmp[j]=genoFreTmp[j]*genoProbM[j]*likelihood(i,j);
					}
					else
					{
						genoFreTmp[j]=genoFreTmp[j]*pcp2Xm[j](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,j);
					}
					for(size_t k=0;k<child[i].size();k++)
					{
						int indChild=child[i][k];
						if(member[i].get_gender()==1)
						{
							if(member[indChild].get_gender()==1)
							{
								genoFreTmp[j]=genoFreTmp[j]*pcp2Xm[genotype[indChild]](genotype[parent[indChild][0]],j);
							}
							else
							{
								genoFreTmp[j]=genoFreTmp[j]*pcp2Xf[genotype[indChild]](genotype[parent[indChild][0]],j);
							}
						}
					}
				}
			}
			else
			{
				for(int j=0;j<3;j++)
				{
					if(parent[i].size()==0)
					{
						genoFreTmp[j]=genoFreTmp[j]*genoProbF[j]*likelihood(i,j);
					}
					else
					{
						genoFreTmp[j]=genoFreTmp[j]*pcp2Xf[j](genotype[parent[i][0]],genotype[parent[i][1]])*likelihood(i,j);
					}
					for(size_t k=0;k<child[i].size();k++)
					{
						int indChild=child[i][k];
						if(member[i].get_gender()==1)
						{
							if(member[indChild].get_gender()==1)
							{
								genoFreTmp[j]=genoFreTmp[j]*pcp2Xm[indChild](j,genotype[parent[indChild][1]]);
							}
							else
							{
								genoFreTmp[j]=genoFreTmp[j]*pcp2Xf[indChild](j,genotype[parent[indChild][1]]);
							}
						}
					}
				}
			}
			double sFre=0;
			for(int j=0;j<3;j++)
			{
				sFre=sFre+genoFreTmp[j];
			}
			if(sFre<=0)
			{
				for(int j=0;j<3;j++)
				{
					genoFreTmp[j]=0;
				}
			}
			else
			{
				for(int j=0;j<3;j++)
				{
					genoFreTmp[j]=genoFreTmp[j]/sFre;
				}
			}
			double rd=double(rand())/double(RAND_MAX);
			if(rd<genoFreTmp[0])
			{
				genotype[i]=0;
			}
			else if(rd>(1.0-genoFreTmp[2]))
			{
				genotype[i]=2;
			}
			else
			{
				genotype[i]=1;
			}

			for(int j=0;j<3;j++)
			{
				genoFre(i,j)=genoFre(i,j)+genoFreTmp[j];
			}
		}
	}
	return true;
}
