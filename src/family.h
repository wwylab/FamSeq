/*
 * family.h
 * individual and family class
 *
 *  Created on: Mar 7, 2012
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#ifndef FAMILY_H_
#define FAMILY_H_

#include <vector>
#include <string>
#include "dMatrix.h"

enum calMethod {MCMC, BN, Peeling};

class individual
{
private:
	//individual's id
	int m_id;
	//gender: 1:male, 2:female, 0: unknown
	int m_gender;
	//age
	int m_age;
	//state
	int m_state;
	//mother's id
	int m_mId;
	//father's id;
	int m_fId;
	//sample name
	std::string m_sName;

public:
	//constructor
	individual(int id, int gender, int age, int state, int mId, int fId, std::string sName);

	//get functions
	int get_id() const { return m_id;}
	int get_gender() const { return m_gender; }
	int get_age() const { return m_age; }
	int get_state() const { return m_state; }
	int get_mId() const { return m_mId; }
	int get_fId() const { return m_fId; }
	std::string get_sName() const {return m_sName; }

	//set functions
	bool set_id( int id);
	bool set_gender(int gender);
	bool set_age(int age);
	bool set_sate(int state);
	bool set_mId(int mId);
	bool set_fId(int fId);
	bool set_sName(std::string sName);
};

class family
{
private:
	//////////////////////////////////
	///     private variables      ///
	//////////////////////////////////

	//family member
	std::vector<individual> member;

	//sometimes we cannot get child's two parents' sequence
	//in this case, add one flat parent(the likelihood of RR, RA and AA is the same) into the family
	//numInd=realNumInd+numFlatInd
	//number of individuals in a family
	unsigned int numInd;

	//real family members
	unsigned int realNumInd;


	//one individual's children (store only the index of "member")
	std::vector<std::vector<int> > child;


	//one individual's parents (store only the index of "member")
	std::vector<std::vector<int> > parent;

	//one individual's spouse (store only the index of "member")
	std::vector<std::vector<int> > spouse;

	//mutation rate
	//0 as default
	double mRate;

	/////////indicators////////

	//indicate whether the family has been initiate
	bool flagInit;

	//indicate whether the likelihood has been set
	bool flagLK;

	//indicate whether postProb has been calculated
	bool flagPB;

	//indicate whether postProbS has been calculated
	bool flagPBS;

	//indicate whether there is gender information
	bool flagGender;

	// likelihood criteria. if likelihood is larger than this value, we don't use pedigree information
	// to calculate posterior probability. We use single method instead.
	double m_lc;



	//Pr(genotype) for new variant
	std::vector<double> genoProbN;

	//Pr(genotype) for known variant
	std::vector<double> genoProbK;

	//Pr(genotype) for new variant on chromosome X for male
	std::vector<double> genoProbXN;

	//Pr(genotype) for known variant on chromosome X for male
	std::vector<double> genoProbXK;


	//child's genotype probability conditional on two parents' genotype
	//pcp2m[i](j,k): probability of child has genotype i under the condition that parents genotypes are j and k
	std::vector<dMatrix<double> > pcp2;

	//female child's genotype probability conditional on two parents' genotype at chromosome X
	//pcp2mXf[i](j,k): probability of female child has genotype i given parents' genotype are j and k
	//Pr(genotype of child=i | genotytpe of mother=j, genotype of father =k)
	std::vector<dMatrix<double> > pcp2Xf;

	//male child's genotype probability conditional on two parents' genotype at chromosome X
	//pcp2mXm[i](j,k): probability of male child has genotype i given parents' genotype are j and k
	//Pr(genotype of child=i | genotytpe of mother=j, genotype of father =k)
	std::vector<dMatrix<double> > pcp2Xm;

	//posterior probability (FamSeqPro)
	dMatrix<double> postProb;

	//posterior probability (Single)
	dMatrix<double> postProbSingle;

	//likelihood
	//Pr(D|G)
	dMatrix<double> likelihood;

	//mapping from ped file to vcf file
	std::vector<int> mapP2V;

	//mapping from vcf file to ped file
	std::vector<int> mapV2P;


	//////////////////////////////////
	///     private functions      ///
	//////////////////////////////////

	//set conditional probability based on mutation rate
	//set pcp2m, pcp2mXf,pcp2mXm
	bool setPCP();

	//get the relationship in the family
	//return false if the family is not fulfilled
	bool setRelation();

	// check pedigree data
	bool checkPed();

	//calculate genotype probability conditional on two parents genotype for one gene only
	//Input:
	//nAllele: number of allels in gene
	//Output:
	//PCP2S: genotype probability of child given the genotype of two parents
	//PCP2S[i](j,k) means the probability that child has genotype i under the condition that parents genotype are j and k
	bool calPCP2S(int nAllele, std::vector<dMatrix<double> > & PCP2S);

	bool calPCP2();

	bool calPCP2Xf();

	bool calPCP2Xm();

	//calculate anterior probabilities in peeling
	//iInd: individual's index
	//iGeno: individual's genotype
	double calAntProb(int iInd, int iGeno, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb);

	//calculate posterior probabilities in peeling
	//iInd: individual's index
	//jInd: individual's index (iInd's spouse)
	//iGeno: individual (iInd) 's genotype
	double calPosProb(int iInd, int iGeno, int jInd, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb);

	//calculate anterior probabilities in peeling for chromosome X
	//iInd: individual's index
	//iGeno: individual's genotype
	double calAntProbX(int iInd, int iGeno, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb);

	//calculate posterior probabilities in peeling for chromosome X
	//iInd: individual's index
	//jInd: individual's index (iInd's spouse)
	//iGeno: individual (iInd) 's genotype
	double calPosProbX(int iInd, int iGeno, int jInd, dMatrix<double> & antProb, std::vector<dMatrix<double> > & posProb);


	bool fulfillFamily();

	bool estGenoProb(std::vector<int> & genotype, dMatrix<double> & genoFre, bool Known, int chrType);

public:
	//constructor
	//mem: member in the family
	//lk: likelihood
	//gProb: genotype probability
	//mutationRate: mutation rate
	family(const std::vector<individual> & mem, double mutationRate=1e-7);

	family();

	family & operator=(const family & fm);

	//initiation
	bool init();

	//~family();

	//////////////////////////////////
	///       get functions        ///
	//////////////////////////////////

	//member
	std::vector<individual> get_member() const{ return member; }

	//numInd
	unsigned int get_numInd() const{ return numInd; }

	//realNumInd
	unsigned int get_realNumInd() const {return realNumInd; }

	//child
	std::vector<std::vector<int> > get_child() const {return child;}

	//parent
	std::vector<std::vector<int> > get_parent() const {return parent;}

	//spouse
	std::vector<std::vector<int> > get_spouse() const {return spouse;}

	//mutation rate
	double get_mRate() const{ return mRate; }

	//posterior probability of genotype (FamSeqPro)
	dMatrix<double> get_postProb(bool flag=true) const;

	//posterior probability of genotype (Single)
	dMatrix<double> get_postProbSingle(bool flag=true) const;

	//likelihood
	dMatrix<double> get_LK() const {return likelihood; }

	//whether the family has been initiate
	bool get_flagInit() const {return flagInit;}

	//flagLK
	bool get_flagLK() const {return flagLK;}

	//flagPB
	bool get_flagPB() const {return flagPB;}

	//flagPBS
	bool get_flagPBS() const {return flagPBS;}

	//flagGender
	bool get_flagGender() const {return flagGender;}

	//genoProbN
	std::vector<double> get_genoProbN() const {return genoProbN;}

	//genoProbK
	std::vector<double> get_genoProbK() const {return genoProbK;}

	//genoProbXN
	std::vector<double> get_genoProbXN() const {return genoProbXN;}

	//genoProbXK
	std::vector<double> get_genoProbXK() const {return genoProbXK;}

	//pcp2
	std::vector<dMatrix<double> > get_pcp2() const {return pcp2;}

	//pcp2Xf
	std::vector<dMatrix<double> > get_pcp2Xf() const {return pcp2Xf;}

	//pcp2Xm
	std::vector<dMatrix<double> > get_pcp2Xm() const {return pcp2Xm;}

	//get the estimate genotype (FamSeqPro)
	//0:RR,1:RA,2:AA
	std::vector<int> get_postRlt() const;

	//get the estimate genotype (Single)
	//0:RR,1:RA,2:AA
	std::vector<int> get_postRltSingle() const;

	//mapP2V
	std::vector<int> get_mapP2V() const {return mapP2V;}

	//mapV2P
	std::vector<int> get_mapV2P() const {return mapV2P;}

	// lc
	double get_lc() const { return m_lc;}


	//////////////////////////////////
	///       set functions        ///
	//////////////////////////////////

	//set mutation rate
	bool set_mRate(double rate);

	//set genotype probability for new variant
	bool set_genoProbN(const std::vector<double> & prob);

	//set genotype probatility for known variant
	bool set_genoProbK(const std::vector<double> & prob);

	//set genotype probability for new variant on chromosome X for male
	bool set_genoProbXN(const std::vector<double> & prob);

	//set genotype probability for known variant on chromosome X for male
	bool set_genoProbXK(const std::vector<double> & prob);

	//set likelihood
	bool set_LK(const dMatrix<double> & lk);

	//set gender flag
	bool set_fg(bool fg);

	//set mapP2V
	bool set_mapP2V(std::vector<int> mP2V);

	//set mapV2P
	bool set_mapV2P(std::vector<int> mV2P);

	// set lc
	bool set_lc(double lc);


	//////////////////////////////////
	///        calculating        ///
	//////////////////////////////////

	//calculate posterior probability by peeling
	//Known: whether it is a known SNP
	//chrType
	//0: normal chromosome
	//1: X chromosome
	bool calPostProbPeeling(bool Known=false, int chrType=0);

	//calculate posterior probability by BN
	//Known: whether it is a known SNP
	//chrType
	//0: normal chromosome
	//1: X chromosome
	bool calPostProbBN(bool Known=false, int chrType=0);

	//calculate posterior probability by Single
	//Known: whether it is a known SNP
	//chrType
	//0: normal
	//1: X chromosome
	bool calPostProbSingle(bool Known=false, int chrType=0);

	//calculate posterior probability by MCMC
	//Known: whether it is a known SNP
	//chrType
	//0: normal
	//1: X chromosome
	bool calPostProbMCMC(int numBurnIn, int numRep,bool Known=false,int chrType=0);
};

#endif /* FAMILY_H_ */
