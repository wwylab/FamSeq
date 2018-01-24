/*
 * checkInput.h
 * check input of the command
 *
 *  Created on: Mar 18, 2012
 *
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#ifndef CHECKINPUT_H_
#define CHECKINPUT_H_

#include <string>
#include <vector>

class inputVCF
{
private:
	//vcf file name
	std::vector<std::string> fileNameVCF;

	//ped file name
	std::string fileNamePed;

	//output file name
	std::string outputName;

	//the file to indicate the position that we want to call the variant
	std::string locationFile;

	//multi sample in one vcf file
	//bool multi;

	//call variant only
	bool varOnly;

	//method
	int method;

	//mutation rate
	double mRate;

	//Pr(genotype) for new variant
	std::vector<double> genoProbN;

	//Pr(genotype) for known variant
	std::vector<double> genoProbK;

	//Pr(genotype) for new variant on chromosome X for male
	std::vector<double> genoProbXN;

	//Pr(genotype) for known variant on chromosome X for male
	std::vector<double> genoProbXK;

	//whether the position in vcf file has been ordered
	bool posOrder;

	//output every line in the vcf file, it will be blocked when varOnly is true
	bool allLine;

	//only output the line with different called genotype with the original vcf file
	bool diffOnly;

	//number of burn in in MCMC
	int numBurnIn;

	//number of repeat in MCMC
	int numRep;

	//whether we need to calculate
	//if likelihood larger than LRCrteria, we don't calculate for that locus
	double LCriteria;

public:
	inputVCF();

	inputVCF(const std::vector<std::string> & FileNameVCF, std::string FileNamePed, std::string OutputName, std::string LocationFile,bool Multi, bool VarOnly, int Method, double MRate, bool PosOrder);

	//////////////////////////////////
	///       get functions        ///
	//////////////////////////////////

	//fileNameVCF
	std::vector<std::string> get_fileNameVCF() const {return fileNameVCF;}

	//fileNamePed
	std::string get_fileNamePed() const {return fileNamePed;}

	//outputName
	std::string get_outputName() const {return outputName;}

	//locationFile
	std::string get_locationFile() const {return locationFile;}

	//multi
	//bool get_multi() const {return multi;}

	//varOnly
	bool get_varOnly() const {return varOnly;}

	//method
	int get_method() const { return method;}

	//mRate
	double get_mRate() const{ return mRate;}

	//genoProbN
	std::vector<double> get_genoProbN() const{ return genoProbN;}

	//genoProbK
	std::vector<double> get_genoProbK() const{ return genoProbK;}

	//genoProbXN
	std::vector<double> get_genoProbXN() const {return genoProbXN;}

	//genoProbXK
	std::vector<double> get_genoProbXK() const {return genoProbXK;}

	//posOrder
	bool get_posOrder() const {return posOrder;}

	//allLine
	bool get_allLine() const {return allLine;}

	//diffOnly
	bool get_diffOnly() const { return diffOnly;}

	//burn in
	int get_numBurinIn() const {return numBurnIn;}

	//repeat times
	int get_numRep() const {return numRep;}

	// LRCriteria
	double get_LCriteria() const{ return LCriteria;}


	//////////////////////////////////
	///       set functions        ///
	//////////////////////////////////

	//fileNameVCF
	bool set_fileNameVCF(const std::vector<std::string> & FileNameVCF);

	//fileNamePed
	bool set_fileNamePed(std::string FileNamePed);

	//outputName
	bool set_outputName(std::string OutputName);

	//locationFile
	bool set_locationFile(std::string LocationFile);

	//multi
	//bool set_multi(bool Multi);

	//varOnly
	bool set_varOnly(bool VarOnly);

	//method
	bool set_method(int Method);

	//mRate
	bool set_mRate(double MRate);

	//genoProbN
	bool set_genoProbN(const std::vector<double> & GenoProbN);

	//genoProbK
	bool set_genoProbK(const std::vector<double> & GenoProbK);

	//genoProbXN
	bool set_genoProbXN(const std::vector<double> & GenoProbXN);

	//genoProbXK
	bool set_genoProbXK(const std::vector<double> & GenoProbXK);

	//posOrder
	bool set_posOrder(bool PosOrder);

	//allLine
	bool set_allLine(bool AllLine);

	//diffOnly
	bool set_diffOnly(bool DiffOnly);

	//number of burn in
	bool set_numBurnIn(int NumBurnIn);

	//number of repeat
	bool set_numRep(int NumRep);

	// LRCriteria
	bool set_LCriteria(double LC);
};

class inputLK
{
private:
	//likelihood file name
	std::string fileNameLK;

	//ped file name
	std::string fileNamePed;

	//output file name
	std::string outputName;

	//method
	int method;

	//mutation rate
	double mRate;

	//likelihood type
	//1: normal
	//2: log10
	//3: ln
	//4: Phred-scaled
	int lkType;

	//Pr(genotype) for new variant
	std::vector<double> genoProbN;

	//Pr(genotype) for known variant
	std::vector<double> genoProbK;

	//Pr(genotype) for new variant on chromosome X for male
	std::vector<double> genoProbXN;

	//Pr(genotype) for known variant on chromosome X for male
	std::vector<double> genoProbXK;

	//number of burn in in MCMC
	int numBurnIn;

	//number of repeat in MCMC
	int numRep;

	//whether we need to calculate
	//if likelihood larger than LRCrteria, we don't calculate for that locus
	double LCriteria;

public:
	inputLK();

	inputLK(const std::string & FileNameLK, const std::string & FileNamePed, const std::string & OutputName, int Method, double MRate);

	//////////////////////////////////
	///       get functions        ///
	//////////////////////////////////

	//fileNameVCF
	std::string get_fileNameLK() const {return fileNameLK;}

	//fileNamePed
	std::string get_fileNamePed() const {return fileNamePed;}

	//outputName
	std::string get_outputName() const {return outputName;}

	//mRate
	double get_mRate() const{ return mRate;}

	//method
	int get_method() const{ return method;}

	//likelihood type
	int get_lkType() const{ return lkType;}

	//genoProbN
	std::vector<double> get_genoProbN() const{ return genoProbN;}

	//genoProbK
	std::vector<double> get_genoProbK() const{ return genoProbK;}

	//genoProbXN
	std::vector<double> get_genoProbXN() const {return genoProbXN;}

	//genoProbXK
	std::vector<double> get_genoProbXK() const {return genoProbXK;}

	//burn in
	int get_numBurinIn() const {return numBurnIn;}

	//repeat times
	int get_numRep() const {return numRep;}

	// LRCriteria
	double get_LCriteria() const{ return LCriteria;}

	//////////////////////////////////
	///       set functions        ///
	//////////////////////////////////

	//fileNameVCF
	bool set_fileNameLK(std::string FileNameLK);

	//fileNamePed
	bool set_fileNamePed(std::string FileNamePed);

	//outputName
	bool set_outputName(std::string OutputName);

	//method
	bool set_method(int Method);

	//mRate
	bool set_mRate(double MRate);

	//lkType
	bool set_lkType(int LKType);

	//genoProbN
	bool set_genoProbN(const std::vector<double> & GenoProbN);

	//genoProbK
	bool set_genoProbK(const std::vector<double> & GenoProbK);

	//genoProbXN
	bool set_genoProbXN(const std::vector<double> & GenoProbXN);

	//genoProbXK
	bool set_genoProbXK(const std::vector<double> & GenoProbXK);

	//number of burn in
	bool set_numBurnIn(int NumBurnIn);

	//number of repeat
	bool set_numRep(int NumRep);

	// LRCriteria
	bool set_LCriteria(double LC);
};

//check the options when input vcf file
//input:
//argc:
//argv:
//output:
//vcfInfo: parameter information
//return
//0: good
//1: warning
//-1: stop
int checkInputVCF(int argc, char* argv[], inputVCF & vcfInfo);

//check the options when input likelihood file
//input:
//argc:
//argv:
//output:
//lkInfo: parameter information
//return
//0: good
//1: warning
//-1: stop
int checkInputLK(int argc, char* argv[], inputLK & lkInfo);

#endif /* CHECKINPUT_H_ */
