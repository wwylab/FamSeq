/*
 * file.h
 *
 * read vcf and other files
 *
 *  Created on: Mar 16, 2012
 *  Author: Gang Peng <gpeng1@mdanderson.org>
 *
 *  FamSeq is free software. You can redistribute and/or modify it under GNU General Public License
 *  of version 3(GPLv3).
 */

#ifndef FILE_H_
#define FILE_H_

#include "family.h"
#include "checkInput.h"

//read pedigree information file
//input: fileName
//output: mem
//return false when
//File cannot be opened
//ID is not set correctly
bool readPed(std::string fileName, std::vector<individual> & mem);

//read the vcf file with only one sample
//the output file are organized according to chromosome
//input:
//fileName
//output:
//header: the header of the vcf file
//pos: position
//id: ID
//ref: reference
//alt: alternative
//dp: depther
//others: other information
//genotype: genotype called by other method
//pl1,pl2,pl3: scaled likelihood
//return false when
//Cannot open vcf file
bool readVCF(std::string fileName, std::string & header, std::vector<std::vector<int> > & pos, std::vector<std::vector<std::string> > & id, std::vector<std::vector<char> > & ref, std::vector<std::vector<char> > & alt, std::vector<std::vector<int> > & dp, std::vector<std::vector<std::vector<std::string> > > & others, std::vector<std::vector<std::string> > & genotype, std::vector<std::vector<int> > & pl1, std::vector<std::vector<int> > & pl2, std::vector<std::vector<int> > & pl3);

//when multi samples are in one vcf file, call the genotype
//input:
//vcfFileName
//outFileName: the file in which the result will be written
//fam: store the family information
//sampleInd: sample indicator. match the sample order in fam and sample order in vcf file
//VarOnly: If true, only output the loci with variant
//return false when
//Cannot open vcf file
//bool callGenoMVCF(std::string vcfFileName, std::string outFileName, family fam, bool VarOnly, int method);
bool callGenoMVCF(const inputVCF & vcfInfo, family fam);

//when input likelihood file, call the genotype
//input:
//inputLK: input information
//fam: family information
//return:
//true: we could call the genotype
//false: something wrong
bool callGenoLK(const inputLK & lkInfo, family fam);

//set family according the information in vcfInfo
//input:
//inputVCF
//output:
//fam
//return
bool setFam(const inputVCF & vcfInfo, family &fam);

//set family according the information in lkInfo
//input:
//inputLK
//output:
//fam
//return
bool setFam(const inputLK & lkInfo, family & fam);

#endif /* FILE_H_ */
