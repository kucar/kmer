/*
 * file_op.h
 *
 *  Created on: May 9, 2017
 *      Author: kucar
 */

#ifndef FILE_OP_H_
#define FILE_OP_H_

#include <fstream>
#include <iostream>
#include <string>

#define SAMPLE_FILE_SIZE 47363704  //sample file 47 mb
#define SAMPLE_FILE_LINENUMBER 799440
//get file size and nucleotide line length
inline unsigned long long getFileSize(std::string & filename,unsigned int& length)
{
	std::streampos begin,end;
	std::string line;
	std::ifstream myfile (filename, std::ios::in);
	//scan through second line
	{
		getline(myfile,line);
		getline(myfile,line);
	}
	//get the number of chars in nucleotide line
	length=line.length();
	//rewind to begining
	myfile.seekg (0, std::ios::beg);
#ifdef DEBUG_M
	std::cout<<"line length to be processed:"<<length<<std::endl;
#endif
	begin = myfile.tellg();
	myfile.seekg (0, std::ios::end);
	end = myfile.tellg();
	myfile.close();
	return (end-begin);
}

inline unsigned long long estimate_line_num(unsigned long long fsize)
{
	int multiple = fsize / SAMPLE_FILE_SIZE;
	int remaing_size=(multiple)?  fsize % (multiple*SAMPLE_FILE_SIZE): fsize;
#ifdef DEBUG_M
	std::cout<<"line num estimate dbg: multiple:"<<multiple<<std::endl;
	std::cout<<"line num estimate dbg:remaing_size:"<<remaing_size<<std::endl;
#endif
	return (multiple*SAMPLE_FILE_LINENUMBER)+(((double)remaing_size/(double)SAMPLE_FILE_SIZE)*SAMPLE_FILE_LINENUMBER);
}
inline unsigned long long estimate_insertions(unsigned long long line_num,int line_length,int kmersize)
{
	auto eff_line_num=line_num/NUMLINES_ENTRY;
	return eff_line_num*(line_length-kmersize);
}

#endif /* FILE_OP_H_ */
