/*
 * nucleotide_defs.h
 *
 *  Created on: Apr 11, 2017
 *      Author: kucar
 */

#ifndef NUCLEOTIDE_DEFS_H_
#define NUCLEOTIDE_DEFS_H_
enum nucleotideCode_t {
	A = 'A',
	C = 'C',
	G = 'G',
	T = 'T'
};


int indexMapper[T + 1] = {0};

char indexRevMapper[4]={0};



#endif /* NUCLEOTIDE_DEFS_H_ */
