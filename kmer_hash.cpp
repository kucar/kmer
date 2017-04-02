#include "kmer_hash.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <fstream>
#include <limits>
#include <ctime>

#define NUMLINES_ENTRY 4
#define VALID_ENTRY    2

KMER_BASE::~KMER_BASE(){};

//starts up the map and result vector(histogram winners)
void KMER_BASE::Init(void)
{
	m_sequencehash.reserve(HASHSIZE);		  //prevent rehashing
	m_topNvector.resize(m_topcount);
	std::string line;
	kint linenum =1 ;
	std::ifstream myfile (m_fastq_file);
	if (myfile.is_open())
	{
		while ( getline (myfile,line) )
		{
			if(linenum % NUMLINES_ENTRY == VALID_ENTRY) //process sequence
			{
				for ( unsigned int index=0;
					  index<=(((m_line_size_static) ? SIZE_LINE : line.length())-m_kmer_size);
					  ++index)
				{
					++m_sequencehash[line.substr(index,m_kmer_size)];
				}
				linenum=(linenum+1);  //next line
			}

			else  //skip this line
			{
				linenum++;
				continue;
			}
		}
		myfile.close();
	}

	else 
	{
		std::cout<<"ERROR: FILE CANNOT BE OPENED::"<<m_fastq_file<<std::endl;
		exit(EXIT_FAILURE);
	}

};

//sorts the map according to the count and prints the top N winners
void KMER_BASE::FindTopN()
{
	std::partial_sort_copy(m_sequencehash.begin(), m_sequencehash.end(),m_topNvector.begin(),m_topNvector.end(),
			[](std::pair<const std::string, kint> const& l,
					std::pair<const std::string, kint> const& r)
					{	return l.second > r.second;	});
	for (auto const& p: m_topNvector)
	{
		std::cout << "{ " << p.first << ", " << p.second << "}\n";
	}
};

void KMER_BASE::Begin()
{
	this->Init();
	this->FindTopN();
}
void KMER_BASE::PrintStats(clock_t _begin, clock_t _end)
{
	std::cout << "hash size = " << m_sequencehash.size() << std::endl;
	std::cout << "bucket count = " << m_sequencehash.bucket_count() << std::endl;
	std::cout << "load_factor = " << m_sequencehash.load_factor() << std::endl;
	std::cout << "max_load_factor = " << m_sequencehash.max_load_factor() << std::endl;
	double elapsed_secs = double(_end - _begin) / CLOCKS_PER_SEC;
	std::cout<<"elapsed sec:"<<elapsed_secs<<std::endl;
	std::cout<<"top vector size: "<<m_topNvector.size()<<std::endl;
}
