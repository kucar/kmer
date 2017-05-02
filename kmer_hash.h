#ifndef KMER_HASH_H
#define KMER_HASH_H
#include "bloom_filter.hpp"
#include <unordered_map>
#include <map>
#include <vector>
#include <string>

#define UNIQUE 1
#define SEQ_READ_BUFFER 10000
class KMER_BASE;
typedef  unsigned int kint;
typedef std::pair<std::string,kint>  kmer_entry;
typedef std::pair<unsigned long long, kint> kmer_entry_zip;
typedef std::vector<std::string> readbuf;
typedef std::unordered_map<unsigned long long,kint> kmerhash;

typedef void (KMER_BASE::*filter_ptr)(std::string&,unsigned int);

enum filter_t {
	shrink,
	bloom
};

//Base class holds basic attributes for kmer processing.
//This class can be derived for other kmer counting implementations in future
//Example :
//std::shared_ptr<KMER_BASE> kmerobj
//									(new KMER_BASE(filename,
//													topcount,
//													linesizestatic,
//													kmersize,
//													stats));
class KMER_BASE 
{
public:

	explicit KMER_BASE(std::string & filename,kint topcount,kint kmersize, bool stats,filter_t filter,bool parallel,int ram):
		m_fastq_file(filename),
		m_topcount(topcount),
		m_kmer_size(kmersize),
		m_stats(stats),
		m_filtertyp(filter),
		m_shrink_cnt(0),
		m_parallel(parallel),
		m_ram(ram)
	{
		m_filter_ptr_arr[0]=&KMER_BASE::Shrink_Insert;
		m_filter_ptr_arr[1]=&KMER_BASE::Bloom_Insert;
	};
	virtual ~KMER_BASE();
	virtual void Begin(void);
	virtual void Init(void);
	virtual void Init(bool parallel);
	virtual void FindTopN(void);
	virtual void PrintStats(clock_t , clock_t );
	virtual void ClearSequenceHash(void);
	virtual void ClearTopVector(void);
	virtual void ShrinkHash(void);
	virtual void BloomInit(void);
	virtual void Bloom_Insert(std::string& ,unsigned int);
	virtual void Shrink_Insert(std::string& ,unsigned int);
	virtual void FileParser(std::string& );
	virtual void ProcessBuffer(void);
protected:
	std::string m_fastq_file;
	kint m_topcount;
	kint m_kmer_size;
	bool m_stats=false;
	std::vector<kmer_entry_zip> m_topNvector_zip;
	kmerhash m_sequencehash_zip;
	filter_t m_filtertyp;
	bloom_filter m_bloom_filter;
	kint m_shrink_cnt;
	readbuf m_buffer;
	bool m_parallel;
	int m_ram;
	filter_ptr m_filter_ptr_arr[2];

};


#endif //KMER_HASH_H
