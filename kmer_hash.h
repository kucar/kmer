#ifndef KMER_HASH_H
#define KMER_HASH_H

#include "bloom_filter.hpp"
#include <unordered_map>
#include <map>
#include <vector>
#include <string>

#define UNIQUE 1
typedef  unsigned int kint;
typedef std::pair<std::string,kint>  kmer_entry;
typedef std::pair<unsigned long long, kint> kmer_entry_zip;


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
	explicit KMER_BASE(std::string & filename,kint topcount,kint kmersize, bool stats,bool bloom):
		m_fastq_file(filename),
		m_topcount(topcount),
		m_kmer_size(kmersize),
		m_stats(stats),
		m_bloom(bloom),
		m_shrink_cnt(0){};
	virtual ~KMER_BASE();
	virtual void Begin(void);
	virtual void Init(void);
	virtual void FindTopN(void);
	virtual void PrintStats(clock_t , clock_t );
	virtual void ClearSequenceHash(void);
	virtual void ClearTopVector(void);
	virtual void ShrinkHash(void);
	virtual void BloomInit(void);
	virtual void Bloomify(std::string& ,unsigned int);
protected:
	std::string m_fastq_file;
	kint m_topcount;
	kint m_kmer_size;
	bool m_stats=false;
	std::vector<kmer_entry_zip> m_topNvector_zip;
	std::unordered_map<unsigned long long,kint> m_sequencehash_zip;
	bool m_bloom;
	bloom_filter m_bloom_filter;
	kint m_shrink_cnt;
};


#endif //KMER_HASH_H
