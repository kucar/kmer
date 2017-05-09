#ifndef KMER_HASH_H
#define KMER_HASH_H
#include "bloom_filter.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>

#define UNIQUE 1
#define SEQ_READ_BUFFER 10000
#define NUMLINES_ENTRY 4
#define VALID_ENTRY    2
#define SAMPLE_FILE_SIZE 47363704
#define SAMPLE_FILE_LINENUMBER 799440
#define HASHSIZE  (1ULL<<25)
//hash size threshold per shrink
#define THRES_PER_SHR ((HASHSIZE>>2)+(HASHSIZE>>4)) //(10*(1ULL<<20))
#define CALLFILTERFUNC(x,y) (this->*m_filter_ptr_arr[m_filtertyp])(x,y)


typedef unsigned long long int u_32bits;
typedef unsigned char uch;
class KMER_COUNTER;
typedef  unsigned int kint;
typedef std::pair<std::string,kint>  kmer_entry;
typedef std::pair<u_32bits, kint> kmer_entry_zip;
typedef std::unordered_map<u_32bits,kint> kmerhash;

typedef void (KMER_COUNTER::*filter_ptr)(std::string&,unsigned int);

enum filter_t {
	shrink,
	bloom,
	nofilter
};

//Base class holds basic attributes for kmer processing.
//Example :
//std::shared_ptr<KMER_COUNTER> kmerobj
//									(new KMER_COUNTER(filename,
//													topcount,
//													kmersize,
//													filtertyp,
//													stats));
class KMER_COUNTER 
{
public:

	explicit KMER_COUNTER(std::string & filename,kint topcount,kint kmersize, bool stats,filter_t filter,int ram,unsigned long long num_insertions,unsigned int linelength):
		m_fastq_file(filename),
		m_topcount(topcount),
		m_kmer_size(kmersize),
		m_stats(stats),
		m_filtertyp(filter),
		m_shrink_cnt(0),
		m_ram(ram),
		m_estimated_insertions(num_insertions),
		m_linelength(linelength)
		{
			m_filter_ptr_arr[0]=&KMER_COUNTER::Shrink_Insert;
			m_filter_ptr_arr[1]=&KMER_COUNTER::Bloom_Insert;
			m_filter_ptr_arr[2]=&KMER_COUNTER::Insert;
		};
	virtual ~KMER_COUNTER();
	virtual void Begin(void);
	virtual void Init(void);
	virtual void FindTopN(void);
	virtual void PrintStats(clock_t , clock_t );
	virtual void ClearSequenceHash(void);
	virtual void ClearTopVector(void);
	virtual void ShrinkHash(void);
	virtual void BloomInit(void);
	virtual void Bloom_Insert(std::string& ,unsigned int);
	virtual void Shrink_Insert(std::string& ,unsigned int);
	virtual void Insert(std::string& ,unsigned int);
protected:
	std::string m_fastq_file;
	kint m_topcount;
	kint m_kmer_size;
	bool m_stats;
	std::vector<kmer_entry_zip> m_topNvector_zip;
	kmerhash m_sequencehash_zip;
	filter_t m_filtertyp;
	bloom_filter m_bloom_filter;
	kint m_shrink_cnt;
	int m_ram;
	filter_ptr m_filter_ptr_arr[3];
	unsigned long long m_estimated_insertions;
	unsigned int m_linelength;
};

#endif //KMER_HASH_H
