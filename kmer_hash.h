#ifndef KMER_HASH_H
#define KMER_HASH_H
#include "bloom_filter.hpp"
#include <unordered_map>
#include <vector>
#include <string>
#include <fstream>
#include <map>

#define UNIQUE         1
#define NUMLINES_ENTRY 4
#define VALID_ENTRY    2

//predefined base hash size
#define HASHSIZE  (1ULL<<25)
//predefined base hash size threshold per shrink
#define THRES_PER_SHR ((HASHSIZE>>2)+(HASHSIZE>>4)) //(10*1024*1024)


typedef unsigned long long int u_64bits;
typedef unsigned char uchar;
class KMER_COUNTER;
typedef  unsigned int kmer_int;
typedef std::pair<std::string,kmer_int>  kmer_entry;
typedef std::pair<u_64bits, kmer_int> kmer_entry_zip;
typedef std::unordered_map<u_64bits,kmer_int> kmerhash;
typedef void (KMER_COUNTER::*filter_ptr)(std::string&,kmer_int);

enum filter_t {
	shrink,
	bloom,
	nofilter
};

enum mem_usage_t{
	minn=0x10,
	low,
	medium,
	high,
	maxx
};
typedef std::map<mem_usage_t,float> loadfactor_thr_t;
typedef std::map<std::string,mem_usage_t> memusage_map_t;



//Base class holds basic attributes and methods for kmer processing.
//Example :
//std::shared_ptr<KMER_COUNTER> kmerobj
//									                 (new KMER_COUNTER(filename,
//													 topcount,
//													 kmersize,
//													 stats,
//													 filtertype,
//													 mem_usage,
//													 expected_number_insertions,
//                                                   nucleotide_line_length_infile
//												     ));
class KMER_COUNTER 
{
public:

	explicit KMER_COUNTER(std::string & filename,kmer_int topcount,
						  kmer_int kmersize, bool stats,
						  filter_t filter,mem_usage_t memusg,
						  u_64bits num_insertions,kmer_int linelength):
						  m_fastq_file(filename),
						  m_topcount(topcount),
						  m_kmer_size(kmersize),
						  m_stats(stats),
						  m_filtertyp(filter),
						  m_shrink_cnt(0),
						  m_memusg(memusg),
						  m_estimated_insertions(num_insertions),
						  m_linelength(linelength),
						  m_empty(false)
	{
		m_filter_ptr_arr[0]=&KMER_COUNTER::Shrink_Insert;
		m_filter_ptr_arr[1]=&KMER_COUNTER::Bloom_Insert;
		m_filter_ptr_arr[2]=&KMER_COUNTER::Insert;
	};
	explicit KMER_COUNTER():
						m_fastq_file(""),
						m_topcount(0),
						m_kmer_size(0),
						m_stats(0),
						m_filtertyp(nofilter),
						m_shrink_cnt(0),
						m_memusg(minn),
						m_filter_ptr_arr(),
						m_estimated_insertions(0),
						m_linelength(0),
						m_empty(true){};
	//move constructer
	KMER_COUNTER(KMER_COUNTER&& other):
				m_fastq_file(other.m_fastq_file),
				m_topcount(other.m_topcount),
				m_kmer_size(other.m_kmer_size),
				m_stats(other.m_stats),
				m_filtertyp(other.m_filtertyp),
				m_shrink_cnt(other.m_shrink_cnt),
				m_memusg(other.m_memusg),
				m_estimated_insertions(other.m_estimated_insertions),
				m_linelength(other.m_linelength),
				m_empty(false)
	{
		m_filter_ptr_arr[0]=&KMER_COUNTER::Shrink_Insert;
		m_filter_ptr_arr[1]=&KMER_COUNTER::Bloom_Insert;
		m_filter_ptr_arr[2]=&KMER_COUNTER::Insert;
		other.ClearSequenceHash();
		other.ClearTopVector();
		other.m_empty=true;
	};
	KMER_COUNTER& operator=(KMER_COUNTER&& other);
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
	kmer_int m_topcount;
	kmer_int m_kmer_size;
	bool m_stats;
	std::vector<kmer_entry_zip> m_topNvector_zip;
	kmerhash m_sequencehash_zip;
	filter_t m_filtertyp;
	bloom_filter m_bloom_filter;
	kmer_int m_shrink_cnt;
	mem_usage_t m_memusg;
	filter_ptr m_filter_ptr_arr[3];
	unsigned long long m_estimated_insertions;
	//assume each and every processed line has same number of chars.
	kmer_int m_linelength;
	bool m_empty;
};

#endif //KMER_HASH_H
