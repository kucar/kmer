#ifndef KMER_HASH_H
#define KMER_HASH_H

#include <unordered_map>
#include <vector>
#include <string>


#define HASHSIZE (1<<25)   //2 pow 25

typedef  unsigned int kint;
typedef std::pair<std::string,kint>  kmer_entry;

#define SIZE_LINE 90


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
	explicit KMER_BASE(std::string & filename,kint topcount,bool linesizestatic,kint kmersize, bool stats):
		m_fastq_file(filename),
		m_topcount(topcount),
		m_line_size_static(linesizestatic),
		m_kmer_size(kmersize),m_stats(stats){};
	virtual ~KMER_BASE();
	virtual void Begin(void);
	virtual void Init(void);
	virtual void FindTopN(void);
	virtual void PrintStats(clock_t , clock_t );
protected:
	std::string m_fastq_file;
	kint m_topcount;
	bool m_line_size_static;
	kint m_kmer_size;
	bool m_stats=false;
	std::vector<kmer_entry> m_topNvector;
	std::unordered_map <std::string,kint> m_sequencehash;

};


#endif //KMER_HASH_H
