#include "kmer_hash.h"
#include "nucleotide_defs.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <fstream>
#include <limits>
#include <ctime>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define NUMLINES_ENTRY 4
#define VALID_ENTRY    2
#define HASHSIZE  (1ULL<<25)
#define HASHSZ_THRESHOLD_PER_SHRINK 250000


inline int getMappedCode(nucleotideCode_t n) {
	return indexMapper[n];
}
inline char getRevMappedCode(int n){
	return indexRevMapper[n];
}

// largest number we're going to hash into. (8 bytes/64 bits/32 nt)
typedef unsigned long long int HashIntoType;
const unsigned char KSIZE_MAX = sizeof(HashIntoType)*4;

// largest size 'k' value for k-mer calculations.  (1 byte/255)
typedef unsigned char WordLength;



HashIntoType _hash(const char * kmer, const WordLength k,
                   HashIntoType& _h)
{
    // sizeof(HashIntoType) * 8 bits / 2 bits/base
    if (k > sizeof(HashIntoType)*4) {
        std::cout<<"Supplied kmer string doesn't match the underlying k-size."<<std::endl;
    }

    if (strlen(kmer) < k) {
        std::cout<<"kmer is too short"<<std::endl;

    }

    HashIntoType h = 0;

//    h |= twobit_repr(kmer[0]);
      h |=getMappedCode((nucleotideCode_t)kmer[0]);

    for (WordLength i = 1; i < k; i++) {
        h = h << 2;
        h |=getMappedCode((nucleotideCode_t)kmer[i]);
    }
    _h = h;
    return h;
}

// _hash: return the maximum of the forward and reverse hash.

HashIntoType _hash(const char * kmer, const WordLength k)
{
    HashIntoType h = 0;
    return _hash(kmer, k, h);
}

// _hash_forward: return the hash from the forward direction only.

HashIntoType _hash_forward(const char * kmer, WordLength k)
{
    HashIntoType h = 0;
 
    _hash(kmer, k, h);
    return h;			// return forward only
}

HashIntoType _hash(const std::string kmer, const WordLength k)
{
    return _hash(kmer.c_str(), k);
}

HashIntoType _hash(const std::string kmer, const WordLength k,
                   HashIntoType& h)
{
    return _hash(kmer.c_str(), k, h);
}

std::string _revhash(HashIntoType hash, WordLength k)
{
    std::string s = "";

    unsigned int val = hash & 3;
      s+=getRevMappedCode(val);

    for (WordLength i = 1; i < k; i++) {
        hash = hash >> 2;
        val = hash & 3;
        s+=getRevMappedCode(val);
    }

    reverse(s.begin(), s.end());

    return s;
}

KMER_BASE::~KMER_BASE(){};

//starts up the map and result vector(histogram winners)
void KMER_BASE::Init(void)
{
	unsigned int shrink_threshold=HASHSZ_THRESHOLD_PER_SHRINK;
	m_sequencehash_zip.rehash(HASHSIZE);
	std::string line;
	std::string subs;
	std::ifstream myfile (m_fastq_file);
	kint linenum =0 ;
	if(m_bloom)
		BloomInit();
	if (myfile.is_open())
	{
		while(getline(myfile,line))
		{
			if(++linenum % NUMLINES_ENTRY == VALID_ENTRY) //process sequence
			{
				unsigned int limit=line.length()-m_kmer_size;
				if (m_bloom)
				{
					Bloomify(line,limit);
				}
				else
				{
					for ( unsigned int index=0; index<=limit;  ++index)
					{
						subs=line.substr(index,m_kmer_size);
						++m_sequencehash_zip[_hash(subs, m_kmer_size)];
					}
				}
			}
			if(!m_bloom)
				if(m_sequencehash_zip.size()>shrink_threshold)
				{
					ShrinkHash();
					++m_shrink_cnt;
					shrink_threshold+=HASHSZ_THRESHOLD_PER_SHRINK;
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

void KMER_BASE::ClearSequenceHash()
{
	m_sequencehash_zip.clear();
	m_sequencehash_zip.rehash(HASHSIZE);
}
void KMER_BASE::ClearTopVector()
{
	m_topNvector_zip.clear();
	m_topNvector_zip.resize(m_topcount);
}

void KMER_BASE::FindTopN()
{
	m_topNvector_zip.resize(m_topcount);
	std::partial_sort_copy(m_sequencehash_zip.begin(), m_sequencehash_zip.end(),m_topNvector_zip.begin(),m_topNvector_zip.end(),
			[](kmer_entry_zip const& l, kmer_entry_zip const& r)
					{	return l.second > r.second;	});
	if(!m_bloom)
		for (auto it:m_topNvector_zip)
		{
			std::cout<<_revhash(it.first,m_kmer_size)<<"---"<<
			it.second<<" (fn rate="<<((long double)m_shrink_cnt)/((long double)it.second)<<
			")"<<std::endl;
		}
	else
		for (auto it:m_topNvector_zip)
			std::cout << _revhash(it.first,m_kmer_size)<<"---"<<it.second<<std::endl;
};
void KMER_BASE::Begin()
{
	this->Init();
	this->FindTopN();
}
void KMER_BASE::PrintStats(clock_t _begin, clock_t _end)
{
	std::cout << "bucket count = " << m_sequencehash_zip.bucket_count() << std::endl;
	std::cout << "load_factor = " << m_sequencehash_zip.load_factor() << std::endl;
	double elapsed_secs = double(_end - _begin) / CLOCKS_PER_SEC;
	std::cout<<"top vector size: "<<m_topNvector_zip.size()<<std::endl;
	std::cout<<"shrink :"<<!m_bloom<<std::endl;
	std::cout<<"bloom : "<<m_bloom<<std::endl;
	if (m_bloom)
	{
		std::cout<<"effective false positive (bloom) "<<m_bloom_filter.effective_fpp()<<std::endl;
	}
	else
	{
		std::cout<<"shrink count (false negatives) :"<<m_shrink_cnt<<std::endl;
		std::cout<<"hash size per shrink threshold:"<<HASHSZ_THRESHOLD_PER_SHRINK<<std::endl;

	}
	std::cout<<"elapsed sec before destructor operation:"<<elapsed_secs<<std::endl;
}

void KMER_BASE::ShrinkHash(void)
{
	for (auto it=m_sequencehash_zip.begin();it!=m_sequencehash_zip.end();)
	{
		if(it->second==UNIQUE)
		{
			it = m_sequencehash_zip.erase(it);
		}
		else 
			it++;
	}
}
void KMER_BASE::BloomInit(void)
{
	bloom_parameters parameters;
	parameters.projected_element_count = 13000000;
	parameters.false_positive_probability = 0.05;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.compute_optimal_parameters();
	bloom_filter filter(parameters);
	m_bloom_filter=filter;
	return;
}
void KMER_BASE::Bloomify(std::string& _line,unsigned int _limit)
{
	for ( unsigned int index=0; index<=_limit;  ++index)
	{
		auto subs=_line.substr(index,m_kmer_size);
		if (m_bloom_filter.contains(subs))
		{
			++m_sequencehash_zip[_hash(subs, m_kmer_size)];
		}
		else
			m_bloom_filter.insert(subs);
	}
}
__attribute__((constructor))
static void _initializeIndexMapper() {
	indexMapper[A] = 0;
	indexMapper[C] = 1;
	indexMapper[G] = 2;
	indexMapper[T] = 3;
	indexRevMapper[0]='A';
	indexRevMapper[1]='C';
	indexRevMapper[2]='G';
	indexRevMapper[3]='T';
}
