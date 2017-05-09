#include "kmer_hash.h"
#include "nucleotide_defs.h"
#include <iostream>
#include <algorithm>
#include <utility>
#include <limits>
#include <ctime>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>






inline uint32_t tune_shrink_threshold(int xgbram)  {
	return (xgbram>=16)? THRES_PER_SHR<<3 : (xgbram>=8) ? THRES_PER_SHR<<2 : (xgbram>=4) ? THRES_PER_SHR<<1 : THRES_PER_SHR ;
}

inline uint32_t tune_hashsize(int xgbram)
{
	return  (xgbram>=16)? HASHSIZE<<3 : (xgbram>=8) ? HASHSIZE<<2 : (xgbram>=4) ? HASHSIZE<<1 : HASHSIZE ;
}
inline int getMappedCode(nucleotideCode_t n) {
	return indexMapper[n];
}
inline char getRevMappedCode(int n){
	return indexRevMapper[n];
}



inline u_32bits char2bit(const char * kmer, const uch k,u_32bits& charsinbits )
{

    u_32bits _charsinbits = 0;

      _charsinbits |=getMappedCode((nucleotideCode_t)kmer[0]);

    for (uch i = 1; i < k; i++) {
        _charsinbits = _charsinbits << 2;
        _charsinbits |=getMappedCode((nucleotideCode_t)kmer[i]);
    }
    charsinbits = _charsinbits;
    return charsinbits;
}

inline u_32bits char2bit(const char * kmer, const uch k)
{
    u_32bits h = 0;
    return char2bit(kmer, k, h);
}

inline u_32bits char2bit(const std::string kmer, const uch k)
{
    return char2bit(kmer.c_str(), k);
}

std::string bit2str(u_32bits charsinbits, uch k)
{
    std::string s = "";

    unsigned int val = charsinbits & 3;
      s+=getRevMappedCode(val);

    for (uch i = 1; i < k; i++) {
        charsinbits = charsinbits >> 2;
        val = charsinbits & 3;
        s+=getRevMappedCode(val);
    }

    reverse(s.begin(), s.end());

    return s;
}


KMER_COUNTER::~KMER_COUNTER(){};

//starts up the map and result vector(histogram winners)
void KMER_COUNTER::Init(void)
{
	uint32_t _size=tune_hashsize(m_ram);
	m_sequencehash_zip.rehash(_size);
	kint limit;
	std::string line;
	std::ifstream myfile (m_fastq_file);
	kint linenum =0 ;
	auto linelength= m_linelength - m_kmer_size;
	if(m_filtertyp==bloom)
		BloomInit();
	if (myfile.is_open())
	{
		while(getline(myfile,line))
			if (++linenum % NUMLINES_ENTRY == VALID_ENTRY) 
			{   //process sequence 
				CALLFILTERFUNC(line,linelength);
			}
		myfile.close();
	}
	else
	{
		std::cout<<"ERROR: FILE CANNOT BE OPENED::"<<m_fastq_file<<std::endl;
		exit(EXIT_FAILURE);
	}
};

void KMER_COUNTER::ClearSequenceHash()
{
	m_sequencehash_zip.clear();
	m_sequencehash_zip.rehash(HASHSIZE);
}
void KMER_COUNTER::ClearTopVector()
{
	m_topNvector_zip.clear();
	m_topNvector_zip.resize(m_topcount);
}

void KMER_COUNTER::FindTopN()
{
	volatile clock_t begin = clock();
	m_topNvector_zip.resize(m_topcount);
	std::partial_sort_copy(m_sequencehash_zip.begin(), m_sequencehash_zip.end(),m_topNvector_zip.begin(),m_topNvector_zip.end(),
			[](kmer_entry_zip const& l, kmer_entry_zip const& r)
					{	return l.second > r.second;	});
	for (auto it:m_topNvector_zip)
			std::cout<<bit2str(it.first,m_kmer_size)<<"---"<<it.second<<std::endl;
	volatile clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
#ifdef DEBUG_M
	std::cout<<"elapsed sec for findTopN:"<<elapsed_secs<<std::endl;
#endif
};
void KMER_COUNTER::Begin()
{
	this->Init();
	this->FindTopN();
}
void KMER_COUNTER::PrintStats(clock_t _begin, clock_t _end)
{
	std::cout << "bucket count = " << m_sequencehash_zip.bucket_count() << std::endl;
	std::cout << "hash size = " << m_sequencehash_zip.size() << std::endl;
	std::cout << "load_factor = " << m_sequencehash_zip.load_factor() << std::endl;
	double elapsed_secs = double(_end - _begin) / CLOCKS_PER_SEC;
	std::cout<<"top vector size: "<<m_topNvector_zip.size()<<std::endl;
	std::string filter_typ=(m_filtertyp==bloom)? "bloom filter" : (m_filtertyp==shrink)? "shrink filter" : "no filter ";
	std::cout<<filter_typ<<" applied "<<std::endl;
	if(m_filtertyp==bloom)
	{
		std::cout<<"effective false positive (bloom) "<<m_bloom_filter.effective_fpp()<<std::endl;
	}
	else if (m_filtertyp==shrink)
	{
		std::cout<<"shrink count (false negatives) :"<<m_shrink_cnt<<std::endl;
		std::cout<<"hash size per shrink threshold:"<<tune_shrink_threshold(m_ram)<<std::endl;
		double max_fn_rate = 0 ; 
		for (auto it: m_topNvector_zip)
		{
			max_fn_rate += (double)m_shrink_cnt/(double)it.second;

		}
		std::cout<<"max avg fn rate:"<<max_fn_rate/m_topNvector_zip.size()<<std::endl;

	}
	else
		std::cout<<"no filter applied"<<std::endl;

	std::cout<<"elapsed sec before destructor operation:"<<elapsed_secs<<std::endl;
}

inline void KMER_COUNTER::ShrinkHash(void)
{
	for (auto it=m_sequencehash_zip.begin();it!=m_sequencehash_zip.end();)
	{
		if(it->second==UNIQUE)
		{
			it = m_sequencehash_zip.erase(it);
		}
		else
		{
			it++;
		}
	}
			
	}

inline void KMER_COUNTER::Shrink_Insert(std::string& _line,unsigned int _limit)
{

	std::string subs;
	unsigned long long shrink_threshold=tune_shrink_threshold(m_ram);
	for ( unsigned int index=0; index<=_limit;  ++index)
	{
		subs=_line.substr(index,m_kmer_size);
		++m_sequencehash_zip[char2bit(subs, m_kmer_size)];
	}
	if( m_sequencehash_zip.size()>=shrink_threshold )
	{
		ShrinkHash();
		++m_shrink_cnt;
		shrink_threshold+=shrink_threshold;
	}
	
}
void KMER_COUNTER::BloomInit(void)
{
	bloom_parameters parameters;
	parameters.projected_element_count = m_estimated_insertions;
	parameters.false_positive_probability = 0.15;
	parameters.random_seed = 0xA5A5A5A5;
	parameters.compute_optimal_parameters();
	bloom_filter filter(parameters);
	m_bloom_filter=filter;
	return;
}
inline void KMER_COUNTER::Bloom_Insert(std::string& _line,unsigned int _limit)
{
	for ( unsigned int index=0; index<=_limit;  ++index)
	{
		auto subs=_line.substr(index,m_kmer_size);
		if (m_bloom_filter.contains(subs))
		{
			++m_sequencehash_zip[char2bit(subs, m_kmer_size)];
		}
		else
			m_bloom_filter.insert(subs);
	}
}
inline void KMER_COUNTER::Insert(std::string& _line, unsigned int _limit)
{
	std::string subs;
	for ( unsigned int index=0; index<=_limit;  ++index)
	{
		subs=_line.substr(index,m_kmer_size);
		++m_sequencehash_zip[char2bit(subs, m_kmer_size)];
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
