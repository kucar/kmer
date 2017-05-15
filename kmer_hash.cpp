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

//predefined hash map load factor threshold per shrink (if shrink filter is selected)
//this is to avoid high rate of collisions
loadfactor_thr_t loadfactor_thr={{minn,0.20},{low,0.25},{medium,0.30},{high,0.35},{maxx,0.40}};


inline uint32_t tune_shrink_threshold(int memusg)  {
    return     (memusg==maxx)    ?  THRES_PER_SHR<<4 :
           (memusg==high)    ?  THRES_PER_SHR<<3 :
           (memusg==medium)  ?  THRES_PER_SHR<<2 :
           (memusg==low)     ?  THRES_PER_SHR<<1 : THRES_PER_SHR ;
}

inline uint32_t tune_hashsize(int memusg)
{
    return  (memusg==maxx)   ?      HASHSIZE<<4 :
            (memusg==high)   ?  HASHSIZE<<3 :
            (memusg==medium) ?  HASHSIZE<<2 :
            (memusg==low)    ?  HASHSIZE<<1 : HASHSIZE ;
}
inline int getMappedCode(nucleotideCode_t n) {
    return indexMapper[n];
}

inline char getRevMappedCode(int n){
    return indexRevMapper[n];
}

inline u_64bits char2bit(const char * kmer, const uchar k,u_64bits& charsinbits )
{

    u_64bits bits = 0;
    bits |=getMappedCode((nucleotideCode_t)kmer[0]);

    for (uchar i = 1; i < k; i++) {
        bits = bits << 2;
        bits |=getMappedCode((nucleotideCode_t)kmer[i]);
    }
    charsinbits = bits;
    return charsinbits;
}

inline u_64bits char2bit(const char * kmer, const uchar k)
{
    u_64bits h = 0;
    return char2bit(kmer, k, h);
}

inline u_64bits char2bit(const std::string kmer, const uchar k)
{
    return char2bit(kmer.c_str(), k);
}

std::string bit2str(u_64bits charsinbits, uchar k)
{
    std::string retval = "";

    kmer_int val = charsinbits & 3;
    retval+=getRevMappedCode(val);

    for (uchar i = 1; i < k; i++) {
        charsinbits = charsinbits >> 2;
        val = charsinbits & 3;
        retval+=getRevMappedCode(val);
    }

    reverse(retval.begin(), retval.end());

    return retval;
}

KMER_COUNTER& KMER_COUNTER::operator=(KMER_COUNTER&& other)
{
    if(this!=&other)
    {
        m_fastq_file=other.m_fastq_file;
        m_topcount  =other.m_topcount;
        m_kmer_size =other.m_kmer_size;
        m_stats= other.m_stats;
        m_topNvector_zip=other.m_topNvector_zip;
        m_sequencehash_zip=other.m_sequencehash_zip;
        m_filtertyp=other.m_filtertyp;
        m_bloom_filter=other.m_bloom_filter;
        m_shrink_cnt=other.m_shrink_cnt;
        m_memusg = other.m_memusg;
        m_filter_ptr_arr[0]=&KMER_COUNTER::Shrink_Insert;
        m_filter_ptr_arr[1]=&KMER_COUNTER::Bloom_Insert;
        m_filter_ptr_arr[2]=&KMER_COUNTER::Insert;
        m_estimated_insertions=other.m_estimated_insertions;
        m_linelength=other.m_linelength;
        m_empty=false;
        other.ClearSequenceHash();
        other.ClearTopVector();
        other.m_empty=true;
    }
    return *this;
}

KMER_COUNTER::~KMER_COUNTER(){};

//process the file and call filter method
//and accumulate the hash.
void KMER_COUNTER::Init(void)
{
    uint32_t sz=tune_hashsize(m_memusg);
    m_sequencehash_zip.reserve((sz>m_estimated_insertions)? m_estimated_insertions : sz);
    std::string line;
    std::ifstream myfile (m_fastq_file);
    kmer_int linenum(0) ;
    auto limit= m_linelength - m_kmer_size;
    if(m_filtertyp==bloom)
        BloomInit();
    if (myfile.is_open())
    {
        while(getline(myfile,line))
            if (++linenum % NUMLINES_ENTRY == VALID_ENTRY) 
            {   //process the line, avoid cpu branching (if-else or switch-case)
                //instead use func ptr for direct access to the destined filter method
                (this->*m_filter_ptr_arr[m_filtertyp])(line,limit);
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
}
void KMER_COUNTER::ClearTopVector()
{
    m_topNvector_zip.clear();
}

void KMER_COUNTER::FindTopN()
{
#ifdef DEBUG_M
    volatile clock_t begin = clock();
#endif
    m_topNvector_zip.resize(m_topcount);
    std::partial_sort_copy(m_sequencehash_zip.begin(), m_sequencehash_zip.end(),m_topNvector_zip.begin(),m_topNvector_zip.end(),
            [](kmer_entry_zip const& l, kmer_entry_zip const& r)
                    {   return l.second > r.second; });
    for (auto it:m_topNvector_zip)
            std::cout<<bit2str(it.first,m_kmer_size)<<"---"<<it.second<<std::endl;
#ifdef DEBUG_M
    volatile clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
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
    std::cout<<"estimated # insertion:"<<m_estimated_insertions<<std::endl;
    std::string filter_typ=(m_filtertyp==bloom)? "bloom filter" : (m_filtertyp==shrink)? "shrink filter" : "no filter ";
    std::cout<<filter_typ<<" applied "<<std::endl;
    if(m_filtertyp==bloom)
    {
        std::cout<<"effective false positive rate (bloom) "<<m_bloom_filter.effective_fpp()<<std::endl;
    }
    else if (m_filtertyp==shrink)
    {
        std::cout<<"shrink count (false negatives) :"<<m_shrink_cnt<<std::endl;
        std::cout<<"hash size per shrink threshold:"<<tune_shrink_threshold(m_memusg)<<std::endl;
        double max_fn_rate = 0 ; 
        for (auto it: m_topNvector_zip)
        {
            max_fn_rate += (double)m_shrink_cnt/(double)it.second;

        }
        std::cout<<"avg fn rate(max):"<<max_fn_rate/m_topNvector_zip.size()<<std::endl;

    }
    else
        std::cout<<"no filter applied"<<std::endl;

    std::cout<<"elapsed sec before destructor operation:"<<elapsed_secs<<std::endl;
}
/*reduce the size of the hashmap by erasing the unique sequences*/
/*TODO: this can be openmp-ed*/
inline void KMER_COUNTER::ShrinkHash(void)
{
#ifdef DEBUG_M
    volatile clock_t begin = clock();
    unsigned int cnt =0;
#endif
    for (auto it=m_sequencehash_zip.begin();it!=m_sequencehash_zip.end();)
    {
        if(it->second==UNIQUE)
        {
            it = m_sequencehash_zip.erase(it);
#ifdef DEBUG_M
                        ++cnt;
#endif
        }
        else
        {
            it++;
        }
    }
#ifdef DEBUG_M
    volatile clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    std::cout<<"elapsed sec for shrinkHash():"<<elapsed_secs<<std::endl;
    std::cout<<"erased count                :"<<cnt<<std::endl;
#endif
            
    }

inline void KMER_COUNTER::Shrink_Insert(std::string& _line,unsigned int _limit)
{

    std::string subs;
    static unsigned long long shrink_threshold=tune_shrink_threshold(m_memusg);
    static float lf_threshold                 =loadfactor_thr[m_memusg];
    auto current_lf=m_sequencehash_zip.load_factor();
    for ( unsigned int index=0; index<=_limit;  ++index)
    {
        subs=_line.substr(index,m_kmer_size);
        ++m_sequencehash_zip[char2bit(subs, m_kmer_size)];
    }
    if( m_sequencehash_zip.size()>=shrink_threshold || current_lf>=lf_threshold )
    {
#ifdef DEBUG_M
        std::cout<<"shrink threshold = :"<<shrink_threshold<<std::endl;
        std::cout<<" current lf      = :"<<current_lf<<std::endl;
#endif
        ShrinkHash();
        ++m_shrink_cnt;
        shrink_threshold*=1.01;
        lf_threshold    *=1.01;
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
