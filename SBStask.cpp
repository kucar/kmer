/* Task : K-mer counting
 * Algorithm: each k-mer of size "kmersize" hashed into a STL hash map
 * Detail   : Input argument "filename" will be modified as only
 * 			  2nd lines of each entry will be placed on the file. File process
 * 			  handled by seperate sed command by system call.Then it is read sequentially
 *            line by line , each line will be processed against the
 *            hash map where kmer string is key and the count number is the value
 *            A hash wrapper object will be created at startup as a smart pointer.
 *            Then the rest of the main program will call init() and findTopN() methods
 *            of the object created.
 *TODO     : the program is open to refactoring based on
 *TODO       efficiency to speed up the process.
 *TODO        **eliminate substr() call - instead strview? could be used c++17 required
 *TODO        **implement a lock free multi threading mechanism ? like in jellyfish
 *TODO        **implement smart string searching algorithms ie  boyer moore  ?
 *TODO        **implement jenkins hash other than std::hash<string> ?
 *TODO     : implement a unit test mechanism like gmock
 *TODO     :imlement doxygen or so for documentation



 *Author : Korcan Ucar korcanucar@gmail.com
 *
 * */

#include "kmer_hash.h"
#include "file_op.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <memory>
#include <ctime>
#include <getopt.h>


using namespace std;

#define SUCCESS 0
#define PROCESS_ERROR 1

//globals for main program arguments

int kmersize(30);   				//default
int topcount(25);  		   		    //default
std::string filename="";
bool stats=false;			        //print stats
filter_t filter = shrink;		    //default
bool parallel=false;
int  ram(3);					    //default ram 3gb


void PrintUsage()
{
	cout<<"Usage: SBStask --filename <FILENAME> --kmersize <[1-31]]> --topcount <int> [--stats --[fp/nofilter] ]\n"<<endl;
	cout<<"--filename       :<string>path of the fastq file"<<endl;
	cout<<"--kmersize       :<integer>size of kmer [1-31]"<<endl;
	cout<<"--topcount       :<integer>top N kmers mostly observed"<<endl;
	cout<<"other options    :"<<endl;
	cout<<"--stats          : display information about hash infrastructure"<<endl;
	cout<<"--fp             : apply bloom filter"<<endl;
	cout<<"--no filter      : apply no filter"<<endl;

}

int ArgParse(int argc, char *argv[])
{
	int opt= 0;
	static struct option long_options[] = {
			{"filename", 		required_argument,  0,  'f' },
			{"kmersize", 		required_argument,  0,  'k' },
			{"topcount", 		required_argument,  0,  't' },
			{"stats", 			no_argument,	    0,  's' },
			{"linesizestatic",  no_argument,	    0,  'l' },
			{"help", 			no_argument,	    0,  'h' },
			{"fp",	 			no_argument,	    0,  'p' },
			{"fn", 				no_argument,	    0,  'n' },
			{"no filter", 		required_argument,	0,  'x' },
			{"ram", 			required_argument,	0,  'r' },
			{0,           		0,                  0,  0   }
	};

	int long_index =0;
	while ((opt = getopt_long(argc, argv,"ktsl:f:",
			long_options, &long_index )) != -1) {
		switch (opt) {
		case 'k' :
			kmersize = stoi(optarg);
			break;
		case 't' :
			topcount = stoi(optarg);
			break;
		case 'f' :
			filename = optarg;
			break;
		case 's' :
			stats = true;
			break;
		case 'p' :
			filter = bloom;
			break;
		case 'x':
			filter= nofilter;
			break;
		case 'r' :
			ram = stoi(optarg);
			break;
		default:
			PrintUsage();
			exit(PROCESS_ERROR);
		}
	}
	if (kmersize <= 0 || topcount <=0 ||filename=="") {
		PrintUsage();
		exit(PROCESS_ERROR);
	}

	return 0;
}


int main(int argc, char *argv[]) {

	ArgParse(argc,argv);
	unsigned int line_length;
	unsigned long long filesize=getFileSize(filename,line_length);
	unsigned long long linenum = estimate_line_num(filesize);
	unsigned long long est_inst= estimate_insertions(linenum,line_length,kmersize);
	std::shared_ptr<KMER_COUNTER> kmerobj (new KMER_COUNTER(filename,topcount,
															kmersize,stats,filter,
															ram,est_inst,line_length));
	volatile clock_t begin = clock();
	kmerobj->Begin();
	volatile clock_t end = clock();
	if(stats)
		kmerobj->PrintStats(begin,end);
	return SUCCESS;

}


