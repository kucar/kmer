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

int kmersize=30;   				//default
int topcount=25;  		   		//default
std::string filename="";
bool stats=false;			   //print stats
bool bloom= false;
bool parallel=false;


void PrintUsage()
{
	cout<<"Usage: SBStask --filename <FILENAME> --kmersize <KMERSZ> --topcount <TOPN> [--stats --linesizestatic]\n"<<endl;

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
			{"multi", 			no_argument,	    0,  'm' },
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
			bloom = true;
			break;
		case 'm' :
			parallel = true;
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
	std::shared_ptr<KMER_BASE> kmerobj (new KMER_BASE(filename,topcount,kmersize,stats,bloom,parallel));
	volatile clock_t begin = clock();
	kmerobj->Begin();
	volatile clock_t end = clock();
	if(stats)
		kmerobj->PrintStats(begin,end);
	return SUCCESS;

}


