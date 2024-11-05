#include <iostream>
#include <string>
#include <cstring>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <sstream>
#include <map>
#include <cstdio>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>
#include <ctime>
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <inttypes.h>
#include <exception>
#include <unistd.h>
#include <thread>
// #include <boost/lockfree/spsc_queue.hpp>
// #include <boost/lockfree/queue.hpp>
// #include <boost/atomic.hpp>
#include "task_queue/concurrentqueue.h"
#include "task_queue/readerwriterqueue.h"
#include "SBT/Kmers.h"
#include "SBT/BF.h"
#include "SBT/BloomTree.h"
#include "SBT/util.h"
#include "SBT/Query.h"


#include "PthreadPool.h"
#include <queue>
#include "executor.h"
#include <sys/mman.h>
#include <fcntl.h>
#include <sys/stat.h> // struct stat 需要的头文件

#include <cstdio>

#define likely(x)       __builtin_expect(!!(x), 1)
#define unlikely(x)     __builtin_expect(!!(x), 0)

// #define ReadChunkSize 2048

// // #define ReadChunkSize 512
// #define InputQueueSize 2048
#define OutputQueueSize 2048
#define KmerSize 8
#define KmerPower 0x3FFF

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif




using namespace std;




typedef uint64_t bwtint_t;
typedef unsigned char ubyte_t;


//BWT-FM index
typedef struct {
	bwtint_t primary; // S^{-1}(0), or the primary index of BWT
	bwtint_t L2[5]; // C(), cumulative count
	bwtint_t seq_len; // sequence length
	bwtint_t bwt_size; // size of bwt, about seq_len/4
	uint32_t *bwt; // BWT
	uint32_t cnt_table[256];
	int sa_intv;
	bwtint_t n_sa;
	bwtint_t *sa;
} bwt_t;

typedef struct {
	int64_t offset;
	int32_t len;
	int32_t n_ambs;
	uint32_t gi;
	char *name, *anno;
} bntann1_t;

typedef struct {
	int64_t offset;
	int32_t len;
	char amb;
} bntamb1_t;

typedef struct {
	int64_t l_pac;
	int32_t n_seqs;
	uint32_t seed;
	bntann1_t *anns; // n_seqs elements
	int32_t n_holes;
	bntamb1_t *ambs; // n_holes elements
	FILE *fp_pac;
} bntseq_t;

typedef struct {
	bwtint_t x[3], info;
} bwtintv_t;

typedef struct {
	bwtintv_t * kmer_hash;
	int kmer_len;
	uint32_t kmer_num;
} bwa_hash_t;

typedef struct {
	bwt_t    *bwt; // FM-index
	bntseq_t *bns; // information on the reference sequences
	uint8_t  *pac; // the actual 2-bit encoded reference sequences with 'N' converted to a random base
	bwa_hash_t *hash;
} bwaidx_t;

// typedef struct
// {
// 	bwtint_t x[3];
// } bwtintv_t;


typedef struct
{
	uint32_t wid; // word id
	uint32_t pos; // occurrence position
} KmerItem_t;

typedef struct
{
	int PosDiff;
	uint32_t rPos;
	uint32_t gPos;
} KmerPair_t;

typedef struct
{
	bool bSimple;
	int rPos; // read position
	int64_t gPos; // genome position
	int rLen; // read block size
	int gLen; // genome block size
	int64_t PosDiff; // gPos-rPos
	int il;
	int interval;
} SeedPair_t;

typedef struct
{
	int Score; // seed score
	int64_t PosDiff;
	int64_t PosDiff_End;
	int PairedAlnCanIdx;
	vector<SeedPair_t> SeedVec;
} AlignmentCandidate_t;

typedef struct
{
	int len;
	int freq;
	bwtint_t* LocArr;
	bwtint_t il;
	int rightset;
	int leftest;
} bwtSearchResult_t;

typedef struct
{
	char* name; // chromosome name
	int64_t FowardLocation;
	int64_t ReverseLocation;
	int64_t len; // chromosome length
} Chromosome_t;

typedef struct
{
	bool bDir; // true:forward, false:reverse
	string CIGAR;
	int64_t gPos;
	int ChromosomeIdx;
} Coordinate_t;

typedef struct
{
	int AlnScore;
	int SamFlag; // sam flag
	int PairedAlnCanIdx;
	Coordinate_t coor;
	int NM ;
	int isPrime;
	int isSecondary;
} AlignmentReport_t;

typedef struct
{
	int rlen;


	char *header;
	char *qual;
	char *seq;

	char *info[4];
	//uint8_t* EncodeSeq;
	// aln report
	int mapq;
	int score;
	int sub_score;
	int CanNum;
	int iBestAlnCanIdx;
	AlignmentReport_t* AlnReportArr;
	int feature;
	int gotN;


	std::set<jellyfish::mer_dna> query_kmers;
    std::vector<const BloomTree*> matching;
    std::vector<float> weight;
} ReadItem_t;

typedef struct{

	ReadItem_t* readArr;

	vector<string> *SamOutputVec;
	vector<AlignmentCandidate_t>* AlignmentVec1;
	vector<SeedPair_t>* SeedPairVec1;
	uint8_t** EncodeSeq;
	int readNum;
	int index ;	
	int load;
	int myUniqueMapping, myUnMapping, myPaired,myDistance;
	int EstDistance;
}QueueItem_t;


// typedef struct{
// 	int index;
// 	int readNum;
// 	vector<string> *SamOutputVec;
// }SamOutputItem_t;
typedef struct {
	size_t offset;
	int sam_block_size;
	int localIndex;
	int readNum;
	int index;
	int myUniqueMapping, myUnMapping, myPaired,myDistance;
	char *samOutputBuffer;
	// vector<string> *SamOutputVec;
}SamOutputItem_t;


// typedef struct {
// 	bwtint_t x[3], info;
// } bwtintv_t;


typedef struct {

	int pen_unpaired;       // phred-scaled penalty for unpaired reads
	int zdrop;              // Z-dropoff

	int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
	int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
	float mask_level_redun;
	int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
	int a, b;               // match score and mismatch penalty
	int o_del, e_del;
	int o_ins, e_ins;
	int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
	int w;                  // band width
	int T;                  // output score threshold; only affecting output
	int seed_len;
	int flag;               // see MEM_F_* macros
	int min_seed_len;       // minimum seed length
	int minSeedLen;
	int errorAllow;
	int refineThreshold ;
	int descendingLength;
	int isGlobal;
	// int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
	int n_threads;          // number of threads
	int chunk_size;         // process chunk_size-bp sequences in a batch 批处理chunk_size-bp序列
	float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
	float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
	float mapQ_coef_len;
	int mapQ_coef_fac;
	int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
} opt_t;

// // 线程函数参数
// typedef struct  {
// 	opt_t *opt;
// 	int copy_comment, actual_chunk_size ,chunk_num;
	
// 	int n_seqs[InputQueueSize] ;
// } threadInputArg;


// 线程函数参数
typedef struct  {
	opt_t *opt;
	int copy_comment, actual_chunk_size,chunk_num;
	int n_seqs ;
	QueueItem_t *temp;
} threadOutputArg;


// 线程函数参数
typedef struct  {
	opt_t *opt;
	int threadId;
	QueueItem_t *temp;
} threadMappingArg;

// 线程信息
typedef struct {
    pthread_t thread_id;
    pthread_mutex_t ioLock;
	pthread_mutex_t runLock;

    pthread_cond_t ioCond;
	pthread_cond_t runCond;

    int runFlag;
	int ioFlag;
    int buffer_flag;
    int shutdown;

	int head;
	int tail;
	int count;
	int full;
	int next;
	int queue_size;
	QueueItem_t * temp;

} thread_info;




typedef struct {
	struct ktp_t *pl;
	int64_t index;
	int step;
	void *data;
} ktp_worker_t;


typedef struct {
	size_t offset1;
	int pos1;
	char *tmp1;
	int right_broundry1;
	int is_end1=0;
} read_block_arg;




struct cmp1   
{
	bool operator() ( QueueItem_t* a, QueueItem_t *b) {
		return a->index > b->index;
	}
};

struct cmp_index   
{
	bool operator() ( SamOutputItem_t* a, SamOutputItem_t *b) {
		return a->index > b->index;
	}
};

// Global variables
extern bwt_t *bwt;
extern bwaidx_t *bwtIdx;
extern string indexPath,readFile1 ,readFile2;
extern unsigned char nst_nt4_table[256];
extern int64_t GenomeSize, TwoGenomeSize;
extern vector<Chromosome_t> ChromosomeVec;
extern vector<string> ReadVec, ReadHeaderVec;
extern vector<int64_t> AccumulationLengthVec, PositionShiftPosVec;

extern int InputQueueSize,ReadChunkSize;
extern const char* VersionStr;
extern map<int64_t, int> ChrLocMap;
extern vector<string> ReadFileNameVec1, ReadFileNameVec2;
extern char *refSeq, *GenomeFileName, *IndexFileName, *OutputFileName,*outputFile ,*indexFile ,*refSeq;;
extern bool bDebugMode, bPairEnd, bPacBioData, gzCompressed, FastQFormat, bMultiHit, bSilent ,bHash,bBlockRead;
extern int maxInsertSize, threadNum, iChromsomeNum, WholeChromosomeNum, ChromosomeNumMinusOne, MaxGaps, MinSeedLength, OutputFileFormat,OutputByOrder;

extern thread_info input_info, output_info, cal_input_info, cal_output_info;
// extern ReadItem_t *ReadArr[InputQueueSize];


// extern threadInputArg *input_arg;
// extern threadOutputArg *output_arg; 
extern threadMappingArg *mappingArg;

extern QueueItem_t **queuePointerForThreads;
extern QueueItem_t *tempQueueItem;

extern PthreadPool pool;


extern size_t fastq_file_size1 ,fastq_file_size2;
extern char* bigBlock1,*bigBlock2;


	/*初试化内存块*/
extern	size_t blocks;
extern	size_t block_size ; // 4K一个块
extern	size_t actual_block_size ;



extern size_t offset1 ,offset2 ; 
extern int pos1,pos2;
extern char *tmp1 ;
extern char *tmp2 ;
extern int right_broundry1,right_broundry2;

extern int is_end1,is_end2;

// extern FILE *rfile_text;
// extern SamOutputItem_t *SamOutputVecBuffer;

// extern priority_queue<QueueItem_t *, vector<QueueItem_t *>, cmp1 > priQueMinFirst;
extern priority_queue<SamOutputItem_t*, vector<SamOutputItem_t*>, cmp_index > priQueMinFirst;


extern moodycamel::ConcurrentQueue<char*> address_array; //用于回收数组

// bwt_fm_index.cpp
extern void RestoreReferenceInfo();
extern void bwa_idx_destroy(bwaidx_t *idx);
extern bwaidx_t *bwa_idx_load(const char *hint);
extern bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop);
extern bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop);
extern bwtSearchResult_t BWT_Search_Forward(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);
extern bwtSearchResult_t BWT_Search_Forward_1(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);
extern bwtSearchResult_t BWT_Search_Forward_2(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);
extern bwtSearchResult_t BWT_Search_Forward_3(uint8_t* seq, int start, int stop,int last_rightest ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);

extern bwtSearchResult_t BWT_Search_Forward_hash(uint8_t* seq, int start, int stop ,int last_rightset,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);

extern void bwt_sa_batch(bwtint_t k,int interval,int match_len,int match_beg,vector<SeedPair_t>& SeedPairVec);
extern bwtSearchResult_t BWT_Only_Search(uint8_t* seq, int start, int stop);
extern bwtSearchResult_t BWT_Search_Backward(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);

extern bwtSearchResult_t BWT_Search_Backward_hash(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);

extern bwtSearchResult_t BWT_Search_BothSide(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);
extern bwtSearchResult_t BWT_Search_BothSide1(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list);
extern int GetNextChunk_Block(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2,size_t& offset1,size_t& offset2,int &is_end,char * resident_content);
extern int GetNextChunk_Block_1(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2,size_t& offset1,size_t& offset2,char *tmp1,char *tmp2,int &pos1,int &pos2,int & right_broundary1,int &right_broundry2 ,int &is_end1,int &is_end2);

extern int GetNextChunk_Block_2(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2);

extern bwtint_t bwt_sa(bwtint_t k);

// GetData.cpp
extern bool CheckReadFormat(const char* filename);
extern bool CheckBWAIndexFiles(string IndexPrefix);
// extern int GetNextChunk(bool bSepLibrary, FILE *file, FILE *file2, ReadItem_t* ReadArr);

// AlignmentRescue.cpp
extern bool RescueUnpairedAlignment(int EstDistance, ReadItem_t& r1, ReadItem_t& r2, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2);


extern int GetNextChunk(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2);

extern int GetNextChunk_Recycle(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2);

extern int IdentifyHeaderBegPos(char* str, int len);
extern int IdentifyHeaderEndPos(char* str, int len);
extern int gzGetNextChunk_block(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr,int chunkSize, int flag);
extern int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr,int chunkSize, int flag);


// Mapping.cpp
extern void Mapping();
extern void process();
// extern bool CheckPairedAlignmentCandidates(int64_t EstiDistance, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2);


// AlignmentCandidates.cpp
extern void RemoveShortSeeds(vector<SeedPair_t>& SeedVec, int thr);
extern bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2);
extern bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2);
extern void IdentifyNormalPairs(int rlen, int glen, vector<SeedPair_t>& SeedVec);
extern vector<SeedPair_t> IdentifySeedPairs_FastMode(int rlen, uint8_t* EncodeSeq);
extern vector<SeedPair_t> IdentifySeedPairs_RescueMode(int rlen, uint8_t* EncodeSeq);
extern vector<SeedPair_t> IdentifySeedPairs_SensitiveMode(int rlen, uint8_t* EncodeSeq);
extern void IdentifySeedPairs_FastMode_P(int rlen, uint8_t* EncodeSeq,vector<SeedPair_t> *SeedPairVec,int rlen_p, uint8_t* EncodeSeq_p,vector<SeedPair_t> *SeedPairVec_p);

extern void GenMappingReport(bool bFirstRead, ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec);
extern vector<AlignmentCandidate_t> GenerateAlignmentCandidateForPacBioSeq(int rlen, vector<SeedPair_t> SeedPairVec);
extern vector<AlignmentCandidate_t> GenerateAlignmentCandidateForIlluminaSeq(int rlen, vector<SeedPair_t> SeedPairVec);


extern void IdentifySeedPairs_FastMode_Recycle(int rlen, uint8_t* EncodeSeq ,vector<SeedPair_t>& SeedPairVec);

extern void IdentifySeedPairs_FastMode_getN(int rlen, uint8_t* EncodeSeq ,vector<SeedPair_t>& SeedPairVec ,int gotN);

extern void GenerateAlignmentCandidateForIlluminaSeq_Recycle(int rlen, vector<SeedPair_t>& SeedPairVec,vector<AlignmentCandidate_t>& AlignmentVec);

extern Coordinate_t GenCoordinateInfo(bool bFirstRead, int64_t gPos, int64_t end_gPos, vector<pair<int, char> >& cigar_vec);


// tools.cpp
extern void ShowSeedLocationInfo(int64_t MyPos);
extern int64_t GetAlignmentBoundary(int64_t gPos);
extern bool CheckFragValidity(SeedPair_t SeedPair);
extern void SelfComplementarySeq(int len, char* rseq);
extern void ShowSeedInfo(vector<SeedPair_t>& SeedPairVec);
extern void GetComplementarySeq(int len, char* seq, char* rseq);
extern bool CheckAlignmentValidity(vector<SeedPair_t>& SeedPairVec);
extern int CalFragPairIdenticalBases(int len, char* frag1, char* frag2);
extern int AddNewCigarElements(string& str1, string& str2, vector<pair<int, char> >& cigar_vec);
extern int ProcessHeadSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern int ProcessTailSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern int ProcessNormalSequencePair(char* seq, SeedPair_t& sp, vector<pair<int, char> >& cigar_vec);
extern void ShowAlignmentCandidateInfo(bool bFirst, char* header, vector<AlignmentCandidate_t>& AlignmentVec);


// KmerAnalysis.cpp
extern vector<KmerItem_t> CreateKmerVecFromReadSeq(int len, char* seq);
extern vector<KmerPair_t> IdentifyCommonKmers(int MaxShift, vector<KmerItem_t>& vec1, vector<KmerItem_t>& vec2);
extern vector<SeedPair_t> GenerateSimplePairsFromCommonKmers(int MinSeedLength, vector<KmerPair_t>& KmerPairVec);
extern vector<SeedPair_t> GenerateSimplePairsFromFragmentPair(int MaxDist, int len1, char* frag1, int len2, char* frag2);

// nw_alignment.cpp
extern void nw_alignment(int m, string& s1, int n, string& s2);

extern void nw_alignment_band(int m, string& s1, int n, string& s2);



