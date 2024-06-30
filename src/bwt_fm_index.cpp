#include "util.h"



#define OCC_Thr 50
#define OCC_INTV_SHIFT 7
#define OCC_INTERVAL   (1LL<<OCC_INTV_SHIFT)
#define OCC_INTV_MASK  (OCC_INTERVAL - 1)

#define __occ_aux4(bwt, b)											\
	((bwt)->cnt_table[(b)&0xff] + (bwt)->cnt_table[(b)>>8&0xff]		\
	 + (bwt)->cnt_table[(b)>>16&0xff] + (bwt)->cnt_table[(b)>>24])

#define bwt_occ_intv(b, k) ((b)->bwt + ((k)>>7<<4))
#define bwt_bwt(b, k) ((b)->bwt[((k)>>7<<4) + sizeof(bwtint_t) + (((k)&0x7f)>>4)])
#define bwt_B0(b, k) (bwt_bwt(b, k)>>((~(k)&0xf)<<1)&3)

#define bwt_set_intv(bwt, c, ik) ((ik).x[0] = (bwt)->L2[(int)(c)]+1, (ik).x[2] = (bwt)->L2[(int)(c)+1]-(bwt)->L2[(int)(c)], (ik).x[1] = (bwt)->L2[3-(c)]+1, (ik).info = 0)


int iChromsomeNum;
map<int64_t, int> ChrLocMap;
int64_t GenomeSize, TwoGenomeSize;
vector<Chromosome_t> ChromosomeVec;




void *IdvLoadReferenceSequences(void *arg)
{
	int base, *my_id;
	int64_t fPos, rPos;

	my_id = (int*)arg;
	for (fPos = *my_id, rPos = TwoGenomeSize - fPos - 1 ; fPos < GenomeSize; fPos += threadNum, rPos-= threadNum)
	{
		base = bwtIdx->pac[fPos >> 2] >> ((~fPos & 3) << 1) & 3;
		switch (base)
		{
		case 0: refSeq[fPos] = 'A'; refSeq[rPos] = 'T'; break;
		case 1: refSeq[fPos] = 'C'; refSeq[rPos] = 'G'; break;
		case 2: refSeq[fPos] = 'G'; refSeq[rPos] = 'C'; break;
		case 3: refSeq[fPos] = 'T'; refSeq[rPos] = 'A'; break;
		default:refSeq[fPos] = refSeq[rPos] = 'N';
		}
	}
	return (void*)(1);
}


void RestoreReferenceSequences()
{
	
	int i, *JobIDArr = new int[threadNum];

	pthread_t *ThreadArr = new pthread_t[threadNum];
	for (i = 0; i < threadNum; i++)
	{
		JobIDArr[i] = i;
		pthread_create(&ThreadArr[i], NULL, IdvLoadReferenceSequences, JobIDArr + i);
	}
	for (i = 0; i < threadNum; i++) pthread_join(ThreadArr[i], NULL);

	delete[] ThreadArr; delete[] JobIDArr;
}


void RestoreReferenceInfo()
{
	int i;
	int64_t iTotalLength = 0;

	GenomeSize = bwtIdx->bns->l_pac; TwoGenomeSize = (GenomeSize << 1);
	iChromsomeNum = bwtIdx->bns->n_seqs; ChromosomeVec.resize(iChromsomeNum);

	fprintf(stderr, "Load the reference sequences...\n");
	fseek(bwtIdx->bns->fp_pac, 0, SEEK_SET);
	(void)fread(bwtIdx->pac, 1, GenomeSize / 4 + 1, bwtIdx->bns->fp_pac);

	for (i = 0; i < iChromsomeNum; i++)
	{
		ChromosomeVec[i].len = bwtIdx->bns->anns[i].len;
		ChromosomeVec[i].name = bwtIdx->bns->anns[i].name;

		ChromosomeVec[i].FowardLocation = iTotalLength; iTotalLength += ChromosomeVec[i].len;
		ChromosomeVec[i].ReverseLocation = TwoGenomeSize - iTotalLength;

		ChrLocMap.insert(make_pair(ChromosomeVec[i].FowardLocation + ChromosomeVec[i].len - 1, i));
		ChrLocMap.insert(make_pair(ChromosomeVec[i].ReverseLocation + ChromosomeVec[i].len - 1, i));

		// fprintf(stdout,"ChromosomeVec[%d]:%s  len =%d FowardLocation:%ld ReverseLocation:%ldn",i,ChromosomeVec[i].name,ChromosomeVec[i].len,ChromosomeVec[i].FowardLocation,ChromosomeVec[i].ReverseLocation);
	}
	refSeq = new char[TwoGenomeSize + 1]; refSeq[TwoGenomeSize] = '\0';
	RestoreReferenceSequences();
	//fprintf(stdout, "\n");
	//for (map<int64_t, int>::iterator iter = ChrLocMap.begin(); iter != ChrLocMap.end(); iter++) printf("chr%d: %lld\n", iter->second, iter->first);
	//for (map<int64_t, int>::iterator iter = ChrLocMap.begin(); iter != ChrLocMap.end(); iter++) printf("Chr: %s [%ld -- %ld]\n", ChromosomeVec[iter->second].name, iter->first - ChromosomeVec[iter->second].len + 1, iter->first);
	// if (bwtIdx->bns->fp_pac){
			
	// } 
	fclose(bwtIdx->bns->fp_pac);
	fprintf(stderr, "Finish loading the reference sequences...\n");
}
void bns_destroy(bntseq_t *bns)
{
	if (bns == 0) return;
	else {
		int i;
		// if (bns->fp_pac){
		// 	fclose(bns->fp_pac);
		// } 
		free(bns->ambs);
		for (i = 0; i < bns->n_seqs; ++i) {
			free(bns->anns[i].name);
			free(bns->anns[i].anno);
		}
		free(bns->anns);
		free(bns);
	}
}

void bwt_destroy(bwt_t *bwt)
{
	if (bwt == 0) return;
	free(bwt->sa); free(bwt->bwt);
	free(bwt);
}

void bwt_hash_destroy(bwa_hash_t * hash)
{
	if (bwt == 0) return;
	delete[] (hash->kmer_hash);
	free(hash);
}

void bwa_idx_destroy(bwaidx_t *idx)
{
	if (idx == 0) return;
	if (idx->bwt) bwt_destroy(idx->bwt);
	if (idx->bns) bns_destroy(idx->bns);
	if (idx->pac) free(idx->pac);
	if (idx->hash) bwt_hash_destroy(idx->hash);
	free(idx);
}

bwtint_t fread_fix(FILE *fp, bwtint_t size, void *a)
{ // Mac/Darwin has a bug when reading data longer than 2GB. This function fixes this issue by reading data in small chunks
	const int bufsize = 0x1000000; // 16M block
	bwtint_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = fread((char*)a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	return offset;
}

void bwt_gen_cnt_table(bwt_t *bwt)
{
	int i, j;
	for (i = 0; i != 256; ++i) {
		uint32_t x = 0;
		for (j = 0; j != 4; ++j)
			x |= (((i&3) == j) + ((i>>2&3) == j) + ((i>>4&3) == j) + (i>>6 == j)) << (j<<3);
		bwt->cnt_table[i] = x;
	}
}

bwt_t *bwt_restore_bwt(const char *fn)
{
	bwt_t *bwt;
	FILE *fp;

	bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
	fp = fopen(fn, "rb");
	fseek(fp, 0, SEEK_END);
	bwt->bwt_size = (ftell(fp) - sizeof(bwtint_t) * 5) >> 2;
	bwt->bwt = (uint32_t*)calloc(bwt->bwt_size, 4);
	fseek(fp, 0, SEEK_SET);
	fread(&bwt->primary, sizeof(bwtint_t), 1, fp);
	fread(bwt->L2 + 1, sizeof(bwtint_t), 4, fp);
	fread_fix(fp, bwt->bwt_size << 2, bwt->bwt);
	bwt->seq_len = bwt->L2[4];
	fclose(fp);
	bwt_gen_cnt_table(bwt);

	return bwt;
}

bwa_hash_t *bwt_restore_hash(const char *fn)
{

	bwa_hash_t *hash;
	FILE *fp;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;
	char *d="\t"; 
	hash = (bwa_hash_t*)calloc(1, sizeof(bwa_hash_t));
	fp = fopen(fn, "rb");

	if((len = getline(&buffer, &size, fp)) != -1){
		hash->kmer_len =atoi(strtok(buffer,d));
		hash->kmer_num =atoi(strtok(NULL,d));
		hash->kmer_hash = new bwtintv_t[hash->kmer_num];
		// printf("%s\n",buffer);
		// printf("p:%d\n",hash->kmer_len);
		// printf("p:%d\n",hash->kmer_num);
	}

	while((len = getline(&buffer, &size, fp)) != -1){
		uint32_t kmer_no =atoi(strtok(buffer,d));
		hash->kmer_hash[kmer_no].x[0] =atol(strtok(NULL,d));
		hash->kmer_hash[kmer_no].x[1] =atol(strtok(NULL,d));
		hash->kmer_hash[kmer_no].x[2] =atol(strtok(NULL,d));
		// printf("kmer_no:%d x0:%ld,x1:%ld,x2:%ld\n",kmer_no,hash->kmer_hash[kmer_no].x[0],hash->kmer_hash[kmer_no].x[1],hash->kmer_hash[kmer_no].x[2]);
	}
	fclose(fp);

	// for(int i=0;i<hash->kmer_num;i++){
	// 	printf("kmer_no:%d x0:%ld,x1:%ld,x2:%ld\n",i,hash->kmer_hash[i].x[0],hash->kmer_hash[i].x[1],hash->kmer_hash[i].x[2]);

	// }


	return hash;
}

void bwt_restore_sa(const char *fn, bwt_t *bwt)
{
	char skipped[256];
	FILE *fp;
	bwtint_t primary;

	fp = fopen(fn, "rb");
	fread(&primary, sizeof(bwtint_t), 1, fp);
	// xassert(primary == bwt->primary, "SA-BWT inconsistency: primary is not the same.");
	fread(skipped, sizeof(bwtint_t), 4, fp); // skip
	fread(&bwt->sa_intv, sizeof(bwtint_t), 1, fp);
	fread(&primary, sizeof(bwtint_t), 1, fp);
	// xassert(primary == bwt->seq_len, "SA-BWT inconsistency: seq_len is not the same.");

	bwt->n_sa = (bwt->seq_len + bwt->sa_intv) / bwt->sa_intv;
	bwt->sa = (bwtint_t*)calloc(bwt->n_sa, sizeof(bwtint_t));
	bwt->sa[0] = -1;

	fread_fix(fp, sizeof(bwtint_t) * (bwt->n_sa - 1), bwt->sa + 1);
	fclose(fp);
}


bwt_t *bwa_idx_load_bwt(const char *prefix)
{
	char *tmp ;
	bwt_t *bwt;

	tmp = (char*)calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".bwt"); // FM-index
	bwt = bwt_restore_bwt(tmp);
	strcat(strcpy(tmp, prefix), ".sa");  // partial suffix array (SA)
	bwt_restore_sa(tmp, bwt);
	free(tmp); 
	return bwt;
}


bwa_hash_t *bwa_idx_load_hash(const char *prefix)
{
	char *tmp ;
	bwa_hash_t *hash;

	tmp = (char*)calloc(strlen(prefix) + 5, 1);
	strcat(strcpy(tmp, prefix), ".hsh"); // Hash of FM-index
	hash = bwt_restore_hash(tmp);

	free(tmp); 
	return hash;
}


void err_fatal(const char *header, const char *fmt, ...)
{
	va_list args;
	va_start(args, fmt);
	fprintf(stderr, "[%s] ", header);
	vfprintf(stderr, fmt, args);
	fprintf(stderr, "n");
	va_end(args);
	exit(EXIT_FAILURE);
}


bntseq_t *bns_restore_core(const char *ann_filename, const char* amb_filename, const char* pac_filename)
{
	char str[8192];
	FILE *fp;
	const char *fname;
	bntseq_t *bns;
	long long xx;
	int i;
	int scanres;
	bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
	{ // read .ann
		fp = fopen(fname = ann_filename, "r");
		scanres = fscanf(fp, "%lld%d%u", &xx, &bns->n_seqs, &bns->seed);
		// printf("%lld -%d -%u ,scanres:%d\n",xx, bns->n_seqs, bns->seed,scanres);
		if (scanres != 3) goto badread;
		bns->l_pac = xx;
		bns->anns = (bntann1_t*)calloc(bns->n_seqs, sizeof(bntann1_t));
		for (i = 0; i < bns->n_seqs; ++i) {
			bntann1_t *p = bns->anns + i;
			char *q = str;
			int c;
			// read gi and sequence name
			scanres = fscanf(fp, "%u%s", &p->gi, str);
			// printf("%u - %s scanres:%d\n",p->gi, str,scanres);
			if (scanres != 2) goto badread;
			p->name = strdup(str);
			// read fasta comments 
			while (q - str < sizeof(str) - 1 && (c = fgetc(fp)) != '\n' && c != EOF) *q++ = c;
			while (c != '\n' && c != EOF) c = fgetc(fp);
			if (c == EOF) {
				scanres = EOF;
				// printf("scanres:%d\n",scanres);
				goto badread;
			}
			*q = 0;
			if (q - str > 1 && strcmp(str, " (null)") != 0) p->anno = strdup(str + 1); // skip leading space
			else p->anno = strdup("");
			// read the rest
			scanres = fscanf(fp, "%lld%d%d", &xx, &p->len, &p->n_ambs);
			// printf("scanres:%d\n",scanres);
			if (scanres != 3) goto badread;
			p->offset = xx;
		}
		fclose(fp);
	}
	{ // read .amb

		int32_t n_seqs;
		fp = fopen(fname = amb_filename, "r");
		scanres = fscanf(fp, "%lld%d%d", &xx, &n_seqs, &bns->n_holes);
		if (scanres != 3) goto badread;

		// xassert(l_pac == bns->l_pac && n_seqs == bns->n_seqs, "inconsistent .ann and .amb files.");
		bns->ambs = bns->n_holes? (bntamb1_t*)calloc(bns->n_holes, sizeof(bntamb1_t)) : 0;
		for (i = 0; i < bns->n_holes; ++i) {
			bntamb1_t *p = bns->ambs + i;
			scanres = fscanf(fp, "%lld%d%s", &xx, &p->len, str);
			if (scanres != 3) goto badread;
			p->offset = xx;
			p->amb = str[0];
		}
		fclose(fp);
	}
	{ // open .pac
		bns->fp_pac = fopen(pac_filename, "rb");
	}
	return bns;

 badread:
	if (EOF == scanres) {
		err_fatal(__func__, "Error reading %s : %sn", fname, ferror(fp) ? strerror(errno) : "Unexpected end of file");
	}
	err_fatal(__func__, "Parse error reading %sn", fname);
}



bntseq_t *bns_restore(const char *prefix)
{  
	char ann_filename[256], amb_filename[256], pac_filename[256];
	bntseq_t *bns;
	strcat(strcpy(ann_filename, prefix), ".ann");
	strcat(strcpy(amb_filename, prefix), ".amb");
	strcat(strcpy(pac_filename, prefix), ".pac");



	bns =bns_restore_core(ann_filename, amb_filename, pac_filename);
	if (bns == 0){
		return 0;
	}
		// if ((fp = fopen(strcat(strcpy(alt_filename, prefix), ".alt"), "r")) != 0) {} // read .alt file if present 
	return bns;
}

bwaidx_t *bwa_idx_load(const char *hint)
{
	bwaidx_t *idx;

	fprintf(stderr, "Load the genome index files...");
	idx = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));
	idx->bwt = bwa_idx_load_bwt(hint);
	idx->bns = bns_restore(hint);
	idx->pac = (uint8_t*)calloc(idx->bns->l_pac/4+1, 1);
	idx->hash = bwa_idx_load_hash(hint);
	fprintf(stderr, "\n");

	return idx;
}


bool CheckReadFormat(const char* filename)
{
	char buf[1];
	gzFile file = gzopen(filename, "rb");
	gzread(file, buf, 1); gzclose(file);

	if (buf[0] == '@') return true; // fastq
	else return false;
}

bool CheckBWAIndexFiles(string IndexPrefix)
{
	fstream file;
	string filename;
	bool bChecked = true;

	// filename = IndexPrefix + ".ann"; file.open(filename.c_str(), ios_base::in);
	// if (!file.is_open()) return false; else file.close();

	// filename = IndexPrefix + ".amb"; file.open(filename.c_str(), ios_base::in);
	// if (!file.is_open()) return false; else file.close();

	// filename = IndexPrefix + ".pac"; file.open(filename.c_str(), ios_base::in);
	// if (!file.is_open()) return false; else file.close();

	return bChecked;
}



int IdentifyHeaderBegPos(char* str, int len)
{
	int i;

	for (i = 1; i < len; i++)
	{
		if (str[i] != '>' && str[i] != '@') return i;
	}
	return len - 1;
}

int IdentifyHeaderBegPos_1(char* str, int len)
{
	int i;

	for (i = 1; i < len; i++)
	{
		if (*(str+i) != '>' && *(str+i) != '@') return i;
	}
	return len - 1;
}



int IdentifyHeaderEndPos_1(char* str, int len){
	int i;

	for (i = 1; i < len; i++)
	{
		if (*(str+i) == ' ' || *(str+i) == '/' || *(str+i) == '\t') return i;
	}
	return len - 1;
}

int IdentifyHeaderEndPos(char* str, int len){
	int i;

	for (i = 1; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '/' || str[i] == '\t') return i;
	}
	return len - 1;
}

int IdentifyHeaderBegPos(string &str, int len){
	int i;

	for (i = 1; i < len; i++)
	{
		if (str[i] != '>' && str[i] != '@') return i;
	}
	return len - 1;
}

int IdentifyHeaderEndPos(string &str, int len){
	int i;

	for (i = 1; i < len; i++)
	{
		if (str[i] == ' ' || str[i] == '/' || str[i] == '\t') return i;
	}
	return len - 1;
}



char GetComplementaryBase(char c)
{
	switch (c)
	{
	case 'A': return 'T';
	case 'a': return 'T';
	case 'C': return 'G';
	case 'c': return 'G';
	case 'G': return 'C';
	case 'g': return 'C';
	case 'T': return 'A';
	case 't': return 'A';
	default:  return 'N';
	}
}

void GetComplementarySeq(int len, char* seq, char* rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i<j; i++, j--){
		rseq[i] = GetComplementaryBase(seq[j]);
		rseq[j] = GetComplementaryBase(seq[i]);
	}
	if (i == j) rseq[i] = GetComplementaryBase(seq[i]);
}

void GetComplementarySeq1(int len,const char* seq, char* rseq)
{
	int i, j;

	for (j = len - 1, i = 0; i<j; i++, j--){
		rseq[i] = GetComplementaryBase(seq[j]);
		rseq[j] = GetComplementaryBase(seq[i]);
	}
	if (i == j) rseq[i] = GetComplementaryBase(seq[i]);
}

void GetComplementarySeq2(int len, char* seq)
{
	int i, j;
	char tmpi,tmpj;
	for (j = len - 1, i = 0; i<j; i++, j--){
		tmpi =seq[i];
		tmpj =seq[j];
		seq[i] = GetComplementaryBase(tmpj);
		seq[j] = GetComplementaryBase(tmpi);
	}
	if (i == j) seq[i] = GetComplementaryBase(seq[i]);
}

void GetReverseSeq(int len, char* seq){
	int i, j;
	char tmpi,tmpj;
	for (j = len - 1, i = 0; i<j; i++, j--){
		tmpi =seq[i];
		seq[i] = seq[j];
		seq[j] = tmpi;
	}
	if (i == j) seq[i] = (seq[i]);
}

ReadItem_t GetNextEntry(FILE *file)
{
	int p1, p2;
	size_t len;
	size_t size = 0;
	ReadItem_t read;
	char *buffer = NULL;

	read.header = read.seq = read.qual = NULL; read.rlen = 0;

	if ((len = getline(&buffer, &size, file)) != -1)
	{
		p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, (buffer + 1), len); read.header[len] = '0';

		if (FastQFormat)
		{
			if ((read.rlen = getline(&buffer, &size, file)) != -1)
			{
				read.seq = new char[read.rlen];
				strncpy(read.seq, buffer, read.rlen);
				getline(&buffer, &size, file); getline(&buffer, &size, file);
				read.qual = new char[read.rlen]; strncpy(read.qual, buffer, read.rlen);
				read.rlen -= 1; read.seq[read.rlen] = '\0'; read.qual[read.rlen] = '\0';
			}
			else read.rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read.rlen = (int)seq.length()) > 0)
			{
				read.seq = new char[read.rlen + 1];
				strcpy(read.seq, (char*)seq.c_str());
				read.seq[read.rlen] = '\0';
			}
		}
	}
	free(buffer);

	return read;
}


size_t getline_block_buffer(FILE *file ,size_t file_size ,size_t& offset,char *& buffer,char *tmp, int &pos ,int & right_broundary,int & is_end){

	int last =pos;
	size_t len =0;
	size_t ret =0;

	while(tmp[pos]!='\n' && pos<right_broundary){
		pos++;
	}
	len =pos-last;


	if(tmp[pos]=='\n'){
		pos++;
		buffer =(tmp+last);
		return len;
	}else{
		if(is_end){
			return -1;
		}
		strncpy(tmp, (tmp+last), len);
		tmp[len]='\0';

		if(offset + actual_block_size >=file_size){

			ret=fread(tmp+len,sizeof(char), actual_block_size, file);
			tmp[len+(file_size - offset)]='\0';
			right_broundary =len+(file_size - offset);
			is_end =1;
			// printf("tmp:%s\n",tmp);
		}else{
			ret=fread(tmp+len, sizeof(char),actual_block_size,  file);
			tmp[len+actual_block_size]='\0';
			right_broundary =len+actual_block_size;
			offset+=actual_block_size;

		}

		pos =len;
		last=0;
		while(tmp[pos]!='\n'){
			pos++;
		}
		len =pos-last;
		pos++;
	}

	buffer =(tmp+last);

	return len;	

}
size_t getline_block_buffer_2(FILE *file  ,char *& buffer){

	int last =pos1;
	size_t len =0;
	size_t ret =0;
	int is_eof;

	while(tmp1[pos1]!='\n' && pos1<right_broundry1){
		pos1++;
	}
	len =pos1-last;


	if(tmp1[pos1]=='\n'){
		pos1++;
		buffer =(tmp1+last);
		return len;
	}else{
		if(is_end1){
			return -1;
		}
		strncpy(tmp1, (tmp1+last), len);
		tmp1[len]='\0';

		// if(offset1 + actual_block_size >=fastq_file_size1){

		// 	ret=fread(tmp1+len,sizeof(char), actual_block_size, file);
		// 	tmp1[len+(fastq_file_size1 - offset1)]='\0';
		// 	right_broundry1 =len+(fastq_file_size1 - offset1);
		// 	is_end1 =1;
		// 	// printf("tmp:%s\n",tmp);
		// }else{
		// 	ret=fread(tmp1+len, sizeof(char),actual_block_size,  file);
		// 	tmp1[len+actual_block_size]='\0';
		// 	right_broundry1 =len+actual_block_size;
		// 	offset1+=actual_block_size;
		// 		// __builtin_prefetch (tmp1, 0);

		// }

		ret=fread(tmp1+len,sizeof(char), actual_block_size, file);
		
		// ret =read(fileno(file),tmp1+len, sizeof(char)* actual_block_size);

		// ret=gzread(file,tmp1+len,sizeof(char)*actual_block_size);
		tmp1[len+ret]='\0';
		right_broundry1 =len+ret;
		is_eof =feof(file);		
		// printf("ret:%d\teof:%d\n",ret,is_eof);
		// printf("%s\n",tmp1);
		if(is_eof){
			is_end1 =1;
		}else{
			offset1+=actual_block_size;			
		}




		pos1 =len;
		last=0;
		while(tmp1[pos1]!='\n'){
			pos1++;
		}
		len =pos1-last;
		pos1++;
		__builtin_prefetch (tmp1+pos1, 0);

	}
	buffer =(tmp1+last);
	return len;	
}

size_t getline_block_buffer_2_pair(FILE *file  ,char *& buffer){

	int last =pos2;
	size_t len =0;
	size_t ret =0;
	int is_eof;
	//
	while(tmp2[pos2]!='\n' && pos2<right_broundry2){
		pos2++;
	}
	len =pos2-last;


	if(tmp2[pos2]=='\n'){
		pos2++;
		buffer =(tmp2+last);
		return len;
	}else{
		if(is_end2){
			return -1;
		}
		strncpy(tmp2, (tmp2+last), len);
		tmp2[len]='\0';



		ret=fread(tmp2+len, sizeof(char),actual_block_size,  file);
		// ret=gzread(file,tmp2+len,sizeof(char)*actual_block_size);
		tmp2[len+ret]='\0';
		right_broundry2 =len+ret;
		is_eof =feof(file);

		// printf("ret:%d\teof:%d\n",ret,is_eof);
		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end2 =1;
		}else{
			offset1+=actual_block_size;			
		}


		pos2 =len;
		last=0;
		while(tmp2[pos2]!='\n'){
			pos2++;
		}
		len =pos2-last;
		pos2++;
				__builtin_prefetch (tmp2+pos2, 0);

	}

	buffer =(tmp2+last);

	return len;	

}


size_t gzgetline_block_buffer_2(gzFile file  ,char *& buffer){

	int last =pos1;
	size_t len =0;
	size_t ret =0;
	int is_eof;

	while(tmp1[pos1]!='\n' && pos1<right_broundry1){
		pos1++;
	}
	len =pos1-last;


	if(tmp1[pos1]=='\n'){
		pos1++;
		buffer =(tmp1+last);
		return len;
	}else{
		if(is_end1){
			return -1;
		}
		strncpy(tmp1, (tmp1+last), len);
		tmp1[len]='\0';




		ret=gzread(file,tmp1+len,sizeof(char)*actual_block_size);
		tmp1[len+ret]='\0';
		right_broundry1 =len+ret;
		is_eof =gzeof(file);		
		// printf("ret:%d\teof:%d\n",ret,is_eof);

		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end1 =1;
		}else{
			offset1+=actual_block_size;			
		}



		pos1 =len;
		last=0;
		while(tmp1[pos1]!='\n'){
			pos1++;
		}
		len =pos1-last;
		pos1++;
		__builtin_prefetch (tmp1+pos1, 0);

	}
	buffer =(tmp1+last);
	return len;	
}


size_t gzgetline_block_buffer_2_pair(gzFile file  ,char *& buffer){


	int last =pos2;
	size_t len =0;
	size_t ret =0;
	int is_eof;

	while(tmp2[pos2]!='\n' && pos2<right_broundry2){
		pos2++;
	}
	len =pos2-last;


	if(tmp2[pos2]=='\n'){
		pos2++;
		buffer =(tmp2+last);
		return len;
	}else{
		if(is_end2){
			return -1;
		}
		strncpy(tmp2, (tmp2+last), len);
		tmp2[len]='\0';

		ret=gzread(file,tmp2+len,sizeof(char)*actual_block_size);
		tmp2[len+ret]='\0';
		right_broundry2 =len+ret;
		is_eof =gzeof(file);

		// printf("ret:%d\teof:%d\n",ret,is_eof);
		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end2 =1;
		}else{
			offset2+=actual_block_size;			
		}


		pos2 =len;
		last=0;
		while(tmp2[pos2]!='\n'){
			pos2++;
		}
		len =pos2-last;
		pos2++;
		__builtin_prefetch (tmp2+pos2, 0);

	}
	buffer =(tmp2+last);
	return len;	

}


size_t gzgetline_block_buffer_qual_2(gzFile file  ,char *& buffer,int rlen){

	int last =pos1;
	size_t len =0;
	size_t ret =0;
	size_t thesize =0;
	int is_eof;

	if(last+rlen <right_broundry1){
		pos1 =last+rlen;
		len =rlen;
		buffer =(tmp1+last);
		pos1++;
		return len;
	}else{
		if(is_end1){
			return -1;
		}
		len =right_broundry1-last;
		strncpy(tmp1, (tmp1+last), len);	
		tmp1[len]='\0';

		ret=gzread(file,tmp1+len,sizeof(char)*actual_block_size);
		tmp1[len+ret]='\0';
		right_broundry1 =len+ret;
		is_eof =gzeof(file);		
		// printf("ret:%d\teof:%d\n",ret,is_eof);

		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end1 =1;
		}else{
			offset1+=actual_block_size;			
		}



		buffer =(tmp1);
		last=0;
		pos1=rlen;
		len =rlen;
		pos1++;
		return len;	
	}

}


size_t gzgetline_block_buffer_qual_2_pair(gzFile file  ,char *& buffer,int rlen){


	int last =pos2;
	size_t len =0;
	size_t ret =0;
	size_t thesize =0;
	int is_eof ;



	if(last+rlen <right_broundry2){
		pos2 =last+rlen;
		len =rlen;
		buffer =(tmp2+last);
		pos2++;
		return len;
	}else{
		if(is_end2){
			return -1;
		}
		len =right_broundry2-last;
		strncpy(tmp2, (tmp2+last), len);	
		tmp2[len]='\0';


		ret=gzread(file,tmp2+len,sizeof(char)*actual_block_size);
		tmp2[len+ret]='\0';
		right_broundry2 =len+ret;
		is_eof =gzeof(file);

		// printf("ret:%d\teof:%d\n",ret,is_eof);
		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end2 =1;
		}else{
			offset2+=actual_block_size;			
		}


		buffer =(tmp2);
		last=0;
		pos2=rlen;
		len =rlen;
		pos2++;
		return len;	
	}

}

size_t getline_block_buffer_qual(FILE *file ,size_t file_size ,size_t& offset,char *& buffer,char *tmp, int &pos ,int & right_broundary,int & is_end,int rlen){

	int last =pos;
	size_t len =0;
	size_t ret =0;
	size_t thesize =0;



	if(last+rlen <right_broundary){
		pos =last+rlen;
		len =rlen;
		buffer =(tmp+last);
		pos++;
		return len;
	}else{
		if(is_end){
			return -1;
		}
		len =right_broundary-last;
		strncpy(tmp, (tmp+last), len);	tmp[len]='\0';
		if(offset + actual_block_size >=file_size){
			ret=fread(tmp+len,sizeof(char), actual_block_size, file);
			tmp[len+(file_size - offset)]='\0';
			right_broundary =len+(file_size - offset);
			
			is_end =1;
		}else{
			ret=fread(tmp+len, sizeof(char),actual_block_size,  file);
			tmp[len+actual_block_size]='\0';
			right_broundary =len+actual_block_size;
			offset+=actual_block_size;
			// printf("tmp:%s\n",tmp);


		}
		buffer =(tmp);
		last=0;
		pos=rlen;
		len =rlen;
		pos++;
		return len;	
	}

}


size_t getline_block_buffer_qual_2(FILE *file  ,char *& buffer,int rlen){

	int last =pos1;
	size_t len =0;
	size_t ret =0;
	size_t thesize =0;
	int is_eof;


	if(last+rlen <right_broundry1){
		pos1 =last+rlen;
		len =rlen;
		buffer =(tmp1+last);
		pos1++;
		return len;
	}else{
		if(is_end1){
			if(last+rlen ==right_broundry1){
				pos1 =last+rlen;
				len =rlen;
				buffer =(tmp1+last);
				pos1++;
				return len;			
			}
			return -1;
		}
		len =right_broundry1-last;
		strncpy(tmp1, (tmp1+last), len);	
		tmp1[len]='\0';

		ret=fread(tmp1+len,sizeof(char), actual_block_size, file);
		// ret=gzread(file,tmp1+len,sizeof(char)*actual_block_size);
		tmp1[len+ret]='\0';
		right_broundry1 =len+ret;
		is_eof =feof(file);		
		// printf("ret:%d\teof:%d\n",ret,is_eof);

		if(is_eof){
			is_end1 =1;
		}else{
			offset1+=actual_block_size;			
		}


		buffer =(tmp1);
		last=0;
		pos1=rlen;
		len =rlen;
		pos1++;
		return len;	
	}

}

size_t getline_block_buffer_qual_2_pair(FILE *file  ,char *& buffer,int rlen){

	int last =pos2;
	size_t len =0;
	size_t ret =0;
	size_t thesize =0;
	int is_eof;


	if(last+rlen <right_broundry2){
		pos2 =last+rlen;
		len =rlen;
		buffer =(tmp2+last);
		pos2++;
		return len;
	}else{
		if(is_end2){
			if(last+rlen ==right_broundry2){
				pos2 =last+rlen;
				len =rlen;
				buffer =(tmp2+last);
				pos2++;
				return len;	
			}
			return -1;
		}
		len =right_broundry2-last;
		strncpy(tmp2, (tmp2+last), len);	
		tmp2[len]='\0';
		ret=fread(tmp2+len, sizeof(char),actual_block_size,  file);
		// ret=gzread(file,tmp2+len,sizeof(char)*actual_block_size);
		tmp2[len+ret]='\0';
		right_broundry2 =len+ret;
		is_eof =feof(file);

		// printf("ret:%d\teof:%d\n",ret,is_eof);
		if(is_eof){
			// printf("ret:%d\teof:%d\n",ret,gzeof(file));
			is_end2 =1;
		}else{
			offset2+=actual_block_size;			
		}
		buffer =(tmp2);
		last=0;
		pos2=rlen;
		len =rlen;
		pos2++;
		return len;	
	}

}


int GetNextEntry_Recycle_block(FILE *file,ReadItem_t *read,size_t file_size ,size_t& offset,char *tmp, int &pos,int & right_broundary,int &is_end){
	int p1, p2;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;

	read->header = read->seq = read->qual = NULL; read->rlen = 0;


	// if ((len = getline_block_buffer(file,file_size,offset,buffer,tmp,pos,right_broundary,is_end)) != -1)	{	
	if ((len =  getline(&buffer, &size, file)) != -1)	{	

		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		
		read->header = new char[len + 1]; 
		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, (buffer + 1), len); read.header[len] = '0';
		if (FastQFormat){
			// if ((read->rlen = getline_block_buffer(file,file_size,offset,buffer,tmp,pos,right_broundary,is_end)) != -1){

			if ((read->rlen  = getline(&buffer, &size, file)) != -1){


				read->seq = new char[read->rlen+1];
				strncpy(read->seq, buffer, read->rlen);
				getline(&buffer, &size, file);
				getline(&buffer, &size, file);
				// getline_block_buffer(file,file_size,offset,buffer,tmp,pos,right_broundary,is_end); 
				// getline_block_buffer_qual(file,file_size,offset,buffer,tmp,pos,right_broundary,is_end,read->rlen);
				read->qual = new char[read->rlen+1]; 

				
				strncpy(read->qual, buffer, read->rlen);
				// read->rlen -= 1; 
				read->seq[read->rlen] = '\0'; read->qual[read->rlen] = '\0';
				// printf("read-seq:%s\n",read->seq);

			}
			else read->rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read->rlen = (int)seq.length()) > 0)
			{
				read->seq = new char[read->rlen + 1];
				strcpy(read->seq, (char*)seq.c_str());
				read->seq[read->rlen] = '\0';
			}
		}
	}

	return read->rlen;
}

int GetNextEntry_Recycle_block_2(FILE *file,ReadItem_t *read){
	int p1, p2;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;

	// read->header = read->seq = read->qual = NULL; 
	read->rlen = 0;


	if ((len = getline_block_buffer_2(file,buffer)) != -1)	{	

		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		
		// read->header = new char[len + 1]; 
		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, (buffer + 1), len); read.header[len] = '0';
		if (FastQFormat){
			if ((read->rlen = getline_block_buffer_2(file,buffer)) != -1){
				// read->seq = new char[read->rlen+1];
				strncpy(read->seq, buffer, read->rlen);
				// printf("wdnmd 0 %d %d is_end1:%d %s\n",pos1,right_broundry1,is_end1,buffer);
				// int wdnmd1 =
				getline_block_buffer_2(file,buffer); 
				// printf("wdnmd1 %d %d %d is_end1:%d %s\n",wdnmd1,pos1,right_broundry1,is_end1,buffer);
				// int wdnmd2 =
				getline_block_buffer_qual_2(file,buffer,read->rlen);
				// read->qual = new char[read->rlen+1];
				// printf("wdnmd2 %d %d %d is_end1:%d %s\n",wdnmd2,pos1,right_broundry1,is_end1,buffer); 
				strncpy(read->qual, buffer, read->rlen);
				// read->rlen -= 1; 
				read->seq[read->rlen] = '\0'; 
				read->qual[read->rlen] = '\0';
				// printf("read-header1:\t%s\n",read->header);
				// printf("read-seq1:\t%s\n",read->seq);
				// printf("read-qual1:\t%s\n",read->qual);

				// fprintf(rfile_text,"%s\n",read->header);
				// fprintf(rfile_text,"%s\n",read->seq);
				// fprintf(rfile_text,"%s\n",read->qual);
			}else read->rlen = 0;
		}else{
			string seq;
			while (true){
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>'){
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read->rlen = (int)seq.length()) > 0){
				read->seq = new char[read->rlen + 1];
				strcpy(read->seq, (char*)seq.c_str());
				read->seq[read->rlen] = '\0';
			}
		}
	}

	return read->rlen;
}


int GetNextEntry_Recycle_block_2_pair(FILE *file,ReadItem_t *read){
	int p1, p2;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;

	// read->header = read->seq = read->qual = NULL; 
	read->rlen = 0;


	if ((len = getline_block_buffer_2_pair(file,buffer)) != -1)	{	
	// if ((len =  getline(&buffer, &size, file)) != -1)	{	

		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		
		// read->header = new char[len + 1]; 
		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, (buffer + 1), len); read.header[len] = '0';
		if (FastQFormat){
			if ((read->rlen = getline_block_buffer_2_pair(file,buffer)) != -1){

			// if ((read->rlen  = getline(&buffer, &size, file)) != -1){

				// printf("%d %d is_end2:%d :%s \n",pos2,right_broundry2, is_end2, buffer);
				// read->seq = new char[read->rlen+1];
				strncpy(read->seq, buffer, read->rlen);
				// getline(&buffer, &size, file);
				// getline(&buffer, &size, file);
				// int wdnmd1 =
				getline_block_buffer_2_pair(file,buffer); 
				// printf("%d %d %d is_end2:%d %s\n",wdnmd1,pos2,right_broundry2,is_end2,buffer);

				// int wdnmd2 =
				getline_block_buffer_qual_2_pair(file,buffer,read->rlen);
				// read->qual = new char[read->rlen+1]; 
				// printf("%d %d %d is_end2:%d %s\n",wdnmd2,pos2,right_broundry2,is_end2,buffer); 

				
				strncpy(read->qual, buffer, read->rlen);
				// read->rlen -= 1; 
				read->seq[read->rlen] = '\0'; read->qual[read->rlen] = '\0';
				// printf("read-seq:%s\n",read->seq);
				// printf("read-header2:\t%s\n",read->header);
				// printf("read-seq2:\t%s\n",read->seq);
				// printf("read-qual2:\t%s\n",read->qual);
				// fprintf(rfile_text,"%s\n",read->header);
				// fprintf(rfile_text,"%s\n",read->seq);
				// fprintf(rfile_text,"%s\n",read->qual);



			}
			else read->rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read->rlen = (int)seq.length()) > 0)
			{
				read->seq = new char[read->rlen + 1];
				strcpy(read->seq, (char*)seq.c_str());
				read->seq[read->rlen] = '\0';
			}
		}
	}

	return read->rlen;
}

int GetNextEntry_Recycle(FILE *file,ReadItem_t *read){
	int p1, p2;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;

	read->header = read->seq = read->qual = NULL; read->rlen = 0;

	if ((len = getline(&buffer, &size, file)) != -1)	{	
		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
		read->header = new char[len + 1]; 
		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';
		//len = IdentifyHeaderBoundary(buffer, len) - 1; read.header = new char[len + 1];
		//strncpy(read.header, (buffer + 1), len); read.header[len] = '0';

		if (FastQFormat)
		{
			if ((read->rlen = getline(&buffer, &size, file)) != -1)
			{
				read->seq = new char[read->rlen];
				strncpy(read->seq, buffer, read->rlen);
				getline(&buffer, &size, file); getline(&buffer, &size, file);
				read->qual = new char[read->rlen]; strncpy(read->qual, buffer, read->rlen);
				read->rlen -= 1; read->seq[read->rlen] = '\0'; read->qual[read->rlen] = '\0';
			}
			else read->rlen = 0;
		}
		else
		{
			string seq;
			while (true)
			{
				if ((len = getline(&buffer, &size, file)) == -1) break;
				if (buffer[0] == '>')
				{
					fseek(file, 0 - len, SEEK_CUR);
					break;
				}
				else
				{
					buffer[len - 1] = '\0'; seq += buffer;
				}
			}
			if ((read->rlen = (int)seq.length()) > 0)
			{
				read->seq = new char[read->rlen + 1];
				strcpy(read->seq, (char*)seq.c_str());
				read->seq[read->rlen] = '\0';
			}
		}
	}
	free(buffer);
	return read->rlen;
}




int GetNextChunk(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2)
{
	char* rseq;
	int iCount = 0;
	while (true)
	{
		if ((chunk[iCount] = GetNextEntry(file)).rlen == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (bSepLibrary) chunk[iCount] = GetNextEntry(file2);
		else chunk[iCount] = GetNextEntry(file);

		if (chunk[iCount].rlen == 0) break;
		if (bPairEnd)
		{
			rseq = new char[chunk[iCount].rlen];
			GetComplementarySeq(chunk[iCount].rlen, chunk[iCount].seq, rseq);
			copy(rseq, rseq + chunk[iCount].rlen, chunk[iCount].seq); delete[] rseq;
			if (FastQFormat)
			{
				string rqual = chunk[iCount].qual; reverse(rqual.begin(), rqual.end());
				copy(rqual.c_str(), rqual.c_str() + chunk[iCount].rlen, chunk[iCount].qual);
			}
		}
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;
}


int GetNextChunk_Recycle(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2)
{
	char* rseq;
	int iCount = 0;
	while (true){
		if ((GetNextEntry_Recycle(file,&chunk[iCount])) == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (bSepLibrary) GetNextEntry_Recycle(file2,&chunk[iCount]);
		else GetNextEntry_Recycle(file,&chunk[iCount]);

		if (chunk[iCount].rlen == 0) break;
		if (bPairEnd)
		{
			rseq = new char[chunk[iCount].rlen];
			GetComplementarySeq(chunk[iCount].rlen, chunk[iCount].seq, rseq);
			copy(rseq, rseq + chunk[iCount].rlen, chunk[iCount].seq); delete[] rseq;
			if (FastQFormat)
			{
				string rqual = chunk[iCount].qual; reverse(rqual.begin(), rqual.end());
				copy(rqual.c_str(), rqual.c_str() + chunk[iCount].rlen, chunk[iCount].qual);
			}
		}
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;
}



// int GetLines(char *tmp  ,int& incomplete  ,int lineSize ,int &lineNum ,ReadItem_t *read,char* resident_content){
	
// 	int i=0, pos =0;
// 	int len =strlen(tmp);
// 	int is_full=0;

// 	int readNo =lineNum/4;
// 	int readInfoNo =lineNum%4;
// 	int str_len =0;
// 	int right_str_len=0;
// 	char * tmp_line=nullptr;


// 	if(incomplete){
// 		while(tmp[i] !='\n' ){
// 			i++;
// 		}

// 		readNo =lineNum/4;
// 		readInfoNo =lineNum%4;

// 		tmp_line=(read+readNo)->info[readInfoNo];

// 		right_str_len =i-pos;
// 		str_len =strlen(tmp_line) +right_str_len;
// 		tmp_line =(char *)realloc(tmp_line,str_len+1);
// 		// char * wdnmd = (char *)malloc(str_len+1);
// 		// strncat(tmp_line, tmp_line, strlen(tmp_line));
// 		strncat(tmp_line, (tmp+pos), right_str_len);
// 		// free(tmp_line) ;
// 		// tmp_line =wdnmd;
// 		tmp_line[str_len] ='\0';
// 		// printf("%d %s \n",str_len,tmp_line);
		

// 		lineNum++;		
// 		incomplete =0;		
// 		if(lineNum==lineSize){

// 			memcpy(resident_content, tmp+(++i), actual_block_size-i);
// 			resident_content[actual_block_size-i]='\0';
// 			// offset_full +=(++i);
// 			is_full =1;
// 			return is_full;
// 		}


// 		pos =++i;
// 	}


// 	for( ;i<len;i++){

// 		if(tmp[i] !='\n'){
// 			continue;
// 		}else{
// 						readNo =lineNum/4;
// 			readInfoNo =lineNum%4;


// 			str_len =i-pos;	
// 			((read+readNo)->info)[readInfoNo] =(char *)malloc(str_len+1);
// 			tmp_line =((read+readNo)->info)[readInfoNo];	
// 			strncpy(tmp_line, (tmp+pos), str_len);
// 			tmp_line[str_len] ='\0';


// 			// printf("readNo:%d readInfoNo:%d %s  str_len:%d\n",readNo,readInfoNo,tmp_line,str_len);

// 			lineNum++;			
// 			if(lineNum==lineSize){
// 				memcpy(resident_content, tmp+(++i), actual_block_size-i);
// 				resident_content[actual_block_size-i]='\0';

// 				// offset_full += (++i);
// 				is_full =1;
// 				return is_full;
// 			}
// 			pos =i+1;
// 		}
// 	}

// 	//最后一个字符不是\n 要么最后一行是断开的 ，要么刚好下一个块的开头就是\n
// 	if(tmp[len-1] !='\n' && !is_full ){
// 		readNo =lineNum/4;
// 		readInfoNo =lineNum%4;

// 		str_len =i-pos;
// 		((read+readNo)->info)[readInfoNo] =(char *)malloc(str_len+1);
// 		tmp_line =((read+readNo)->info)[readInfoNo];	
// 		strncpy(tmp_line, (tmp+pos), str_len);
// 		tmp_line[str_len] ='\0';


// 		incomplete =1;
// 	}

// 	return is_full;	
// }

// void  GetLines_resident(char *tmp  ,int& incomplete  ,int &lineNum ,ReadItem_t *read){
	
// 	int i=0, pos =0;
// 	int len =strlen(tmp);

// 	int readNo =lineNum/4;
// 	int readInfoNo =lineNum%4;
// 	int str_len =0;

// 	for( ;i<len;i++){
// 		if(tmp[i] =='\n'){
// 			readNo =lineNum/4;
// 			readInfoNo =lineNum%4;
// 			// if(readInfoNo==2){
// 			// 	lineNum++;			

// 			// 	pos =i+1;
// 			// 	continue;
// 			// }
// 			str_len =i-pos;	
// 			((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 			strncpy(((read+readNo)->info)[readInfoNo], (tmp+pos), str_len);
// 			((read+readNo)->info)[readInfoNo][str_len] ='\0';

// 			lineNum++;			

// 			pos =i+1;
// 		}
// 	}

// 	//最后一个字符不是\n 要么最后一行是断开的 ，要么刚好下一个块的开头就是\n
// 	if(tmp[len-1] !='\n'  ){
// 		readNo =lineNum/4;
// 		readInfoNo =lineNum%4;
// 		str_len =i-pos;


// 		((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 		strncpy(((read+readNo)->info)[readInfoNo], (tmp+pos), str_len);
// 		((read+readNo)->info)[readInfoNo][str_len] ='\0';
// 		incomplete =1;
// 	}

// }


// int GetLines_Paired_1(char *tmp  ,int& incomplete ,size_t &offset_full ,int lineSize ,int &lineNum ,ReadItem_t *read){
	
// 	int i=0, pos =0;
// 	int len =strlen(tmp);
// 	int is_full=0;

// 	int readNo =(lineNum/4)*2;
// 	int readInfoNo =lineNum%4;
// 	int str_len =0;
// 	int right_str_len=0;
// 	char * tmp_line=nullptr;


// 	if(incomplete){
// 		while(tmp[i] !='\n' ){
// 			i++;
// 		}

// 		readNo =(lineNum/4)*2;
// 		readInfoNo =lineNum%4;
// 		tmp_line=(read+readNo)->info[readInfoNo];


// 		right_str_len =i-pos;
// 		str_len =strlen(tmp_line) +right_str_len;
// 				// printf("%s %d %d %d\n",tmp_line ,strlen(tmp_line),right_str_len,str_len);

// 		// (read+readNo)->info[readInfoNo] =(char *)realloc(tmp_line,str_len+1);
		
// 		// strncat(tmp_line, (tmp+pos), right_str_len);

// 		char * wdnmd =new char[str_len+1];
// 		strncat(wdnmd, tmp_line, strlen(tmp_line));
// 		strncat(tmp_line, (tmp+pos), right_str_len);
// 		delete tmp_line;
// 		tmp_line =wdnmd;
// 		tmp_line[str_len] ='\0';



// 		lineNum++;		
// 		incomplete =0;		
// 		if(lineNum==lineSize){
// 			offset_full +=(++i);
// 			is_full =1;
// 			return is_full;
// 		}
// 		pos =++i;
// 	}


// 	for( ;i<len;i++){
// 		if(tmp[i] =='\n'){
// 			readNo =(lineNum/4)*2;
// 			readInfoNo =lineNum%4;

// 			str_len =i-pos;	
// 			((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 			tmp_line =((read+readNo)->info)[readInfoNo];	
// 			strncpy(tmp_line, (tmp+pos), str_len);
// 			tmp_line[str_len] ='\0';


// 			// printf("readNo:%d readInfoNo:%d %s  str_len:%d\n",readNo,readInfoNo,tmp_line,str_len);

// 			lineNum++;			
// 			if(lineNum==lineSize){
// 				offset_full += (++i);
// 				is_full =1;
// 				return is_full;
// 			}
// 			pos =i+1;
// 		}
// 	}

// 	//最后一个字符不是\n 要么最后一行是断开的 ，要么刚好下一个块的开头就是\n
// 	if(tmp[len-1] !='\n' && !is_full ){
// 		readNo =(lineNum/4)*2;
// 		readInfoNo =lineNum%4;
// 		str_len =i-pos;


// 		((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 		tmp_line =((read+readNo)->info)[readInfoNo];	
// 		strncpy(tmp_line, (tmp+pos), str_len);
// 		tmp_line[str_len] ='\0';
// 				// printf("%s %d %d\n",tmp_line ,strlen(tmp_line),str_len);

// 		incomplete =1;
// 	}

// 	return is_full;	
// }

// int GetLines_Paired_2(char *tmp  ,int& incomplete ,size_t &offset_full ,int lineSize ,int &lineNum ,ReadItem_t *read){
	
// 	int i=0, pos =0;
// 	int len =strlen(tmp);
// 	int is_full=0;

// 	int readNo =(lineNum/4)*2+1;
// 	int readInfoNo =lineNum%4;
// 	int str_len =0;
// 	int right_str_len=0;
// 	char * tmp_line=nullptr;


// 	if(incomplete){
// 		while(tmp[i] !='\n' ){
// 			i++;
// 		}

// 		readNo =(lineNum/4)*2+1;
// 		readInfoNo =lineNum%4;
// 		tmp_line=(read+readNo)->info[readInfoNo];


// 		right_str_len =i-pos;
// 		str_len =strlen(tmp_line) +right_str_len;
// 		// printf("%s %d %d\n",tmp_line ,strlen(tmp_line),right_str_len);
// 		char * wdnmd =new char[str_len+1];
// 		strncat(wdnmd, tmp_line, strlen(tmp_line));
// 		strncat(tmp_line, (tmp+pos), right_str_len);
// 		delete tmp_line;
// 		tmp_line =wdnmd;
// 		tmp_line[str_len] ='\0';
		

// 		lineNum++;		
// 		incomplete =0;		
// 		if(lineNum==lineSize){
// 			offset_full +=(++i);
// 			is_full =1;
// 			return is_full;
// 		}
// 		pos =++i;
// 	}


// 	for( ;i<len;i++){
// 		if(tmp[i] =='\n'){
// 			readNo =(lineNum/4)*2+1;
// 			readInfoNo =lineNum%4;

// 			str_len =i-pos;	
// 			((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 			tmp_line =((read+readNo)->info)[readInfoNo];	
// 			strncpy(tmp_line, (tmp+pos), str_len);
// 			tmp_line[str_len] ='\0';


// 			// printf("readNo:%d readInfoNo:%d %s  str_len:%d\n",readNo,readInfoNo,tmp_line,str_len);

// 			lineNum++;			
// 			if(lineNum==lineSize){
// 				offset_full += (++i);
// 				is_full =1;
// 				return is_full;
// 			}
// 			pos =i+1;
// 		}
// 	}

// 	//最后一个字符不是\n 要么最后一行是断开的 ，要么刚好下一个块的开头就是\n
// 	if(tmp[len-1] !='\n' && !is_full ){
// 		readNo =(lineNum/4)*2+1;
// 		readInfoNo =lineNum%4;
// 		str_len =i-pos;


// 		((read+readNo)->info)[readInfoNo] =new char[str_len+1];
// 		tmp_line =((read+readNo)->info)[readInfoNo];	
// 		strncpy(tmp_line, (tmp+pos), str_len);
// 		tmp_line[str_len] ='\0';

// 		incomplete =1;
// 	}

// 	return is_full;	
// }


// int GetNextBlock(size_t &offset_full,int chunkSize,size_t file_size,size_t & offset,FILE *file,ReadItem_t *read ,int &readNum ,char* resident_content){

// 	// char *tmp = nullptr;
// 	int incomplete =0;
// 	int end =0;
// 	int lineSize =chunkSize*4;
// 	size_t ret;
// 	char *tmp =(char *)malloc(block_size);
// 	tmp[block_size - 1] ='\0';

// 	int lineNum =0;


// 	GetLines_resident( resident_content , incomplete ,lineNum ,read);



// 	// while (!address_array.try_dequeue(tmp));
// 	while(1){
// 		// printf("%d\n",lineNum);
// 		if(offset + actual_block_size >=file_size){
// 			memset(tmp, 0, actual_block_size+1);
// 			ret = fread(tmp, file_size - offset, 1, file);
// 			// memcpy(tmp, mapped , file_size - offset);
// 			GetLines( tmp , incomplete ,lineSize,lineNum ,read,resident_content);
// 			end =1;
// 			break;
// 		}else{
// 			ret = fread(tmp, actual_block_size, 1, file);
// 			// memcpy(tmp, mapped, actual_block_size);
// 			//读取直到获取到对应chunkSize
// 			if(!GetLines( tmp , incomplete ,lineSize,lineNum,read, resident_content)){
// 				offset+=actual_block_size;
// 			}else{
// 				offset+=actual_block_size;
// 				// offset+=offset_full;
// 				// fseek(file, -(actual_block_size-offset_full), SEEK_CUR);
// 				break;
// 			}
// 		}
// 	}
// 	// printf("lineNums:%d chunkSize:%d\n",lineNum,chunkSize);

// 	readNum =lineNum/4;

// 	// for(int i =0;i<readNum;i++){
// 	// 	printf("%s\n",(read+i)->info[0]);
// 	// }

// 	// address_array.enqueue(tmp);
// 	free(tmp);
// 	return end;

// }



// int GetNextBlock_Paired_1(size_t &offset_full,int chunkSize,size_t file_size,size_t & offset,FILE *file,ReadItem_t *read ,int &readNum){

// 	// char *tmp = nullptr;
// 	int incomplete =0;
// 	int end =0;
// 	int lineSize =chunkSize*4;
// 	size_t ret;
// 	char *tmp =(char *)malloc(block_size);
// 	tmp[block_size - 1] ='\0';

// 	int lineNum =0;

// 	//pair=0则是不用取反向互补
// 	//pair=1则是需要取反向互补

// 	// while (!address_array.try_dequeue(tmp));
// 	while(1){
// 		if(offset + actual_block_size >=file_size){
// 			memset(tmp, 0, actual_block_size+1);
// 			ret = fread(tmp, file_size - offset, 1, file);
// 			// memcpy(tmp, mapped , file_size - offset);
// 			GetLines_Paired_1( tmp , incomplete ,offset_full,lineSize,lineNum ,read);
// 			end =1;
// 			break;
// 		}else{
// 			ret = fread(tmp, actual_block_size, 1, file);
// 			// memcpy(tmp, mapped, actual_block_size);
// 			//读取直到获取到对应chunkSize
// 			if(!GetLines_Paired_1( tmp , incomplete ,offset_full,lineSize,lineNum,read)){
// 				offset+=actual_block_size;
// 			}else{
// 				offset+=offset_full;
// 				fseek(file, -(actual_block_size-offset_full), SEEK_CUR);
// 				break;
// 			}
// 		}
// 	}
// 	// printf("lineNums:%d chunkSize:%d\n",lineNum,chunkSize);

// 	readNum =lineNum/4;

// 	// for(int i =0;i<readNum;i++){
// 	// 	printf("%s\n",(read+i)->info[0]);
// 	// }

// 	// address_array.enqueue(tmp);
// 	free(tmp);
// 	return end;

// }



// int GetNextBlock_Paired_2(size_t &offset_full,int chunkSize,size_t file_size,size_t & offset,FILE *file,ReadItem_t *read ,int &readNum){

// 	// char *tmp = nullptr;
// 	int incomplete =0;
// 	int end =0;
// 	int lineSize =chunkSize*4;
// 	size_t ret;
// 	char *tmp =(char *)malloc(block_size);
// 	tmp[block_size - 1] ='\0';

// 	int lineNum =0;

// 	//pair=0则是不用取反向互补
// 	//pair=1则是需要取反向互补

// 	// while (!address_array.try_dequeue(tmp));
// 	while(1){
// 		if(offset + actual_block_size >=file_size){
// 			memset(tmp, 0, actual_block_size+1);
// 			ret = fread(tmp, file_size - offset, 1, file);
// 			// memcpy(tmp, mapped , file_size - offset);
// 			GetLines_Paired_2( tmp , incomplete ,offset_full,lineSize,lineNum ,read);
// 			end =1;
// 			break;
// 		}else{
// 			ret = fread(tmp, actual_block_size, 1, file);
// 			// memcpy(tmp, mapped, actual_block_size);
// 			//读取直到获取到对应chunkSize
// 			if(!GetLines_Paired_2( tmp , incomplete ,offset_full,lineSize,lineNum,read)){
// 				offset+=actual_block_size;
// 			}else{
// 				offset+=offset_full;
// 				fseek(file, -(actual_block_size-offset_full), SEEK_CUR);
// 				break;
// 			}
// 		}
// 	}
// 	// printf("lineNums:%d chunkSize:%d\n",lineNum,chunkSize);

// 	readNum =lineNum/4;

// 	// for(int i =0;i<readNum;i++){
// 	// 	printf("%s\n",(read+i)->info[0]);
// 	// }

// 	// address_array.enqueue(tmp);
// 	free(tmp);
// 	return end;

// }

void GetRead(vector<string> &res1,ReadItem_t *tmp_Read,int &i,int is_paired){

	int p1, p2,len,rlen;
	len =res1[i].size();
	p1 = IdentifyHeaderBegPos(res1[i], len); 
	p2 = IdentifyHeaderEndPos(res1[i], len); len = p2 - p1;
	tmp_Read->header = new char[len + 1]; 	
	strncpy(tmp_Read->header, (res1[i].c_str() + p1), len); tmp_Read->header[len] = '\0';i++;


	if (is_paired){
		len =res1[i].size();
		tmp_Read->rlen =len;
		tmp_Read->seq = new char[len+1];
		GetComplementarySeq1(len, res1[i].c_str(), tmp_Read->seq);
		tmp_Read->seq[len] ='\0';i+=2;

		len =res1[i].size();
		tmp_Read->qual = new char[len+1];
		// if (FastQFormat){
		// 	reverse(res1[i].begin(), res1[i].end());
		// }			
		reverse(res1[i].begin(), res1[i].end());
		strncpy(tmp_Read->qual, res1[i].c_str(), len);tmp_Read->qual[len] ='\0';i++;

	}else{
		len =res1[i].size();
		tmp_Read->rlen =len;
		tmp_Read->seq = new char[len+1];
		strncpy(tmp_Read->seq, res1[i].c_str(), len);tmp_Read->seq[len] ='\0';i+=2;

		len =res1[i].size();
		tmp_Read->qual = new char[len+1];
		strncpy(tmp_Read->qual, res1[i].c_str(), len);tmp_Read->qual[len] ='\0';i++;
	}

}

// int GetNextEntry_Block(FILE *file,size_t & offset,ReadItem_t *read ,int chunkSize ,size_t file_size ,int &readNum,char * resident_content){
// 	// char *tmp = nullptr;

// 	size_t offset_full =0;
// 	int incomplete =0;
// 	int end =0;


// 	int lineSize =chunkSize*4;
// 	size_t ret;
// 	char *tmp =(char *)malloc(block_size);
// 	tmp[block_size - 1] ='\0';

// 	int lineNum =0;

// 	if(strlen(resident_content) >0){
// 		GetLines_resident( resident_content , incomplete ,lineNum ,read);
// 	}

// 	// printf("lineNum:%d incomplete:%d len:%d\n",lineNum,incomplete ,strlen(resident_content));

// 	// while (!address_array.try_dequeue(tmp));
// 	while(1){
// 		// printf("%d\n",lineNum);
// 		if(offset + actual_block_size >=file_size){
// 			memset(tmp, 0, actual_block_size+1);
// 			ret = fread(tmp, file_size - offset, 1, file);
// 			// memcpy(tmp, mapped , file_size - offset);
// 			GetLines( tmp , incomplete ,lineSize,lineNum ,read,resident_content);
// 			end =1;
// 			break;
// 		}else{
// 			ret = fread(tmp, actual_block_size, 1, file);
// 			// memcpy(tmp, mapped, actual_block_size);
// 			//读取直到获取到对应chunkSize
// 			if(!GetLines( tmp , incomplete ,lineSize,lineNum,read, resident_content)){
// 				offset+=actual_block_size;
// 			}else{
// 				offset+=actual_block_size;
// 				break;
// 			}
// 		}
// 	}
// 	// printf("lineNums:%d chunkSize:%d\n",lineNum,chunkSize);

// 	readNum =lineNum/4;

// 	// for(int i =0;i<readNum;i++){
// 	// 	printf("%s\n",(read+i)->info[0]);
// 	// }

// 	// address_array.enqueue(tmp);
// 	free(tmp);



// 	return end;
// }




// int GetNextEntry_Block_Pair(int &readNum1,int &readNum2,int chunkSize,ReadItem_t *read,FILE *file,size_t & offset1, size_t file_size1, 
// 																		FILE *file2,size_t & offset2,size_t file_size2 ){
// 	char *tmp = nullptr;

// 	size_t offset_full1 =0,offset_full2=0;
// 	int incomplete =0;
// 	int i=0,j=0,end1,end2;
// 	int lineNum =0;



// 	end1 =GetNextBlock_Paired_1( offset_full1, chunkSize/2, file_size1, offset1,file,read,readNum1);
// 	end2 =GetNextBlock_Paired_2( offset_full2, chunkSize/2, file_size2, offset2,file2,read,readNum2);


// 	// for(i=0;i<chunkSize;i++){
// 	// 	int p1, p2,len;
// 	// 	char* tmp_header =(read+i)->info[0];
// 	// 	len =strlen(tmp_header);
// 	// 	p1 = IdentifyHeaderBegPos(tmp_header, len); 
// 	// 	p2 = IdentifyHeaderEndPos(tmp_header, len); len = p2 - p1;
// 	// 	memcpy(tmp_header, tmp_header+p1, len);
// 	// 	tmp_header[len]='\n';
// 	// 	// tmp_Read->header = new char[len + 1]; 	
// 	// 	// strncpy(tmp_Read->header, (res1[i].c_str() + p1), len); tmp_Read->header[len] = '\0';i++;
// 	// }

// 	// for(int x=0;x<res2.size();x++){
// 	// 	printf("fastq2[%d]:%s\n",x,res2[x].c_str());
// 	// }
// 	// while(i<lineNum){
// 	// 	// printf("header1[%d]:%s\n",i,res1[i].c_str());
// 	// 	// printf("header2[%d]:%s\n",j,res2[j].c_str());
// 	// 	ReadItem_t *tmp_Read1  =(read++);
// 	// 	ReadItem_t *tmp_Read2 =(read++);
// 	// 	GetRead(res1, tmp_Read1,i,0);readNum++;
// 	// 	GetRead(res2, tmp_Read2,j,bPairEnd);readNum++;
// 	// }

// 	return end1;
// }


//  int GetNextChunk_Block(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2,size_t& offset1,size_t& offset2,int &is_end,char * resident_content){

// 		int readNum=0,i;
// 	int readNum1=0,readNum2=0,tmp_len =0;

// 	// printf("bSepLibrary:%d bSepLibrary:%d\n",bSepLibrary ,FastQFormat);
// 	//获取一批reads,数量是chunkSize
// 	if (bSepLibrary){

// 		is_end =GetNextEntry_Block_Pair(readNum1,readNum2,chunkSize,chunk,file,offset1,fastq_file_size1,file2,offset2,fastq_file_size2);
// 		readNum =readNum1+readNum2;
// 		for(i=1;i<readNum;i+=2 ){
// 			char* tmp_seq =(chunk+i)->info[1];
// 			tmp_len =strlen(tmp_seq);
// 			GetComplementarySeq2(tmp_len,tmp_seq);
// 			char* tmp_qual =(chunk+i)->info[3];
// 			tmp_len =strlen(tmp_qual);
// 			GetReverseSeq(tmp_len,tmp_qual);		
// 		}

// 	}else{
// 		is_end =GetNextEntry_Block(file,offset1,chunk ,chunkSize,fastq_file_size1,readNum, resident_content);
		
// 	}

// 	for( i=0;i<readNum;i++){
// 		int p1, p2,len,headerLen;
// 		char* tmp_header =(chunk+i)->info[0];
// 		len =strlen(tmp_header);
// 		p1 = IdentifyHeaderBegPos(tmp_header, len); 
// 		p2 = IdentifyHeaderEndPos(tmp_header, len); headerLen = p2 - p1;
// 		memcpy(tmp_header, tmp_header+p1, headerLen);
// 		tmp_header[headerLen]='\0';
// 		// memset(tmp_header+headerLen,0,len-headerLen);
// 		// printf("i:%d %s\n",i,tmp_header);
// 	}

// 	// ReadItem_t *temp =chunk;
// 	// printf("%d\n",readNum);
// 	// for( i=0;i<readNum;i++){
// 	// 	printf("header:%s\n",(chunk+i)->info[0]);
// 	// 	printf("seq:%s\n",(chunk+i)->info[1]);
// 	// 	printf("qual:%s\n",(chunk+i)->info[3]);
// 	// }

// 	return readNum;
//  }



int GetNextChunk_Block_1(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2,
						size_t& offset1,size_t& offset2,char * tmp1,char *tmp2,
						int &pos1,int &pos2,int & right_broundary1,int &right_broundary2,int &is_end1,int &is_end2){

	char* rseq;
	int iCount = 0;
	// printf("file_size:%ld pos:%d offset:%d\n",fastq_file_size1,pos1,offset1);

	while (true){
		if ((GetNextEntry_Recycle_block(file,&chunk[iCount],fastq_file_size1,offset1,tmp1,pos1,right_broundary1,is_end1)) == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (bSepLibrary) GetNextEntry_Recycle_block(file2,&chunk[iCount],fastq_file_size2,offset2,tmp2,pos2,right_broundary2,is_end2);
		else GetNextEntry_Recycle_block(file,&chunk[iCount],fastq_file_size1,offset1,tmp1,pos1,right_broundary1,is_end1);

		if (chunk[iCount].rlen == 0) break;
		if (bPairEnd){

			GetComplementarySeq2(chunk[iCount].rlen,chunk[iCount].seq);

			// rseq = new char[chunk[iCount].rlen];
			// GetComplementarySeq(chunk[iCount].rlen, chunk[iCount].seq, rseq);
			// copy(rseq, rseq + chunk[iCount].rlen, chunk[iCount].seq); delete[] rseq;
			if (FastQFormat)
			{	
							GetReverseSeq(chunk[iCount].rlen,chunk[iCount].qual);		

				// string rqual = chunk[iCount].qual; reverse(rqual.begin(), rqual.end());
				// copy(rqual.c_str(), rqual.c_str() + chunk[iCount].rlen, chunk[iCount].qual);
			}
		}
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;

}
int GetNextChunk_Block_2(ReadItem_t *chunk,int bSepLibrary,int chunkSize, int flag,FILE *file, FILE *file2){

	char* rseq;
	int iCount = 0;
	// printf("file_size:%ld pos:%d offset:%d\n",fastq_file_size1,pos1,offset1);

	while (true){
		if ((GetNextEntry_Recycle_block_2(file,&chunk[iCount])) == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		//if (bPacBioData) bp_sum += ReadArr[iCount].rlen;
		iCount++;
		if (bSepLibrary) GetNextEntry_Recycle_block_2_pair(file2,&chunk[iCount]);
		else GetNextEntry_Recycle_block_2(file,&chunk[iCount]);

		if (chunk[iCount].rlen == 0) break;
		if (bPairEnd){
			GetComplementarySeq2(chunk[iCount].rlen,chunk[iCount].seq);
			if (FastQFormat){	
				GetReverseSeq(chunk[iCount].rlen,chunk[iCount].qual);		
			}
		}

		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;

}

int gzGetNextEntry_block(gzFile file, ReadItem_t *read){

	int p1, p2;
	size_t len;
	size_t size = 0;

	
	// ReadItem_t read;
	char *buffer = NULL;

	// read->header = read->seq = read->qual = NULL; 
	read->rlen = 0;



	if ((len = gzgetline_block_buffer_2(file,buffer)) != -1){
		len = strlen(buffer); 
		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len);
		len = p2 - p1;

		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';


		if (FastQFormat){
			if ((read->rlen = gzgetline_block_buffer_2(file,buffer)) != -1){
				// read->seq = new char[read->rlen+1];
				strncpy(read->seq, buffer, read->rlen);
				gzgetline_block_buffer_2(file,buffer); 
				gzgetline_block_buffer_qual_2(file,buffer,read->rlen);
				// read->qual = new char[read->rlen+1]; 
				strncpy(read->qual, buffer, read->rlen);
				// read->rlen -= 1; 
				read->seq[read->rlen] = '\0'; read->qual[read->rlen] = '\0';
				// printf("read-header1:\t%s\n",read->header);
				// printf("read-seq1:\t%s\n",read->seq);
				// printf("read-qual1:\t%s\n",read->qual);
			}else{
				read->rlen = 0;
			} 
		}
		// }else{
		// 	// string seq;
		// 	// while (true){
		// 	// 	if ((len = getline(&buffer, &size, file)) == -1) break;
		// 	// 	if (buffer[0] == '>'){
		// 	// 		fseek(file, 0 - len, SEEK_CUR);
		// 	// 		break;
		// 	// 	}
		// 	// 	else{
		// 	// 		buffer[len - 1] = '\0'; seq += buffer;
		// 	// 	}
		// 	// }
		// 	// if ((read->rlen = (int)seq.length()) > 0){
		// 	// 	read->seq = new char[read->rlen + 1];
		// 	// 	strcpy(read->seq, (char*)seq.c_str());
		// 	// 	read->seq[read->rlen] = '\0';
		// 	// }
		// }

	}


	return read->rlen;




}


int gzGetNextEntry_block_pair(gzFile file, ReadItem_t *read)
{
	int p1, p2;
	size_t len;
	size_t size = 0;

	// ReadItem_t read;
	char *buffer = NULL;

	// read->header = read->seq = read->qual = NULL; 
	read->rlen = 0;

	if ((len = gzgetline_block_buffer_2_pair(file,buffer)) != -1)
	{
		len = strlen(buffer); 
		p1 = IdentifyHeaderBegPos(buffer, len); 
		p2 = IdentifyHeaderEndPos(buffer, len);
		len = p2 - p1;

		strncpy(read->header, (buffer + p1), len); 
		read->header[len] = '\0';


		if (FastQFormat){
			if ((read->rlen = gzgetline_block_buffer_2_pair(file,buffer)) != -1){
				// read->seq = new char[read->rlen+1];
				strncpy(read->seq, buffer, read->rlen);
				gzgetline_block_buffer_2_pair(file,buffer); 
				gzgetline_block_buffer_qual_2_pair(file,buffer,read->rlen);
				// read->qual = new char[read->rlen+1]; 
				strncpy(read->qual, buffer, read->rlen);
				// read->rlen -= 1; 
				read->seq[read->rlen] = '\0'; 
				read->qual[read->rlen] = '\0';
				// printf("read-header2:\t%s\n",read->header);
				// printf("read-seq2:\t%s\n",read->seq);
				// printf("read-qual2:\t%s\n",read->qual);

			}else read->rlen = 0;
		}


	}


	return read->rlen;

}

ReadItem_t gzGetNextEntry(gzFile file)
{
	int p1, p2;
	char* buffer;
	ReadItem_t read;
	int len, buf_size;
	
	if (bPacBioData) buf_size = 1000000;
	else buf_size = 1000;

	buffer = new char[buf_size];

	read.header = read.seq = read.qual = NULL; read.rlen = 0;

	if (gzgets(file, buffer, buf_size) != NULL)
	{
		len = strlen(buffer); p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len);
		len = p2 - p1;
		if (len > 0 && (buffer[0] == '@' || buffer[0] == '>'))
		{
			//p1 = IdentifyHeaderBegPos(buffer, len); p2 = IdentifyHeaderEndPos(buffer, len); len = p2 - p1;
			//read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
			read.header = new char[len + 1]; strncpy(read.header, (buffer + p1), len); read.header[len] = '\0';
			gzgets(file, buffer, buf_size); read.rlen = strlen(buffer) - 1; read.seq = new char[read.rlen + 1]; read.seq[read.rlen] = '\0';
			strncpy(read.seq, buffer, read.rlen);

			if (FastQFormat)
			{
				gzgets(file, buffer, buf_size); gzgets(file, buffer, buf_size);
				read.qual = new char[read.rlen + 1]; read.qual[read.rlen] = '\0';
				strncpy(read.qual, buffer, read.rlen);
			}
		}
	}
	delete[] buffer;

	return read;
}


int gzGetNextChunk_block(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* chunk,int chunkSize, int flag)
{
	char* rseq;
	int iCount = 0;

	while (true)
	{
		if (( gzGetNextEntry_block(file,&chunk[iCount])) == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iCount++;
		if (bSepLibrary){
			gzGetNextEntry_block_pair(file2,&chunk[iCount]);
		}  
		else {
			gzGetNextEntry_block(file,&chunk[iCount]);
		} 

		if (chunk[iCount].rlen == 0) break;

		if (bPairEnd)
		{
			GetComplementarySeq2(chunk[iCount].rlen,chunk[iCount].seq);
			if (FastQFormat){	
				GetReverseSeq(chunk[iCount].rlen,chunk[iCount].qual);		
			}
		}
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;

}



int gzGetNextChunk(bool bSepLibrary, gzFile file, gzFile file2, ReadItem_t* ReadArr,int chunkSize, int flag)
{
	char* rseq;
	int iCount = 0;

	while (true)
	{
		if ((ReadArr[iCount] = gzGetNextEntry(file)).rlen == 0) break;
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iCount++;

		if (bSepLibrary) ReadArr[iCount] = gzGetNextEntry(file2);
		else ReadArr[iCount] = gzGetNextEntry(file);

		if (ReadArr[iCount].rlen == 0) break;

		if (bPairEnd)
		{
			rseq = new char[ReadArr[iCount].rlen];
			GetComplementarySeq(ReadArr[iCount].rlen, ReadArr[iCount].seq, rseq);
			copy(rseq, rseq + ReadArr[iCount].rlen, ReadArr[iCount].seq);
			delete[] rseq;
			if (FastQFormat)
			{
				string rqual = ReadArr[iCount].qual; reverse(rqual.begin(), rqual.end());
				copy(rqual.c_str(), rqual.c_str() + ReadArr[iCount].rlen, ReadArr[iCount].qual);
			}
		}
		//ReadArr[iCount].EncodeSeq = new uint8_t[ReadArr[iCount].rlen];
		//for (i = 0; i != ReadArr[iCount].rlen; i++) ReadArr[iCount].EncodeSeq[i] = nst_nt4_table[(int)ReadArr[iCount].seq[i]];
		iCount++;
		if (iCount == ReadChunkSize || (bPacBioData && iCount == 10)) break;
	}
	return iCount;
}



static inline int __occ_aux(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = ((c&2)? y : ~y) >> 1 & ((c&1)? y : ~y) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

bwtint_t bwt_occ(const bwt_t *bwt, bwtint_t k, ubyte_t c)
{
	bwtint_t n;
	uint32_t *p, *end;

	if (k == bwt->seq_len) return bwt->L2[c+1] - bwt->L2[c];
	if (k == (bwtint_t)(-1)) return 0;
	k -= (k >= bwt->primary); // because $ is not in bwt

	// retrieve Occ at k/OCC_INTERVAL
	n = ((bwtint_t*)(p = bwt_occ_intv(bwt, k)))[c];
	p += sizeof(bwtint_t); // jump to the start of the first BWT cell

	// calculate Occ up to the last k/32
	end = p + (((k>>5) - ((k&~OCC_INTV_MASK)>>5))<<1);
	for (; p < end; p += 2) n += __occ_aux((uint64_t)p[0]<<32 | p[1], c);

	// calculate Occ
	n += __occ_aux(((uint64_t)p[0]<<32 | p[1]) & ~((1ull<<((~k&31)<<1)) - 1), c);
	if (c == 0) n -= ~k&31; // corrected for the masked bits

	return n;
}

void bwt_occ4(const bwt_t *bwt, bwtint_t k, bwtint_t cnt[4])
{
	bwtint_t x;
	uint32_t *p, tmp, *end;
	if (k == (bwtint_t)(-1)) {
		memset(cnt, 0, 4 * sizeof(bwtint_t));
		return;
	}
	k -= (k >= bwt->primary); // because $ is not in bwt
	p = bwt_occ_intv(bwt, k);
	memcpy(cnt, p, 4 * sizeof(bwtint_t));
	p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
	end = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4)); // this is the end point of the following loop
	for (x = 0; p < end; ++p) x += __occ_aux4(bwt, *p);
	tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
	x += __occ_aux4(bwt, tmp) - (~k&15);
	cnt[0] += x&0xff; cnt[1] += x>>8&0xff; cnt[2] += x>>16&0xff; cnt[3] += x>>24;
}

void bwt_2occ4(const bwt_t *bwt, bwtint_t k, bwtint_t l, bwtint_t cntk[4], bwtint_t cntl[4])
{
	bwtint_t _k, _l;
	_k = k - (k >= bwt->primary);
	_l = l - (l >= bwt->primary);
	if (_l >> OCC_INTV_SHIFT != _k >> OCC_INTV_SHIFT || k == (bwtint_t)(-1) || l == (bwtint_t)(-1)) {
		bwt_occ4(bwt, k, cntk);
		bwt_occ4(bwt, l, cntl);
	} else {
		bwtint_t x, y;
		uint32_t *p, tmp, *endk, *endl;
		k -= (k >= bwt->primary); // because $ is not in bwt
		l -= (l >= bwt->primary);
		p = bwt_occ_intv(bwt, k);
		memcpy(cntk, p, 4 * sizeof(bwtint_t));
		p += sizeof(bwtint_t); // sizeof(bwtint_t) = 4*(sizeof(bwtint_t)/sizeof(uint32_t))
		// prepare cntk[]
		endk = p + ((k>>4) - ((k&~OCC_INTV_MASK)>>4));
		endl = p + ((l>>4) - ((l&~OCC_INTV_MASK)>>4));
		for (x = 0; p < endk; ++p) x += __occ_aux4(bwt, *p);
		y = x;
		tmp = *p & ~((1U<<((~k&15)<<1)) - 1);
		x += __occ_aux4(bwt, tmp) - (~k&15);
		// calculate cntl[] and finalize cntk[]
		for (; p < endl; ++p) y += __occ_aux4(bwt, *p);
		tmp = *p & ~((1U<<((~l&15)<<1)) - 1);
		y += __occ_aux4(bwt, tmp) - (~l&15);
		memcpy(cntl, cntk, 4 * sizeof(bwtint_t));
		cntk[0] += x&0xff; cntk[1] += x>>8&0xff; cntk[2] += x>>16&0xff; cntk[3] += x>>24;
		cntl[0] += y&0xff; cntl[1] += y>>8&0xff; cntl[2] += y>>16&0xff; cntl[3] += y>>24;
	}
}

static inline bwtint_t bwt_invPsi(const bwt_t *bwt, bwtint_t k) // compute inverse CSA
{
	bwtint_t x = k - (k > bwt->primary);
	x = bwt_B0(bwt, x);
	x = bwt->L2[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary? 0 : x;
}

bwtint_t bwt_sa(bwtint_t k)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;
	while (k & mask) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	return sa + bwt->sa[k/bwt->sa_intv];
}


void bwt_sa_batch(bwtint_t k,int interval,int match_len,int match_beg,vector<SeedPair_t> &SeedPairVec)
{
	bwtint_t sa = 0, mask = bwt->sa_intv - 1;

		// printf("%d\n",bwt->sa_intv);
		uint64_t il =k;
		for(uint64_t i=0;i<interval;i++){
			sa =0;
			k =il;
			while (k & mask) {
				++sa;
				k = bwt_invPsi(bwt, k);
			}
			SeedPair_t SeedPair;
			SeedPair.bSimple = true;
			// pos_list[i] =(sa + bwt->sa[k/bwt->sa_intv]);
			SeedPair.gPos =sa + bwt->sa[k/bwt->sa_intv];
			SeedPair.rPos = match_beg;
			SeedPair.PosDiff = SeedPair.gPos -SeedPair.rPos;
			SeedPair.rLen = SeedPair.gLen = match_len;


			SeedPairVec.push_back(SeedPair);			

			il++;
		}

	/* without setting bwt->sa[0] = -1, the following line should be
	   changed to (sa + bwt->sa[k/bwt->sa_intv]) % (bwt->seq_len + 1) */
	// return sa + bwt->sa[k/bwt->sa_intv];
}


/*********************
 * Bidirectional BWT *
 *********************/

void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
{
	bwtint_t tk[4], tl[4];
	int i;
	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
	for (i = 0; i != 4; ++i) {
		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
		ok[i].x[2] = tl[i] - tk[i];
	}
	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
}

// bwtSearchResult_t bwtExactMatchBackward(const bwt_t *bwt,int query_begin, int len, const ubyte_t *str, bwtint_t *sa_begin, bwtint_t *sa_end,bwtintv_t *ik,bwtintv_t *ok)
// {

// 	bwtSearchResult_t bwtSearchResult;
// 	int i, pos, p;
// 	bwtintv_t *ik, ok[4];
// 	bwt_set_intv(bwt, str[query_begin+len - 1], *ik); // the initial interval of a single base
// 	ubyte_t c ;
// 	ik->info =1;
// 	// printf("%d ",c);
// 	// printf("k:%ld intv:%ld\n",ik->x[0],ik->x[2]);

// 	// bwtint_t k, l, ok, ol,lastK,lastL;
// 	int i;
// 	// k = 0; l = bwt->seq_len;
// 	for (i = query_begin+len - 2; i >= query_begin  ; --i) {
// 		c = str[i]; // complement of q[i]
// 		// printf("%d ",c);

// 		//这个应该是x[2]interval
// 		if (ik->x[2] < 0) { // an interval small enough
// 			break;
// 		} else if (str[i] < 4) { // an A/C/G/T base
			
// 			bwt_extend(bwt, ik, ok, 1);
// 			if (ok[c].x[2] != ik->x[2]) { // change of the interval size
// 				if (ok[c].x[2] <= 0) {
// 					break; // the interval size is too small to be extended further
// 				}
// 			}
// 			ik->x[0] = ok[c].x[0]; 
// 			ik->x[1] = ok[c].x[1]; 
// 			ik->x[2] = ok[c].x[2]; 			
// 			// ik->info = i + 1;
// 			ik->info ++;
// 		} else { // an ambiguous base
// 			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
// 		}
// 	}

// 	while(ik->x[2] >20 && i>=0){
		
// 		c = str[i];
// 		// printf("%d ",c);

// 		if (str[i] < 4) { // an A/C/G/T base
// 			bwt_extend(bwt, ik, ok, 1);
// 			if (ok[c].x[2] != ik->x[2]) { // change of the interval size
// 				if (ok[c].x[2] <= 0) {
// 					break; // the interval size is too small to be extended further
// 				}
// 			}
// 			ik->x[0] = ok[c].x[0]; 
// 			ik->x[1] = ok[c].x[1]; 
// 			ik->x[2] = ok[c].x[2]; 			
// 			ik->info ++;
// 			i--;
// 		} else { // an ambiguous base
// 			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
// 		}
// 	}
// 		// printf("k:%ld intv:%ld\n",ik->x[0],ik->x[2]);
// 	return  ik->x[2];
// 	// if (sa_begin) *sa_begin = k;
// 	// if (sa_end)   *sa_end = l;
// 	// return l - k + 1;
// }


bwtSearchResult_t BWT_Search(uint8_t* seq, int start, int stop)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;


	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}


	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len = pos - start) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{
			bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}
	return bwtSearchResult;
}


bwtSearchResult_t BWT_Only_Search(uint8_t* seq, int start, int stop)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;


	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}
	bwtSearchResult.len = pos - start;
	bwtSearchResult.freq = (int)ik.x[2];
	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	return bwtSearchResult;
}


bwtSearchResult_t BWT_Search_Forward(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];



	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;

	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}
	
	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len = pos - start) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(start);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}
		bwtSearchResult.rightset =pos;

	// if(stop-pos >=MinSeedLength){
	// 	int com_start =pos/10*10;
	// 	int com_pos ;
	// 	int match_len;
	// 	p = (int)seq[com_start];
	// 	ik.x[0] = bwt->L2[p] + 1;
	// 	ik.x[1] = bwt->L2[3 - p] + 1;
	// 	//interval
	// 	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	// 	for (com_pos = com_start + 1; com_pos < stop; com_pos++){
	// 		if (seq[com_pos] > 3) break;// ambiguous base
	// 		//正向匹配
	// 		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
	// 		for (i = 0; i != 4; ++i) {
	// 			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
	// 			ok[i].x[2] = tl[i] - tk[i];
	// 		}
	// 		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
	// 		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
	// 		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
	// 		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

	// 		i = 3 - seq[com_pos];
	// 		//interval
	// 		if (ok[i].x[2] == 0){
	// 			break; // extension ends
	// 		} else{
	// 			ik = ok[i];
	// 		} 
	// 	}
	// 	match_len =com_pos -com_start;
	// 	if(match_len >=MinSeedLength && (int)ik.x[2] <= OCC_Thr){
	// 			il_list->push_back(ik.x[0]);
	// 			interval_list->push_back((int)ik.x[2]);
	// 			match_len_list->push_back(match_len);
	// 			match_beg_list->push_back(com_start);

	// 	}
	// 	bwtSearchResult.rightset =com_start+match_len;

	// }



	return bwtSearchResult;
}


bwtSearchResult_t BWT_Search_Forward_1(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];



	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;

	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}
	
	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len = pos - start) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(start);

			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else {
			il_list->push_back(ik.x[0]);
			interval_list->push_back(OCC_Thr);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(start);
		}
	}
	bwtSearchResult.rightset =pos;

	if(stop-pos >=MinSeedLength){
		int com_start =pos/10*10;
		int com_pos ;
		int match_len;
		int rightest=0;
		p = (int)seq[com_start];
		ik.x[0] = bwt->L2[p] + 1;
		ik.x[1] = bwt->L2[3 - p] + 1;
		//interval
		ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

		for (com_pos = com_start + 1; com_pos < stop; com_pos++){
			if (seq[com_pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

			i = 3 - seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		rightest =com_pos;

		for (com_pos = com_start-1; com_pos >=0 ; com_pos--){
			if (seq[com_pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 

		}

		match_len =rightest -(com_pos+1);
		if(match_len >=MinSeedLength ){
				// fprintf(stderr,"%d-%d\n",com_pos,match_len);

				il_list->push_back(ik.x[0]);
				if( (int)ik.x[2] <= OCC_Thr){
					interval_list->push_back((int)ik.x[2]);
				}else{
					interval_list->push_back(OCC_Thr);
				}
				
				match_len_list->push_back(match_len);
				match_beg_list->push_back(com_pos);

		}
		bwtSearchResult.rightset =rightest;

	}

	return bwtSearchResult;
}




bwtSearchResult_t BWT_Search_Forward_2(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p,rightest ,leftest;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];



	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;
	bwtSearchResult.leftest=0;
	bwtSearchResult.len=0;
	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}

	rightest =pos;
	for (pos = start-1; pos >=0 ; pos--){
		if (seq[pos] > 3) break;// ambiguous base
		//反向匹配
		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
		i =  seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 

	}
	leftest =pos+1;
	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len =rightest -leftest ) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		bwtSearchResult.rightset =rightest;
		bwtSearchResult.leftest =leftest;
		fprintf(stderr,"wdnmd:start:%d %d-%d-%ld\n",start ,bwtSearchResult.leftest,bwtSearchResult.leftest+bwtSearchResult.len,ik.x[0]);

		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr){	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}else {
			il_list->push_back(ik.x[0]);
			interval_list->push_back(OCC_Thr);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
		}
	}
	if(stop-rightest >=MinSeedLength  ){
		
		int com_start =   bwtSearchResult.rightset%10==0 ?  ((bwtSearchResult.rightset-1)/10)*10   :  bwtSearchResult.rightset/10*10;
		int com_pos ;
		int com_match_len=0;
		int com_rightest=0;
		int com_leftest =0;
		p = (int)seq[com_start];
		ik.x[0] = bwt->L2[p] + 1;
		ik.x[1] = bwt->L2[3 - p] + 1;
		//interval
		ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

		for (com_pos = com_start + 1; com_pos < stop; com_pos++){
			if (seq[com_pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			i = 3 - seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_rightest =com_pos;

		for (com_pos = com_start-1; com_pos >=0 ; com_pos--){
			if (seq[com_pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_leftest =com_pos+1;
		com_match_len =com_rightest -com_leftest;
			fprintf(stderr,"nmsl:com_start:%d %d-%d-%ld\n",com_start,com_leftest,com_leftest+com_match_len,ik.x[0]);
		//
		if(com_match_len >=MinSeedLength  && (com_rightest >bwtSearchResult.rightset )){
			il_list->push_back(ik.x[0]);
			if( (int)ik.x[2] <= OCC_Thr){
				interval_list->push_back((int)ik.x[2]);
			}else{
				interval_list->push_back(OCC_Thr);
			}
			match_len_list->push_back(com_match_len);
			match_beg_list->push_back(com_leftest);
		}
		bwtSearchResult.rightset =com_rightest;
	}
	
	return bwtSearchResult;
}

bwtSearchResult_t BWT_Search_Forward_3(uint8_t* seq, int start, int stop ,int last_rightset,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p,rightest ,leftest;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];



	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;
	bwtSearchResult.leftest=0;
	bwtSearchResult.len=0;
	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}

	rightest =pos;
	for (pos = start-1; pos >=0 ; pos--){
		if (seq[pos] > 3) break;// ambiguous base
		//反向匹配
		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
		i =  seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 

	}
	leftest =pos+1;
	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len =rightest -leftest ) < MinSeedLength  || rightest <=last_rightset){
		bwtSearchResult.rightset =last_rightset;
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		bwtSearchResult.rightset =rightest;
		bwtSearchResult.leftest =leftest;
		// fprintf(stderr,"wdnmd:start:%d %d-%d-%d --%ld\n",start ,leftest,rightest,last_rightset,ik.x[0]);

		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr){	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}else {
			il_list->push_back(ik.x[0]);
			interval_list->push_back(OCC_Thr);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
		}
	}
	if(stop-rightest >=MinSeedLength){
		int com_start =   rightest%10==0 ?  ((rightest-1)/10)*10   :  rightest/10*10;
		int com_pos ;
		int com_match_len=0;
		int com_rightest=0;
		int com_leftest =0;
		p = (int)seq[com_start];
		ik.x[0] = bwt->L2[p] + 1;
		ik.x[1] = bwt->L2[3 - p] + 1;
		//interval
		ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

		for (com_pos = com_start + 1; com_pos < stop; com_pos++){
			if (seq[com_pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			i = 3 - seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_rightest =com_pos;

		for (com_pos = com_start-1; com_pos >=0 ; com_pos--){
			if (seq[com_pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_leftest =com_pos+1;
		com_match_len =com_rightest -com_leftest;
		//	fprintf(stderr,"nmsl:com_start:%d %d-%d-%d  --%ld\n",com_start,com_leftest,com_leftest+com_match_len,bwtSearchResult.rightset,ik.x[0]);

		if(com_match_len >=MinSeedLength  && (com_rightest >bwtSearchResult.rightset )){
			il_list->push_back(ik.x[0]);
			if( (int)ik.x[2] <= OCC_Thr){
				interval_list->push_back((int)ik.x[2]);
			}else{
				interval_list->push_back(OCC_Thr);
			}
			match_len_list->push_back(com_match_len);
			match_beg_list->push_back(com_leftest);
		}
		bwtSearchResult.rightset =com_rightest;
	}
	
	return bwtSearchResult;
}




bwtSearchResult_t BWT_Search_Forward_hash(uint8_t* seq, int start, int stop ,int last_rightset,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p,rightest ,leftest,temp_start;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;
	// printf("forward33333333\n");
	uint32_t kmer_no =(int)seq[start] & 0b11;
	// printf("%d",kmer_no);
	for(int i=1;i<10;i++){
			// printf("%d",(int)seq[start+i]);
		kmer_no = (kmer_no<<2) +((int)seq[start+i] & 0b11);
	}
	temp_start =start;
	start +=9;

	// printf("\n");
	// printf("kmer_no:%d\n",kmer_no);

	// printf("%ld %ld %ld\n",bwtIdx->hash->kmer_hash[kmer_no].x[0],bwtIdx->hash->kmer_hash[kmer_no].x[1],bwtIdx->hash->kmer_hash[kmer_no].x[2]);

	// p = (int)seq[start];
	// ik.x[0] = bwt->L2[p] + 1;
	// ik.x[1] = bwt->L2[3 - p] + 1;
	// //interval
	// ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];
	// 		printf("%d %ld %ld %ld\n",seq[start],ik.x[0],ik.x[1],ik.x[2]);

	ik.x[0] = bwtIdx->hash->kmer_hash[kmer_no].x[0];
	ik.x[1] = bwtIdx->hash->kmer_hash[kmer_no].x[1];
	ik.x[2] = bwtIdx->hash->kmer_hash[kmer_no].x[2];//interval

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;
	bwtSearchResult.leftest=0;
	bwtSearchResult.len=0;

	if(ik.x[2] !=0){


		
		for (pos = start + 1; pos < stop; pos++){
			if (seq[pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			i = 3 - seq[pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
				// printf("pos:%d %d %ld %ld %ld\n",pos,seq[pos],ik.x[0],ik.x[1],ik.x[2]);

			} 
		}

		rightest =pos;
		for (pos = temp_start-1; pos >=0 ; pos--){
			if (seq[pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 

		}
		leftest =pos+1;
	}


	


	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
		// printf("bwt_search: pos=%d len=%d, freq=%d\n", start, rightest -leftest, (int)ik.x[2]);
	if ((bwtSearchResult.len =rightest -leftest ) < MinSeedLength  || rightest <=last_rightset){
		bwtSearchResult.rightset =last_rightset;
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		bwtSearchResult.rightset =rightest;
		bwtSearchResult.leftest =leftest;
		// fprintf(stderr,"wdnmd:start:%d %d-%d-%d --%ld\n",start ,leftest,rightest,last_rightset,ik.x[0]);

		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr){	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}else {
			il_list->push_back(ik.x[0]);
			interval_list->push_back(OCC_Thr);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.leftest);
		}
	}


	if(stop-rightest >=MinSeedLength){
		int com_start =   rightest%10==0 ?  ((rightest-1)/10)*10   :  rightest/10*10;
		int tmp_Start =com_start;
		int com_pos ;
		int com_match_len=0;
		int com_rightest=0;
		int com_leftest =0;

		kmer_no =(int)seq[com_start] & 0b11;
		// printf("%d",kmer_no);
		for(int i=1;i<10;i++){
				// printf("%d",(int)seq[start+i]);
			kmer_no = (kmer_no<<2) +((int)seq[com_start+i] & 0b11);
		}
		com_start +=9;
		ik.x[0] = bwtIdx->hash->kmer_hash[kmer_no].x[0];
		ik.x[1] = bwtIdx->hash->kmer_hash[kmer_no].x[1];
		ik.x[2] = bwtIdx->hash->kmer_hash[kmer_no].x[2];//interval

		// p = (int)seq[com_start];
		// ik.x[0] = bwt->L2[p] + 1;
		// ik.x[1] = bwt->L2[3 - p] + 1;
		// //interval
		// ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

		for (com_pos = com_start + 1; com_pos < stop; com_pos++){
			if (seq[com_pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			i = 3 - seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_rightest =com_pos;

		for (com_pos = tmp_Start-1; com_pos >=0 ; com_pos--){
			if (seq[com_pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		com_leftest =com_pos+1;
		com_match_len =com_rightest -com_leftest;
		//	fprintf(stderr,"nmsl:com_start:%d %d-%d-%d  --%ld\n",com_start,com_leftest,com_leftest+com_match_len,bwtSearchResult.rightset,ik.x[0]);

		if(com_match_len >=MinSeedLength  && (com_rightest >bwtSearchResult.rightset )){
			il_list->push_back(ik.x[0]);
			if( (int)ik.x[2] <= OCC_Thr){
				interval_list->push_back((int)ik.x[2]);
			}else{
				interval_list->push_back(OCC_Thr);
			}
			match_len_list->push_back(com_match_len);
			match_beg_list->push_back(com_leftest);
		}
		bwtSearchResult.rightset =com_rightest;
	}
	
	return bwtSearchResult;
}






bwtSearchResult_t BWT_Search_BothSide(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;
	bwtSearchResult.len=0;
	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		bwtSearchResult.len ++;

		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}

	bwtSearchResult.rightset =pos;
	for (pos = start-1; pos >=0 ; pos--){
		if (seq[pos] > 3) break;// ambiguous base
		//反向匹配
		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
		bwtSearchResult.len ++;
		i =  seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 

	}

	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len ) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.rightset-bwtSearchResult.len);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}

	return bwtSearchResult;
}

bwtSearchResult_t BWT_Search_BothSide1(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;

	p = (int)seq[start];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.rightset=0;
	bwtSearchResult.len=0;
	for (pos = start + 1; pos < stop; pos++){
		if (seq[pos] > 3) break;// ambiguous base
		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		bwtSearchResult.len ++;

		i = 3 - seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}
	
	bwtSearchResult.rightset =pos;
	for (pos = start-1; pos >=0 ; pos--){
		if (seq[pos] > 3) break;// ambiguous base
		//反向匹配
		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
		bwtSearchResult.len ++;
		i =  seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 

	}

	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
	if ((bwtSearchResult.len ) < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.rightset-bwtSearchResult.len);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}else{
			il_list->push_back(ik.x[0]);
			interval_list->push_back(OCC_Thr);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(bwtSearchResult.rightset-bwtSearchResult.len);
		} 
	}


	if(stop-bwtSearchResult.rightset >=MinSeedLength){
		int com_start =bwtSearchResult.rightset/10*10;
		int com_pos ;
		int match_len;
		int rightest=0;
		p = (int)seq[com_start];
		ik.x[0] = bwt->L2[p] + 1;
		ik.x[1] = bwt->L2[3 - p] + 1;
		//interval
		ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

		for (com_pos = com_start + 1; com_pos < stop; com_pos++){
			if (seq[com_pos] > 3) break;// ambiguous base
			//正向匹配
			bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[1] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];

			i = 3 - seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		rightest =com_pos;

		for (com_pos = com_start-1; com_pos >=0 ; com_pos--){
			if (seq[com_pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			i =  seq[com_pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 

		}

		match_len =rightest -com_pos;
		if(match_len >=MinSeedLength ){
				il_list->push_back(ik.x[0]);
				if( (int)ik.x[2] <= OCC_Thr){
				interval_list->push_back((int)ik.x[2]);					
				}else{
					interval_list->push_back(OCC_Thr);
				}

				match_len_list->push_back(match_len);
				match_beg_list->push_back(com_pos);

		}
		bwtSearchResult.rightset =rightest;

	}















	return bwtSearchResult;
}

// bwtSearchResult_t BWT_Search_Backward(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
// {
// 	int i, pos, p;
// 	bwtintv_t ik, ok[4];
// 	bwtint_t tk[4], tl[4];
// 	bwtSearchResult_t bwtSearchResult;
// 	int right =start;

// 	p = (int)seq[right-1];
// 	ik.x[0] = bwt->L2[p] + 1;
// 	ik.x[1] = bwt->L2[3 - p] + 1;
// 	//interval
// 	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

// 	bwtSearchResult.freq = 0; 
// 	bwtSearchResult.LocArr = NULL;
	

// 	for (pos = right-2; pos >=0 ; pos--){
// 		if (seq[pos] > 3) break;// ambiguous base
// 		//反向匹配
// 		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
// 		for (i = 0; i != 4; ++i) {
// 			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
// 			ok[i].x[2] = tl[i] - tk[i];
// 		}
// 		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
// 		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
// 		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
// 		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];

// 		i =  seq[pos];
// 		//interval
// 		if (ok[i].x[2] == 0){
// 			break; // extension ends
// 		} else{
// 			ik = ok[i];
// 		} 
// 	}
	
// 	//if (bDebugMode) printf("bwt_search: pos=%d len=%d, freq=%d\n", start, pos - start, (int)ik.x[2]);
// 	if ((bwtSearchResult.len = start - pos) < MinSeedLength){
// 		bwtSearchResult.freq = 0;
// 	} else{
// 		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
// 		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
// 		{	
// 			il_list->push_back(ik.x[0]);
// 			interval_list->push_back(bwtSearchResult.freq);
// 			match_len_list->push_back(bwtSearchResult.len);
// 			match_beg_list->push_back(start);
// 			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
// 			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
// 		}
// 		else bwtSearchResult.freq = 0;
// 	}
// 	return bwtSearchResult;
// }

bwtSearchResult_t BWT_Search_Backward_hash(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;
	// int right =start;

	// p = (int)seq[right-1];
	// ik.x[0] = bwt->L2[p] + 1;
	// ik.x[1] = bwt->L2[3 - p] + 1;
	// //interval
	// ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.len =1;	

	int right =start-10;
	uint32_t kmer_no =(int)seq[right] & 0b11;		// printf("%d",kmer_no);
	for(int i=1;i<10;i++){
		// printf("%d",(int)seq[start+i]);
		kmer_no = (kmer_no<<2) +((int)seq[right+i] & 0b11);
	}

	ik.x[0] = bwtIdx->hash->kmer_hash[kmer_no].x[0];
	ik.x[1] = bwtIdx->hash->kmer_hash[kmer_no].x[1];
	ik.x[2] = bwtIdx->hash->kmer_hash[kmer_no].x[2];//interval
	if(ik.x[2] >0){
		bwtSearchResult.freq = ik.x[2]; 
		bwtSearchResult.LocArr = NULL;
		bwtSearchResult.len =10;	


		for (pos = right-1; pos >=0 ; pos--){
			if (seq[pos] > 3) break;// ambiguous base
			//反向匹配
			bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
			for (i = 0; i != 4; ++i) {
				ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
				ok[i].x[2] = tl[i] - tk[i];
			}
			ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
			ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
			ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
			ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
			bwtSearchResult.len ++;
			i =  seq[pos];
			//interval
			if (ok[i].x[2] == 0){
				break; // extension ends
			} else{
				ik = ok[i];
			} 
		}
		pos = pos <0 ? 0: pos;
		bwtSearchResult.il =(ik.x[0]);
		bwtSearchResult.freq = (int)ik.x[2];
	}

	
	if (bwtSearchResult.len  < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(start-bwtSearchResult.len);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}

	return bwtSearchResult;
}




bwtSearchResult_t BWT_Search_Backward(uint8_t* seq, int start, int stop ,vector<bwtint_t>* il_list, vector<int>* interval_list,vector<int>* match_len_list,vector<int>* match_beg_list)
{
	int i, pos, p;
	bwtintv_t ik, ok[4];
	bwtint_t tk[4], tl[4];
	bwtSearchResult_t bwtSearchResult;
	int right =start;

	p = (int)seq[right-1];
	ik.x[0] = bwt->L2[p] + 1;
	ik.x[1] = bwt->L2[3 - p] + 1;
	//interval
	ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];

	bwtSearchResult.freq = 0; 
	bwtSearchResult.LocArr = NULL;
	bwtSearchResult.len =1;	

	for (pos = right-2; pos >=0 ; pos--){
		if (seq[pos] > 3) break;// ambiguous base
		//反向匹配
		bwt_2occ4(bwt, ik.x[0] - 1, ik.x[0] - 1 + ik.x[2], tk, tl);
		for (i = 0; i != 4; ++i) {
			ok[i].x[0] = bwt->L2[i] + 1 + tk[i];
			ok[i].x[2] = tl[i] - tk[i];
		}
		ok[3].x[1] = ik.x[1] + (ik.x[0] <= bwt->primary && ik.x[0] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[1] = ok[3].x[1] + ok[3].x[2];
		ok[1].x[1] = ok[2].x[1] + ok[2].x[2];
		ok[0].x[1] = ok[1].x[1] + ok[1].x[2];
		bwtSearchResult.len ++;
		i =  seq[pos];
		//interval
		if (ok[i].x[2] == 0){
			break; // extension ends
		} else{
			ik = ok[i];
		} 
	}
	pos = pos <0 ? 0: pos;
	bwtSearchResult.il =(ik.x[0]);
	bwtSearchResult.freq = (int)ik.x[2];
	
	if (bwtSearchResult.len  < MinSeedLength){
		bwtSearchResult.freq = 0;
	} else{
		//if ((bwtSearchResult.freq = (int)ik.x[2]) <= (bPacBioData ? 1000 : OCC_Thr))
		if ((bwtSearchResult.freq = (int)ik.x[2]) <= OCC_Thr)
		{	
			il_list->push_back(ik.x[0]);
			interval_list->push_back(bwtSearchResult.freq);
			match_len_list->push_back(bwtSearchResult.len);
			match_beg_list->push_back(start-bwtSearchResult.len);
			// bwtSearchResult.LocArr = new bwtint_t[bwtSearchResult.freq];
			// for (i = 0; i < bwtSearchResult.freq; i++) bwtSearchResult.LocArr[i] = bwt_sa(ik.x[0] + i);
		}
		else bwtSearchResult.freq = 0;
	}

	return bwtSearchResult;
}





// void bwt_extend(const bwt_t *bwt, const bwtintv_t *ik, bwtintv_t ok[4], int is_back)
// {
// 	bwtint_t tk[4], tl[4];
// 	int i;
// 	bwt_2occ4(bwt, ik->x[!is_back] - 1, ik->x[!is_back] - 1 + ik->x[2], tk, tl);
// 	for (i = 0; i != 4; ++i) {
// 		ok[i].x[!is_back] = bwt->L2[i] + 1 + tk[i];
// 		ok[i].x[2] = tl[i] - tk[i];
// 	}
// 	ok[3].x[is_back] = ik->x[is_back] + (ik->x[!is_back] <= bwt->primary && ik->x[!is_back] + ik->x[2] - 1 >= bwt->primary);
// 	ok[2].x[is_back] = ok[3].x[is_back] + ok[3].x[2];
// 	ok[1].x[is_back] = ok[2].x[is_back] + ok[2].x[2];
// 	ok[0].x[is_back] = ok[1].x[is_back] + ok[1].x[2];
// }


bwtSearchResult_t bwtExactMatchForward(uint8_t* str, int start, int rlen)
{
	bwtSearchResult_t bwtSearchResult;
	int i, pos, p;
	bwtintv_t *ik, ok[4];
	// bwtint_t tk[4], tl[4];
	//如果使用正向搜索，第一个则需要是原版序列的第一个
	//这里是建立第一个碱基的interval
	bwt_set_intv(bwt,str[start], *ik); // the initial interval of a single base
	ik->info =1;
	ubyte_t c ;

	// int i;
	// k = 0; l = bwt->seq_len;
	for (i = start+1; i < rlen   ; ++i) {
		c =3 - str[i]; // complement of q[i]
		//这个应该是x[2]interval
		if (ik->x[2] < 0) { // an interval small enough
			break;
		} else if (str[i] < 4) { // an A/C/G/T base
			bwt_extend(bwt, ik, ok, 0);
			if (ok[c].x[2] != ik->x[2]) { // change of the interval size
				if (ok[c].x[2] <= 0) {
					break; // the interval size is too small to be extended further
				}
			}
			ik->x[0] = ok[c].x[0]; 
			ik->x[1] = ok[c].x[1]; 
			ik->x[2] = ok[c].x[2]; 			
			ik->info ++;
		} else { // an ambiguous base
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	while(ik->x[2] >20 && i<rlen){
		c =3 - str[i];
		if (str[i] < 4) { // an A/C/G/T base
			bwt_extend(bwt, ik, ok, 0);
			if (ok[c].x[2] != ik->x[2]) { // change of the interval size
				if (ok[c].x[2] <= 0) {
					// printf("nmsllll:%d\n",len);
					break; // the interval size is too small to be extended further
				}
			}
			ik->x[0] = ok[c].x[0]; 
			ik->x[1] = ok[c].x[1]; 
			ik->x[2] = ok[c].x[2]; 			
			ik->info ++;
			i++;
		} else { // an ambiguous base
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	return  bwtSearchResult;
	// if (sa_begin) *sa_begin = k;
	// if (sa_end)   *sa_end = l;
	// return l - k + 1;
}



