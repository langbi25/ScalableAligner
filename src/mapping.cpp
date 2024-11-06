#include "util.h"
#include "htslib/htslib/sam.h"
#include "htslib/htslib/kseq.h"
#include "htslib/htslib/kstring.h"

#define MAPQ_COEF 30
#define Max_MAPQ  60

// float QUERY_THRESHOLD = 0.9;
// HashPair* hashes;
// UncompressedBF* bf;
FILE *bffile;

BloomTree* root;

time_t StartProcessTime;
FILE *sam_out = 0;
thread_info input_info, output_info;
samFile *bam_out = 0;
bool bSepLibrary = false;
FILE *ReadFileHandler1, *ReadFileHandler2;

gzFile gzReadFileHandler1, gzReadFileHandler2;

int64_t iDistance = 0;
int64_t iTotalReadNum = 0, iUniqueMapping = 0, iUnMapping = 0, iPaired = 0;

bam_hdr_t *header = NULL;
moodycamel::ConcurrentQueue<QueueItem_t*> spsc_queue_input;
moodycamel::ConcurrentQueue<SamOutputItem_t*> spsc_queue_output;
moodycamel::ConcurrentQueue<SamOutputItem_t*> spsc_queue_output_addr;


int ReadChunkSize,InputQueueSize;


std::atomic<int> process_num{0};



std::mutex run_mutex_t;
std::condition_variable run_cv_t;


size_t fastq_file_size1=0 ,fastq_file_size2=0;

/*初试化内存块*/
size_t blocks = 4096 ;
// size_t block_size = 1048576*4; // 4K一个块
size_t block_size = 4096; // 4K一个块

size_t actual_block_size = block_size - 1;

size_t output_block_size = 1048576; // 8M一个块


size_t block_offset =0 ,output_len =0,fwrite_size =1048576*16;
	// char *block_buff;
// FILE *rfile_text = fopen("test_content.txt", "w");

int sam_out_fd;	



struct aiocb rd1,rd2;


size_t offset1=0 ,offset2 =0; 
int pos1=0,pos2=0;
char *tmp1 ;
char *tmp2 ;
char *tmp1_block_1;
char *tmp1_block_2;

size_t aio_pointer_1=0;
size_t aio_pointer_2=0;

int right_broundry1=0,right_broundry2=0;

int is_end1=0,is_end2=0;


priority_queue<SamOutputItem_t*, vector<SamOutputItem_t*>, cmp_index > priQueMinFirst;

void init_resources(int n, ...) {
    va_list arg_ptr ;
    int i;
    va_start(arg_ptr, n);
    thread_info * tmp_info = NULL;

    for(i = 0; i < n; i++) {
        tmp_info = va_arg(arg_ptr, thread_info *);
        pthread_mutex_init(&(tmp_info->ioLock), NULL);
        pthread_mutex_init(&(tmp_info->runLock), NULL);

        pthread_cond_init(&(tmp_info->ioCond), NULL);
        pthread_cond_init(&(tmp_info->runCond), NULL);


       
    }

    va_end(arg_ptr);
}

void free_resources(int n, ...) {
    va_list arg_ptr;
    int i;
    va_start(arg_ptr, n);
    thread_info * tmp_info = NULL;
    for(i = 0; i < n; i++) {
        tmp_info = va_arg(arg_ptr, thread_info *);
        pthread_mutex_destroy(&(tmp_info->ioLock));
        pthread_mutex_destroy(&(tmp_info->runLock));
        pthread_cond_destroy(&(tmp_info->ioCond));
        pthread_cond_destroy(&(tmp_info->runCond));


    }
    va_end(arg_ptr);
}



int CountLines(string filename)
{
    ifstream ReadFile;
    int n=0;
    string tmp;
    ReadFile.open(filename,ios::in);//ios::in 表示以只读的方式读取文件
    if(ReadFile.fail())//文件打开失败:返回0
    {
        return 0;
    }
    else//文件存在
    {
        while(getline(ReadFile,tmp,'\n'))
        {
            n++;
        }
        ReadFile.close();
        return n;
    }
}




// static void *process_inPut(void * args){
// 	threadInputArg * input_arg = (threadInputArg *) args;

// 	int i =0,tail,head,chunkSize;

//     while(1) {
// 		QueueItem_t *queueItem;
// 		if(!spsc_queue_input.try_dequeue(queueItem)){
// 			    fprintf(stderr, "腾不出空间给写线程\n"); 

// 			continue;
// 		}else{


// 			// queueItem->readNum = GetNextChunk(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);
// 			queueItem->readNum =GetNextChunk_Recycle(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);
// 			// queueItem->SamOutputVec =new vector<string>[queueItem->readNum];
// 			// queueItem->AlignmentVec1 =new vector<AlignmentCandidate_t>[queueItem->readNum];
// 			// queueItem->SeedPairVec1 =new vector<SeedPair_t>[queueItem->readNum];
// 			// queueItem->EncodeSeq =new uint8_t*[queueItem->readNum];

// 			if (iPaired >= 1000)
// 			{
// 				queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
// 				queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
// 			}else{
// 				queueItem->EstDistance = maxInsertSize;
// 			}  

// 			if(queueItem->readNum==0){
// 				input_info.shutdown=1;
// 				break;
// 			}else{
// 				queueItem->index =i++;
// 				// while(!address_array.enqueue(queueItem));   
// 				address_array.enqueue(queueItem) ;
// 				// fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld\n ",queueItem->index,queueItem->readNum , queueItem);

// 			}
// 		}
//     }

//     return NULL;
// }



// void process_read(Executor<QueueItem_t *> *executor){

// 	int i =0,tail,head,chunkSize;
// 	// fprintf(stderr, "读线程\n"); 

//     while(1) {
// 		QueueItem_t *queueItem;
// 		if(!spsc_queue_input.try_dequeue(queueItem)){
// 			// fprintf(stderr, "腾不出空间给读线程\n"); 
// 			continue;
// 		}else{
// 			// queueItem->readNum = GetNextChunk(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);
// 			queueItem->readNum =GetNextChunk_Recycle(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);

// 			if (iPaired >= 1000){
// 				queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
// 				queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
// 			}else{
// 				queueItem->EstDistance = maxInsertSize;
// 			}  

// 			if(queueItem->readNum==0){
// 				break;
// 			}else{
// 				queueItem->index =i++;
// 				executor->enqueue(queueItem);
// 				process_num.fetch_add(1);
// 				// fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld -%ld\n ",queueItem->index,queueItem->readNum , queueItem,queueItem->readArr);
// 			}
// 		}
//     }
// }


// void process_read(){

// 	int i =0,tail,head,chunkSize=ReadChunkSize;
// 	// fprintf(stderr, "读线程\n"); 




//     while(1) {
// 		QueueItem_t *queueItem=new QueueItem_t;
// 		queueItem->readArr =new ReadItem_t[chunkSize];
// 			// queueItem->readNum = GetNextChunk(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);
// 		queueItem->readNum =GetNextChunk_Recycle(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);
// 		// for(int x=0;x<queueItem->readNum;x++){
// 		// 					printf("%s\n",queueItem->readArr[x].header);
// 		// 		printf("%s\n",queueItem->readArr[x].qual);
// 		// 		printf("%s\n",queueItem->readArr[x].seq);
// 		// 	// if(strlen(queueItem->readArr[x].seq)!=250    ||strlen(queueItem->readArr[x].qual)!=250 ){
							
// 		// 	// 	printf("%s\n",queueItem->readArr[x].header);

// 		// 	// 	printf("%s\n",queueItem->readArr[x].seq);
// 		// 	// 	printf("%s\n",queueItem->readArr[x].qual);

// 		// 	// 	}

// 		// }
// 		if (iPaired >= 1000){
// 			queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
// 			queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
// 		}else{
// 			queueItem->EstDistance = maxInsertSize;
// 		}  

// 		if(queueItem->readNum==0){
// 			printf("读线程退出\n");
// 			break;
// 		}else{
// 			queueItem->index =i++;
// 			process_num.fetch_add(1);
// 			fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld -%ld\n ",queueItem->index,queueItem->readNum , queueItem,queueItem->readArr);
// 		}

//     }
// }



void   process_read_block_aio(Executor<QueueItem_t *> *executor){

	int i =0,tail,head,pos =0;

	int chunkSize=ReadChunkSize;



	time_t	StartReadTime = time(NULL);




	if(gzCompressed){

		int ret1 =gzread(gzReadFileHandler1,tmp1,sizeof(char)* actual_block_size);
		tmp1[actual_block_size] ='\0';
		right_broundry1=ret1;
		offset1+=ret1;
		// fprintf(stderr, "tmp1:%s\n",tmp1);


		if(bSepLibrary){
			int ret2 =gzread(gzReadFileHandler2,tmp2, sizeof(char)* actual_block_size);
			tmp2[ret2] ='\0';
			right_broundry2=ret2;
			offset2+=ret2;
		// fprintf(stderr, "tmp2:%s\n",tmp2);

		}	
	}else{
		//将rd结构体清空
		bzero(&rd1,sizeof(struct aiocb));
		//为rd.aio_buf分配空间
		rd1.aio_buf = tmp1_block_1+(block_size *((aio_pointer_1++)%2));
		//填充rd结构体
		rd1.aio_fildes = fileno(ReadFileHandler1);
		rd1.aio_nbytes =  block_size;
		rd1.aio_offset = 0;

		//将rd结构体清空
		bzero(&rd2,sizeof(struct aiocb));
		//为rd.aio_buf分配空间
		rd2.aio_buf = tmp1_block_2+(block_size *((aio_pointer_2++)%2));
		//填充rd结构体
		rd2.aio_fildes = fileno(ReadFileHandler2);
		rd2.aio_nbytes =  block_size;
		rd2.aio_offset = 0;
		


		int ret1=aio_read(&rd1);
		while ( aio_error( &rd1 ) == EINPROGRESS ) ;
		if ((ret1 = aio_return( &rd1 )) > 0) {
			
			
			strncpy(tmp1, (const char *)rd1.aio_buf, ret1);
			tmp1[ret1]='\0';
			right_broundry1 =ret1;
			
			rd1.aio_offset +=ret1;
			rd1.aio_buf = tmp1_block_1+(block_size *((aio_pointer_1++)%2));
			ret1 = aio_read(&rd1);
		} else {
			fprintf(stderr, "wrong1!!!\n"); 
		}



		if(bSepLibrary){
			int ret2 =aio_read(&rd2);
			while ( aio_error( &rd2 ) == EINPROGRESS ) ;
			if ((ret2 = aio_return( &rd2 )) > 0) {
				strncpy(tmp2, (const char *)rd2.aio_buf, ret2);
				tmp2[ret2]='\0';
				right_broundry2 =ret2;

				rd2.aio_offset +=ret2;
				rd2.aio_buf = tmp1_block_2+(block_size *((aio_pointer_2++)%2));
				ret2 = aio_read(&rd2);			
			} else {
				fprintf(stderr, "wrong2!!!\n"); 
			}
		}		
	}



    while(1) {
		QueueItem_t *queueItem;
		if(!spsc_queue_input.try_dequeue(queueItem)){
			// fprintf(stderr, "腾不出空间给读线程\n"); 
			continue;
		}else{
			if(gzCompressed){
				queueItem->readNum = gzGetNextChunk_block(bSepLibrary, gzReadFileHandler1, gzReadFileHandler2, queueItem->readArr,chunkSize,tail);

			} else{
				queueItem->readNum =GetNextChunk_Block_aio_2(queueItem->readArr,bSepLibrary ,chunkSize,tail ,ReadFileHandler1,ReadFileHandler2);
			}
			
			// queueItem->readNum = GetNextChunk(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);


			if (iPaired >= 1000){
				queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
				queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
			}else{
				queueItem->EstDistance = maxInsertSize;
			}  
			//最后一批

			if(queueItem->readNum==0){
				break;
			}else{
				queueItem->index =i++;
				executor->enqueue(queueItem);
				process_num.fetch_add(1);
				// fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld -%ld\n ",queueItem->index,queueItem->readNum , queueItem,queueItem->readArr);
			}

		}
    }
	free(tmp1);
	free(tmp2);

	free(tmp1_block_1);
	free(tmp1_block_2);

	if (gzCompressed){
		if (gzReadFileHandler1 != NULL) gzclose(gzReadFileHandler1);
		if (gzReadFileHandler2 != NULL) gzclose(gzReadFileHandler2);
	}else{
		if (ReadFileHandler1 != NULL) fclose(ReadFileHandler1);
		if (ReadFileHandler2 != NULL) fclose(ReadFileHandler2);
	}
	long long readTime =(long long)(time(NULL) - StartReadTime);
	fprintf(stderr, "\rAll the  reads have been read in %lld seconds.\taio\t%d\n",readTime,block_size);


	// precision, recall
	FILE *rfile = fopen("io_exp_2.csv", "a");

	fprintf(rfile, "单线程异步:%d\t%s: %s\t%d\t%.2lld\n",block_size,readFile1.c_str(),outputFile,threadNum, readTime);

	fclose(rfile);

	// fclose(rfile_text);
	printf("读线程结束\n");
}


void process_read(Executor<QueueItem_t *> *executor){

	int i =0,tail,head,pos =0;
	// fprintf(stderr, "读线程\n"); 
	int chunkSize=ReadChunkSize;

    while(1) {
		QueueItem_t *queueItem;
		if(!spsc_queue_input.try_dequeue(queueItem)){
			// fprintf(stderr, "腾不出空间给读线程\n"); 
			continue;
		}else{

			if(gzCompressed){
				queueItem->readNum =gzGetNextChunk(bSepLibrary,gzReadFileHandler1,gzReadFileHandler2,queueItem->readArr ,chunkSize,tail);
			} else{
				queueItem->readNum =GetNextChunk_Recycle(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);

			} 


			if (iPaired >= 1000){
				queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
				queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
			}else{
				queueItem->EstDistance = maxInsertSize;
			}  

			if(queueItem->readNum==0){
				break;
			}else{
				queueItem->index =i++;
				executor->enqueue(queueItem);
				process_num.fetch_add(1);
				// fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld -%ld\n ",queueItem->index,queueItem->readNum , queueItem,queueItem->readArr);
			}

			
		}
    }

	free(tmp1);
	free(tmp2);
	if (gzCompressed){
		if (gzReadFileHandler1 != NULL) gzclose(gzReadFileHandler1);
		if (gzReadFileHandler2 != NULL) gzclose(gzReadFileHandler2);
	}else{
		if (ReadFileHandler1 != NULL) fclose(ReadFileHandler1);
		if (ReadFileHandler2 != NULL) fclose(ReadFileHandler2);
	}
	printf("读线程结束\n");


}




void   process_read_block(Executor<QueueItem_t *> *executor){

	int i =0,tail,head,pos =0;

	int chunkSize=ReadChunkSize;
	time_t	StartReadTime = time(NULL);



	if(gzCompressed){
		int ret1 =gzread(gzReadFileHandler1,tmp1,sizeof(char)* actual_block_size);
		tmp1[actual_block_size] ='\0';
		right_broundry1=actual_block_size;
		offset1+=actual_block_size;


		if(bSepLibrary){
			int ret2 =gzread(gzReadFileHandler2,tmp2, sizeof(char)* actual_block_size);
			tmp2[actual_block_size] ='\0';
			right_broundry2=actual_block_size;
			offset2+=actual_block_size;

		}	
	}else{
		int ret1=fread(tmp1,sizeof(char), actual_block_size,  ReadFileHandler1);

		tmp1[ret1]='\0';
		right_broundry1 =ret1;
		offset1+=ret1;

		if(bSepLibrary){
			int ret2 =fread(tmp2, sizeof(char), actual_block_size, ReadFileHandler2);

			tmp2[ret2]='\0';
			right_broundry2 =ret2;
			offset2+=ret2;

		}		
	}



    while(1) {
		QueueItem_t *queueItem;
		if(!spsc_queue_input.try_dequeue(queueItem)){
			// fprintf(stderr, "腾不出空间给读线程\n"); 
			continue;
		}else{
			if(gzCompressed){
				queueItem->readNum = gzGetNextChunk_block(bSepLibrary, gzReadFileHandler1, gzReadFileHandler2, queueItem->readArr,chunkSize,tail);

			} else{
				queueItem->readNum =GetNextChunk_Block_2(queueItem->readArr,bSepLibrary ,chunkSize,tail ,ReadFileHandler1,ReadFileHandler2);
			}
			
			// queueItem->readNum = GetNextChunk(queueItem->readArr,bSepLibrary ,chunkSize,tail,ReadFileHandler1,ReadFileHandler2);


			if (iPaired >= 1000){
				queueItem->EstDistance = (int)(iDistance / (iPaired >> 2));
				queueItem->EstDistance = queueItem->EstDistance + (queueItem->EstDistance >> 1);
			}else{
				queueItem->EstDistance = maxInsertSize;
			}  
			//最后一批

			if(queueItem->readNum==0){
				break;
			}else{
				queueItem->index =i++;
				executor->enqueue(queueItem);
				process_num.fetch_add(1);
				// fprintf(stderr, "读线程：data%d完成填装%d条序列  addr:%ld -%ld\n ",queueItem->index,queueItem->readNum , queueItem,queueItem->readArr);
			}

		}
    }
	free(tmp1);
	free(tmp2);

	if (gzCompressed){
		if (gzReadFileHandler1 != NULL) gzclose(gzReadFileHandler1);
		if (gzReadFileHandler2 != NULL) gzclose(gzReadFileHandler2);
	}else{
		if (ReadFileHandler1 != NULL) fclose(ReadFileHandler1);
		if (ReadFileHandler2 != NULL) fclose(ReadFileHandler2);
	}
	long long readTime =(long long)(time(NULL) - StartReadTime);
	fprintf(stderr, "\rAll the  reads have been read in %lld seconds.\n",readTime);

	// fclose(rfile_text);
	// precision, recall
	// FILE *rfile = fopen("io_exp_1.csv", "a");

	// fprintf(rfile, "%s: %s\t%d\t%.2lld\n",readFile1.c_str(),outputFile,threadNum, readTime);

	// fclose(rfile);


	printf("读线程结束\n");
}




void updateStatistic(SamOutputItem_t * samItem){

		iDistance +=samItem->myDistance;
		iPaired +=samItem->myPaired;
		iTotalReadNum += samItem->readNum; 
		iUniqueMapping += samItem->myUniqueMapping; 
		iUnMapping += samItem->myUnMapping;

		samItem->myDistance=0;
		samItem->myPaired=0;
		samItem->readNum=0;
		samItem->myUniqueMapping=0;
		samItem->myUnMapping=0;
		samItem->offset=0;
}

static bam_hdr_t *sam_hdr_sanitise(bam_hdr_t *h) 
{
	if (!h) return NULL;
	if (h->l_text == 0) return h;

	uint32_t i;
	char *cp = h->text;
	for (i = 0; i < h->l_text; i++) {
		// NB: l_text excludes terminating nul.  This finds early ones.
		if (cp[i] == 0) break;
	}
	if (i < h->l_text) { // Early nul found.  Complain if not just padding.
		uint32_t j = i;
		while (j < h->l_text && cp[j] == '\0') j++;
	}
	return h;
}

bam_hdr_t *SamHdr2BamHdr(kstring_t *str)
{
	bam_hdr_t *h = NULL;
	h = sam_hdr_parse(str->l, str->s);
	h->l_text = str->l; h->text = str->s;

	return sam_hdr_sanitise(h);
}


 
void writeHeader(FILE *fd){
	int len,i;
	char buffer[1024];
	kstring_t str = { 0, 0, NULL };


	len = sprintf(buffer, "@PG\tID:kart\tPN:Kart\tVN:%s\n", VersionStr);

	if (OutputFileFormat == 0) {
		// fwrite(buffer,len,1,fd);

		write(sam_out_fd, buffer, len);
		// fprintf(sam_out, "%s", buffer);
	}else{
		kputsn(buffer, len, &str);
	} 

		

	for (i = 0; i < iChromsomeNum; i++){
		len = sprintf(buffer, "@SQ\tSN:%s\tLN:%lld\n", ChromosomeVec[i].name, (long long)ChromosomeVec[i].len);
		if (OutputFileFormat == 0){
			// fwrite(buffer,len,1,fd);
			write(sam_out_fd, buffer, len);
			// memcpy((mmapped+offset), buffer, len);
			// offset+=len;
			// fprintf(sam_out, "%s", buffer);
		} else{
			kputsn(buffer, len, &str);
		} 
	}

	if (OutputFileFormat == 1){
		header = SamHdr2BamHdr(&str);
		sam_hdr_write(bam_out, header);
	}
}

void process_write(){
    int i,output_count=0 ,output_index =0,offset =0;
	int tmp_offset;
	// size_t block_offset =0 ,output_len =0,fwrite_size =65536;
	// char block_buff[fwrite_size];

	// memset( buff, '\0', sizeof( buff ));

    fprintf(stderr, "写线程\n");

    SamOutputItem_t *tempQueueItem;
	// string temp_res="";
	// int sam_fd=0;
    // if (OutputFileFormat == 0){
	// 	sam_out = fopen(outputFile, "w");
		// setvbuf(sam_out, block_buff, _IOFBF, fwrite_size);
	// } 
	// // sam_out = fopen(outputFile, "w");
	// if (sam_out == NULL ){
	// 	fprintf(stderr, "Error! Cannot open file [%s]\n", outputFile);
	// 	exit(1);
	// }

	// writeHeader(sam_out);
	

	SamOutputItem_t * minQueue;
   while(1) {
        if(output_info.shutdown && spsc_queue_output.size_approx()==0 && process_num.load()==0 ) { 
			// std::unique_lock<std::mutex> lock(run_mutex_t);
			run_cv_t.notify_one();
			output_info.runFlag =0;
			// lock.unlock();
            fprintf(stderr, "写线程退出循环\n"); 
            break;
        }

        // bool w =spsc_queue.pop(tempQueueItem);
        // QueueItem_t * temp;    

        if(!spsc_queue_output.try_dequeue(tempQueueItem)){ 
            continue;
        }else{
            // fprintf(stderr, "写线程:tempQueueItem->index:%d addr:%ld\n",tempQueueItem->index,tempQueueItem);
            priQueMinFirst.push(tempQueueItem);
            //9766
            while(priQueMinFirst.top()->index ==output_index ){
                minQueue =priQueMinFirst.top();
                priQueMinFirst.pop();

				tmp_offset=0;
				while(tmp_offset+fwrite_size <minQueue->offset){

					write(sam_out_fd, minQueue->samOutputBuffer+tmp_offset, fwrite_size);
					// fwrite(minQueue->samOutputBuffer+tmp_offset, fwrite_size , 1, sam_out );fflush(sam_out);

					tmp_offset+=fwrite_size;
				}
				write(sam_out_fd, minQueue->samOutputBuffer+tmp_offset, minQueue->offset - tmp_offset);
				// flush(sam_out_fd);
				// fwrite(minQueue->samOutputBuffer+tmp_offset, minQueue->offset - tmp_offset, 1, sam_out );
				// fwrite(minQueue->samOutputBuffer, minQueue->offset , 1, sam_out );

				// fsync(fileno(sam_out));	
				// fflush(sam_out);
				fprintf(stderr, "写线程:输出第%d(%d)的%d条reads 已接收:%d\n",output_index++,minQueue->index,minQueue->readNum,output_count);
				// if(output_index%1000==0){
				// 	fsync(fileno(sam_out));			
				// }
				updateStatistic(minQueue);
				spsc_queue_output_addr.enqueue(minQueue);

            }

			// fprintf(stderr,"fflush:%d\n",fflush(sam_out)) ;
            output_count++; 
        }
    }

    fprintf(stderr, "写线程:已输出%d 已接收：%d\n",output_index,output_count);
	// if (OutputFileFormat == 0) {
	// 	printf("%p %p\n",sam_out->_IO_buf_base,sam_out->_IO_buf_end);
	// 	// fsync(fileno(sam_out));
	// 	// munmap(mmapped, offset);
	// 	fprintf(stderr, "写线程:fclose\n");
	// 	fclose(sam_out);
	// }
}




void process_write_disorder(){
    int i,output_count=0 ,offset =0 ,output_read_num =0;;
	int total_processed_read =0;
	int tmp_offset=0;
	// int fwrite_size =65536;
	size_t blockSize =1048676;
	char* block_buff=(char*)malloc(blockSize);
 

    fprintf(stderr, "写线程启动\n");



    SamOutputItem_t *tempQueueItem;
   	while(1) {


        if(output_info.shutdown && spsc_queue_output.size_approx()==0 && process_num.load()==0 ) { 
			write(sam_out_fd, block_buff, tmp_offset );

			std::unique_lock<std::mutex> lock(run_mutex_t);
			run_cv_t.notify_one();
			output_info.runFlag =0;

            fprintf(stderr, "写线程退出循环\n");
            break;
        }else{
			if(!spsc_queue_output.try_dequeue(tempQueueItem)){ 
				continue;

			}else{

				output_count++;	
				output_read_num += tempQueueItem->readNum;
	
				if(tmp_offset +  tempQueueItem->offset >blockSize){
					while(tmp_offset +  tempQueueItem->offset >blockSize){
						blockSize*=2;
					}
					// printf("before sam_block_size:%d %d\n",blockSize,tmp_offset +  tempQueueItem->offset);
					block_buff =(char *)realloc(block_buff,blockSize);

				}				
				memcpy((block_buff)+tmp_offset, tempQueueItem->samOutputBuffer, tempQueueItem->offset);
				tmp_offset+=tempQueueItem->offset;


				if (OutputFileFormat == 0){

					if(output_read_num >=16384){
						write(sam_out_fd, block_buff, tmp_offset );
						tmp_offset =0;
						output_read_num=0;
					}

				}else{

					bam1_t *b = bam_init1();
					kstring_t str = { 0, 0, NULL };
					str.s = block_buff; 
					str.l = tmp_offset;
					if (sam_parse1(&str, header, b) >= 0){
						sam_write1(bam_out, header, b);
					} 
					bam_destroy1(b);
					tmp_offset =0;
					output_read_num=0;

				}					

					
					
					
					// write(sam_out_fd, tempQueueItem->samOutputBuffer, tempQueueItem->offset );
				if(output_count %1000 ==0){
					fprintf(stderr, "写线程:输出第%d(%d)的%d条reads\n",output_count,tempQueueItem->index,tempQueueItem->readNum);
				}

				// fprintf(stderr, "写线程:输出第%d(%d)的%d条reads\n",output_index++,tempQueueItem->index,total_processed_read+=tempQueueItem->readNum);
				updateStatistic(tempQueueItem);
				// fprintf(stderr,"spsc_queue_output.size_approx:%d process_num.load:%d\n",spsc_queue_output.size_approx(),process_num.load());

				// fflush(stderr);
				spsc_queue_output_addr.enqueue(tempQueueItem);
				process_num.fetch_add(-1);
				
			}
		}



    }
    fprintf(stderr, "写线程:已输出%d \n",output_count);

	free(block_buff);
}


vector<SeedPair_t> IdentifySeedPairs_SensitiveMode(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, stop_pos, end_pos;
	vector<SeedPair_t> SeedPairVec;
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; pos = 0, stop_pos = 30; end_pos = rlen - MinSeedLength;
	while (pos < end_pos)
	{
		if (EncodeSeq[pos] > 3) pos++, stop_pos++;
		else
		{
			bwtSearchResult = BWT_Search(EncodeSeq, pos, stop_pos);
			//if (bDebugMode) printf("Pos=%d, Freq=%d, Len=%d\n", pos, bwtSearchResult.freq, bwtSearchResult.len);
			if (bwtSearchResult.freq > 0)
			{
				SeedPair.rPos = pos; SeedPair.rLen = SeedPair.gLen = bwtSearchResult.len;
				for (i = 0; i != bwtSearchResult.freq; i++)
				{
					SeedPair.PosDiff = (SeedPair.gPos = bwtSearchResult.LocArr[i]) - SeedPair.rPos;
					SeedPairVec.push_back(SeedPair);
				}
				delete[] bwtSearchResult.LocArr;
				//pos += 30; stop_pos += 30;
				pos += bwtSearchResult.len; stop_pos += bwtSearchResult.len;
			}
			else
			{
				pos += MinSeedLength; stop_pos += MinSeedLength;
			}
			if (stop_pos > rlen) stop_pos = rlen;
		}
	}
	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByGenomePos);

	return SeedPairVec;
}


vector<AlignmentCandidate_t> GenerateAlignmentCandidateForPacBioSeq(int rlen, vector<SeedPair_t> SeedPairVec)
{
	bool *TakenArr;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	if ((num = (int)SeedPairVec.size()) > 0)
	{
		//thr = (int)(rlen*0.05);
		thr = 0;
		//printf("Raw seeds:\n"); ShowSeedInfo(SeedPairVec);
		AlignmentCandidate.PairedAlnCanIdx = -1; TakenArr = new bool[num]();
		i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
		for (; i < num; i++)
		{
			if (TakenArr[i]) continue;
			AlignmentCandidate.Score = SeedPairVec[i].rLen; TakenArr[i] = true;
			AlignmentCandidate.SeedVec.clear(); AlignmentCandidate.SeedVec.push_back(SeedPairVec[i]);
			//if (bDebugMode) printf("Master seed: r[%d-%d] g[%lld-%lld], len=%d\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen);
			for (j = i, k = i + 1; k < num; k++){
				if (TakenArr[k]) continue;
				if (abs(SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff) < 300){
					if (SeedPairVec[k].rPos > SeedPairVec[j].rPos)
					{
						//if (bDebugMode) printf("add seed: r[%d-%d] g[%lld-%lld], len=%d\n", SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen);
						AlignmentCandidate.Score += SeedPairVec[k].rLen;
						AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
						TakenArr[(j = k)] = true;
					}
				}
				else if (SeedPairVec[k].gPos - SeedPairVec[j].gPos > 1000) break;
			}
			if (AlignmentCandidate.Score >= thr)
			{
				thr = AlignmentCandidate.Score;
				AlignmentCandidate.PosDiff = SeedPairVec[i].PosDiff;
				if (AlignmentCandidate.PosDiff < 0) AlignmentCandidate.PosDiff = 0;

				AlignmentVec.push_back(AlignmentCandidate);
				//if (bDebugMode)
				//{
				//	printf("\n\nCandidate score = %d\n", AlignmentCandidate.Score);
				//	ShowSeedLocationInfo(AlignmentCandidate.PosDiff);
				//	ShowSeedInfo(AlignmentCandidate.SeedVec);
				//}
			}
		}
		delete[] TakenArr;
	}
	return AlignmentVec;
}


void EnCodeReadSeq(int rlen, char* seq, uint8_t* EncodeSeq ,int * gotN)
{
    int window =8;
    uint32_t final_feature=0 ;
	for (int i = 0; i < rlen; i++){
		if(nst_nt4_table[(int)seq[i]] >3){
			*gotN =1;
		}
       EncodeSeq[i] = nst_nt4_table[(int)seq[i]]; 
    }


    // float c = 0;
    // unsigned n = 0;
    // bool weighted = 0;
	// const std::set<jellyfish::mer_dna> & q =kmers_in_string(seq);
    // for ( auto & m : q) {
    //     // DEBUG: std::cout << "checking: " << m.to_str()<<endl;
    //     // if (bf->contains(m)) c++;
    //     //DEBUG: std::cout << c << std::endl;
	// 	// fprintf(stdout,"%s\n",m.to_str().c_str());
    //     if (bf->contains(m)){
	// 		n++;
	// 	} 
	// 	// c+=weight;
    // }
    // // // std::cerr << root->name() <<std::endl;
    // // std::cerr << n << " " << QUERY_THRESHOLD << " " << q.size()<<" " <<(float)n/q.size() << std::endl;
	// int wdnmd =(int)q.size();
	// float nmsl= (float)n/q.size();

}


void EnCodeReadSeq_gotN(int rlen, char* seq, uint8_t* EncodeSeq ,int & gotN)
{
    int window =8;
    uint32_t final_feature=0 ;
	for (int i = 0; i < rlen; i++){
		if(nst_nt4_table[(int)seq[i]] >3){
			gotN =1;
		}
       EncodeSeq[i] = nst_nt4_table[(int)seq[i]]; 
    }
}



void EnCodeReadSeq_Pair(int rlen, char* seq, uint8_t* EncodeSeq ,int rlen_p, char* seq_p, uint8_t* EncodeSeq_p)
{

	for (int i = 0; i < rlen; i++){
       EncodeSeq[i] = nst_nt4_table[(int)seq[i]]; 
    }
	for (int i = 0; i < rlen_p; i++){
       EncodeSeq_p[i] = nst_nt4_table[(int)seq_p[i]]; 
    }
    // for(int i =0;i+window-1<rlen;i++){
    //     uint32_t kmerNumber=0 ;
    //     for(int j=0;j<window;j++){

    //         kmerNumber =(kmerNumber<<2)|EncodeSeq[i+j];
    //     }
    //     final_feature +=kmerNumber;
    // }
    // *feature =final_feature;
}


void RemoveRedundantCandidates(vector<AlignmentCandidate_t>& AlignmentVec)
{
	int thr;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return;
	else
	{
		int score1, score2;

		score1 = score2 = 0;
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
		{
			if (iter->Score > score2)
			{
				if (iter->Score >= score1)
				{
					score2 = score1;
					score1 = iter->Score;
				}
				else score2 = iter->Score;
			}
		}
		if (bPacBioData || score1 == score2 || score1 - score2 > 20) thr = score1;
		else thr = score2;

		//if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
		// for(int x=0;x<AlignmentVec.size();x++){
		// 	fprintf(stderr,"score:%d\n",AlignmentVec[x].Score);
		// }
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++) if (iter->Score < thr) iter->Score = 0;
	}
}

void SetSingleAlignmentFlag(ReadItem_t& read)
{
	int i;

	if (read.score > read.sub_score) // unique mapping
	{
		i = read.iBestAlnCanIdx;
		if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].SamFlag = 0x10;
		else read.AlnReportArr[i].SamFlag = 0;
	}else if(read.score > 0){
		for (i = 0; i < read.CanNum; i++){
			if (read.AlnReportArr[i].AlnScore > 0)
			{
				if (read.AlnReportArr[i].coor.bDir == false) read.AlnReportArr[i].SamFlag = 0x10;
				else read.AlnReportArr[i].SamFlag = 0;
			}
		}
	}else {
		read.AlnReportArr[0].SamFlag = 0x4;
	}


	if (read.score == 0 || read.score == read.sub_score) {
		read.mapq = 0;
	}
	else{
		if (bPacBioData)
		{
			float fScale = 85.0*(int)(ceil(read.rlen / 100 + 0.5));
			if (fScale > 2000) fScale = 2000;
			read.mapq = (int)(Max_MAPQ * (read.score / fScale));
		}
		else if (read.sub_score == 0 || read.score - read.sub_score > 5) read.mapq = Max_MAPQ;
		else read.mapq = (int)(MAPQ_COEF * (1 - (float)(read.score - read.sub_score) / read.score)*log(read.score) + 0.4999);
		if (read.mapq > Max_MAPQ) read.mapq = Max_MAPQ;
	}

}

void EvaluateMAPQ(ReadItem_t& read)
{
	if (read.score == 0 || read.score == read.sub_score) read.mapq = 0;
	else
	{
		if (bPacBioData)
		{
			float fScale = 85.0*(int)(ceil(read.rlen / 100 + 0.5));
			if (fScale > 2000) fScale = 2000;
			read.mapq = (int)(Max_MAPQ * (read.score / fScale));
		}
		else if (read.sub_score == 0 || read.score - read.sub_score > 5) read.mapq = Max_MAPQ;
		else read.mapq = (int)(MAPQ_COEF * (1 - (float)(read.score - read.sub_score) / read.score)*log(read.score) + 0.4999);
		if (read.mapq > Max_MAPQ) read.mapq = Max_MAPQ;
	}
}

void OutputSingledAlignments(ReadItem_t& read, char* buffer, int& myUniqueMapping, int& myUnMapping, vector<string>& SamOutputVec)
{
	int len;
	string rqual;

	if (read.score == 0){
		myUnMapping++;
		len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read.header, read.AlnReportArr[0].SamFlag, read.seq, (FastQFormat ? read.qual : "*"));
		buffer[len] = '\0'; SamOutputVec.push_back(buffer);
	}
	else
	{
		int i;
		char *seq, *rseq;

		if (read.mapq == Max_MAPQ) myUniqueMapping++;

		seq = read.seq; rseq = NULL;
		for (i = read.iBestAlnCanIdx; i < read.CanNum; i++)
		{
			if (read.AlnReportArr[i].AlnScore == read.score){
				if (read.AlnReportArr[i].coor.bDir == false && rseq == NULL)
				{
					rseq = new char[read.rlen + 1]; rseq[read.rlen] = '\0'; GetComplementarySeq(read.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read.header, read.AlnReportArr[i].SamFlag, ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read.AlnReportArr[i].coor.gPos, read.mapq, read.AlnReportArr[i].coor.CIGAR.c_str(), (read.AlnReportArr[i].coor.bDir? seq: rseq), (FastQFormat ? (read.AlnReportArr[i].coor.bDir ? read.qual : rqual.c_str()) : "*"), read.rlen - read.score, read.score, read.sub_score);
				buffer[len] = '\0'; SamOutputVec.push_back(buffer);

				if (!bMultiHit) break;
			}
		}
		if (rseq != NULL)
		{
			delete[] rseq;
			rseq = NULL;
		}
	}
}

void OutputSingledAlignments_block(ReadItem_t& read, char* buffer, SamOutputItem_t * & samItem){

	// printf("OutputSingledAlignments_block:%p\n",samItem);
	int line_len;
	string rqual;
	size_t tmp_offset =samItem->offset;
	SamOutputItem_t * samItem_New;
	if (read.score == 0){
		samItem->myUnMapping++;
		line_len = sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read.header, read.AlnReportArr[0].SamFlag, read.seq, (FastQFormat ? read.qual : "*"));
		// buffer[len] = '\0';SamOutputVec.push_back(buffer);

		if(samItem->offset +  line_len >output_block_size){

			while (!spsc_queue_output_addr.try_dequeue(samItem_New)){
			//    cout << " try直到获取到内存空间! " << endl;
			}
			// fprintf(stderr,"samItem_New:%p\n",samItem_New);
			samItem_New->offset =samItem->offset-tmp_offset;
			samItem_New->index=samItem->index;
			samItem_New->localIndex=samItem->localIndex+samItem->readNum;
			memcpy((samItem_New->samOutputBuffer), (samItem->samOutputBuffer)+tmp_offset, samItem_New->offset);


			samItem->offset =tmp_offset;
			spsc_queue_output.enqueue(samItem); 
			// fprintf(stderr,"samItem_Old:%p %d %d %d\n",samItem ,samItem->index,samItem->localIndex,samItem->readNum);

			samItem =samItem_New;
			tmp_offset =0;
		}				
		memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
		samItem->offset+=line_len;
	}else{
		int i;
		char *seq, *rseq;

		if (read.mapq == Max_MAPQ) samItem->myUniqueMapping++;

		seq = read.seq; rseq = NULL;
		for (i = read.iBestAlnCanIdx; i < read.CanNum; i++){
			if (read.AlnReportArr[i].AlnScore == read.score){
				if (read.AlnReportArr[i].coor.bDir == false && rseq == NULL){
					rseq = new char[read.rlen + 1]; rseq[read.rlen] = '\0'; GetComplementarySeq(read.rlen, seq, rseq);
					if (FastQFormat){
						rqual = read.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				line_len = sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read.header, read.AlnReportArr[i].SamFlag, ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read.AlnReportArr[i].coor.gPos, read.mapq, read.AlnReportArr[i].coor.CIGAR.c_str(), (read.AlnReportArr[i].coor.bDir? seq: rseq), (FastQFormat ? (read.AlnReportArr[i].coor.bDir ? read.qual : rqual.c_str()) : "*"), read.rlen - read.score, read.score, read.sub_score);
				// buffer[len] = '\0'; SamOutputVec.push_back(buffer);
				if(samItem->offset +  line_len >output_block_size){

					while (!spsc_queue_output_addr.try_dequeue(samItem_New)){
					//    cout << " try直到获取到内存空间! " << endl;
					}
					// fprintf(stderr,"samItem_New:%p\n",samItem_New);
					samItem_New->offset =samItem->offset-tmp_offset;
					samItem_New->index=samItem->index;
					samItem_New->localIndex=samItem->localIndex+samItem->readNum;
					memcpy((samItem_New->samOutputBuffer), (samItem->samOutputBuffer)+tmp_offset, samItem_New->offset);


					samItem->offset =tmp_offset;
					spsc_queue_output.enqueue(samItem); 
					// fprintf(stderr,"samItem_Old:%p %d %d %d\n",samItem ,samItem->index,samItem->localIndex,samItem->readNum);

					samItem =samItem_New;
					tmp_offset =0;
				}				
				memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
				samItem->offset+=line_len;

				if (!bMultiHit) break;
			}
		}
		if (rseq != NULL){
			delete[] rseq;
			rseq = NULL;
		}
	}

	samItem->readNum++;
}


void CheckPairedFinalAlignments(ReadItem_t& read1, ReadItem_t& read2)
{
	bool bMated;
	int i, j, s;

	//printf("BestIdx1=%d, BestIdx2=%d\n", read1.iBestAlnCanIdx + 1, read2.iBestAlnCanIdx + 1);
	if (read1.iBestAlnCanIdx != -1 && read2.iBestAlnCanIdx != -1){
		bMated = read1.AlnReportArr[read1.iBestAlnCanIdx].PairedAlnCanIdx == read2.iBestAlnCanIdx ? true : false;
	} 
	else bMated = false;

	if (!bMultiHit && bMated) return;
	if (!bMated && read1.score > 0 && read2.score > 0) // identify mated pairs
	{	
		//找最高分的mate pair
		//score和iBestAlnCanIdx都是最高分mate pair的
		for (s = 0, i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore > 0 && (j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0)
			{
				bMated = true;
				if (s < read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore)
				{
					s = read1.AlnReportArr[i].AlnScore + read2.AlnReportArr[j].AlnScore;
					read1.iBestAlnCanIdx = i; read1.score = read1.AlnReportArr[i].AlnScore;
					read2.iBestAlnCanIdx = j; read2.score = read2.AlnReportArr[j].AlnScore;
				}
			}
		}
	}
	//其他的应该保留prime和secondary
	if(bMated)
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].AlnScore != read1.score || ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore != read2.score))
			{	
				if(!read1.AlnReportArr[i].isPrime && !read1.AlnReportArr[i].isSecondary){
					read1.AlnReportArr[i].AlnScore = 0;
					read1.AlnReportArr[i].PairedAlnCanIdx = -1;
					continue;					
				}

			}
		}
	}
	else // remove all mated info
	{
		for (i = 0; i != read1.CanNum; i++)
		{
			if (read1.AlnReportArr[i].PairedAlnCanIdx != -1) read1.AlnReportArr[i].PairedAlnCanIdx = -1;
			if (read1.AlnReportArr[i].AlnScore > 0 && read1.AlnReportArr[i].AlnScore != read1.score){
				if(!read1.AlnReportArr[i].isPrime && !read1.AlnReportArr[i].isSecondary){
					read1.AlnReportArr[i].AlnScore = 0;
				}
			} 
		}
		for (j = 0; j != read2.CanNum; j++)
		{
			if (read2.AlnReportArr[j].PairedAlnCanIdx != -1) read2.AlnReportArr[j].PairedAlnCanIdx = -1;
			if (read2.AlnReportArr[j].AlnScore > 0 && read2.AlnReportArr[j].AlnScore != read2.score){
				if(!read2.AlnReportArr[j].isPrime && !read2.AlnReportArr[j].isSecondary){
					// read1.AlnReportArr[i].AlnScore = 0;
					read2.AlnReportArr[j].AlnScore = 0;
				}
			} 
		}
	}
}


bool CheckPairedAlignmentCandidates(int64_t EstiDistance, vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int64_t dist;
	bool bPairing = false;
	int i, j, best_mate, s, num1, num2;

	num1 = (int)AlignmentVec1.size(); 
	num2 = (int)AlignmentVec2.size();

	// if (num1*num2 > 1000)
	// {
	// 	RemoveRedundantCandidates(AlignmentVec1);
	// 	RemoveRedundantCandidates(AlignmentVec2);
	// }

	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].Score == 0) continue;

		for (best_mate = -1, s = 0, j = 0; j != num2; j++)
		{
			if (AlignmentVec2[j].Score == 0 || AlignmentVec2[j].PosDiff < AlignmentVec1[i].PosDiff) continue;

			dist = AlignmentVec2[j].PosDiff - AlignmentVec1[i].PosDiff;
			//printf("#%d (s=%d) and #%d (s=%d) (dist=%lld / %d)\n", i+1, AlignmentVec1[i].Score, j+1, AlignmentVec2[j].Score, dist, EstiDistance), fflush(stderr);
			if (dist < EstiDistance)
			{
				if (AlignmentVec2[j].Score > s)
				{
					best_mate = j;
					s = AlignmentVec2[j].Score;
				}else if (AlignmentVec2[j].Score == s){
					best_mate = -1;
				} 
			}
		}
		//best_mate是有唯一的mate
		if (s > 0 && best_mate != -1)
		{
			j = best_mate;
			if (AlignmentVec2[j].PairedAlnCanIdx == -1){
				bPairing = true;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}else if (AlignmentVec1[i].Score > AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].Score){
				AlignmentVec1[AlignmentVec2[j].PairedAlnCanIdx].PairedAlnCanIdx = -1;
				AlignmentVec1[i].PairedAlnCanIdx = j;
				AlignmentVec2[j].PairedAlnCanIdx = i;
			}
		}
	}
	return bPairing;
}


void RemoveUnMatedAlignmentCandidates(vector<AlignmentCandidate_t>& AlignmentVec1, vector<AlignmentCandidate_t>& AlignmentVec2)
{
	int i, j, num1, num2;

	num1 = (int)AlignmentVec1.size(); 
	num2 = (int)AlignmentVec2.size();

	for (i = 0; i != num1; i++)
	{
		if (AlignmentVec1[i].PairedAlnCanIdx == -1){
			AlignmentVec1[i].Score = 0;
		} 
		else{
			j = AlignmentVec1[i].PairedAlnCanIdx;
			AlignmentVec1[i].Score = AlignmentVec2[j].Score = AlignmentVec1[i].Score + AlignmentVec2[j].Score;
		}
	}
	for (j = 0; j != num2; j++){
		if (AlignmentVec2[j].PairedAlnCanIdx == -1){
			AlignmentVec2[j].Score = 0;
		} 
	} 

	if (bDebugMode)
	{
		for (i = 0; i != num1; i++)
		{
			if ((j = AlignmentVec1[i].PairedAlnCanIdx) != -1)
				printf("#%d(s=%d) and #%d(s=%d) are pairing\n", i + 1, AlignmentVec1[i].Score, j + 1, AlignmentVec2[j].Score);
		}
	}
}

void SetPairedAlignmentFlag(ReadItem_t& read1, ReadItem_t& read2)
{
	int i, j;

	//printf("read1:[%d, %d]:#%d, read2:[%d, %d] #%d\n", read1.score, read1.sub_score, read1.iBestAlnCanIdx, read2.score, read2.sub_score, read2.iBestAlnCanIdx); fflush(stderr);
	if (read1.score > read1.sub_score && read2.score > read2.sub_score) // unique mapping
	{
		i = read1.iBestAlnCanIdx;
		read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair

		j = read2.iBestAlnCanIdx;
		read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair

		if (j == read1.AlnReportArr[i].PairedAlnCanIdx) // reads are mapped in a proper pair
		{
			read1.AlnReportArr[i].SamFlag |= 0x2;
			read2.AlnReportArr[j].SamFlag |= 0x2;
		}
		read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
		read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);
	}
	else
	{
		if (read1.score > read1.sub_score) // unique mapping or bMultiHit=false
		{
			i = read1.iBestAlnCanIdx;
			read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);
			if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0){
				read1.AlnReportArr[i].SamFlag |= 0x2;// reads are mapped in a proper pair
			} 
			else{
				read1.AlnReportArr[i].SamFlag |= 0x8; // next segment unmapped
			} 
		}
		else if(read1.score > 0)
		{
			for (i = 0; i < read1.CanNum; i++)
			{
				if (read1.AlnReportArr[i].AlnScore > 0){
					read1.AlnReportArr[i].SamFlag = 0x41; // read1 is the first read in a pair
					read1.AlnReportArr[i].SamFlag |= (read1.AlnReportArr[i].coor.bDir ? 0x20 : 0x10);

					if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0){
						read1.AlnReportArr[i].SamFlag |= 0x2;// reads are mapped in a proper pair
					} 
					else{
						read1.AlnReportArr[i].SamFlag |= 0x8; // next segment unmapped
					} 
				}
			}
		}
		else
		{
			read1.AlnReportArr[0].SamFlag = 0x41; // read1 is the first read in a pair
			read1.AlnReportArr[0].SamFlag |= 0x4;

			if (read2.score == 0) read1.AlnReportArr[0].SamFlag |= 0x8; // next segment unmapped
			else read1.AlnReportArr[0].SamFlag |= (read2.AlnReportArr[read2.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}

		if (read2.score > read2.sub_score) // unique mapping or bMultiHit=false
		{
			j = read2.iBestAlnCanIdx;
			read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

			if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].SamFlag |= 0x2;// reads are mapped in a proper pair
			else read2.AlnReportArr[j].SamFlag |= 0x8; // next segment unmapped
		}
		else if (read2.score > 0)
		{
			for (j = 0; j < read2.CanNum; j++)
			{
				if (read2.AlnReportArr[j].AlnScore > 0)
				{
					read2.AlnReportArr[j].SamFlag = 0x81; // read2 is the second read in a pair
					read2.AlnReportArr[j].SamFlag |= (read2.AlnReportArr[j].coor.bDir ? 0x20 : 0x10);

					if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0) read2.AlnReportArr[j].SamFlag |= 0x2;// reads are mapped in a proper pair
					else read2.AlnReportArr[j].SamFlag |= 0x8; // next segment unmapped
				}
			}
		}
		else
		{
			read2.AlnReportArr[0].SamFlag = 0x81; // read2 is the second read in a pair
			read2.AlnReportArr[0].SamFlag |= 0x4; // segment unmapped
			if (read1.score == 0) read2.AlnReportArr[0].SamFlag |= 0x8; // next segment unmapped
			else read2.AlnReportArr[0].SamFlag |= (read1.AlnReportArr[read1.iBestAlnCanIdx].coor.bDir ? 0x10 : 0x20);
		}
	}
	for(i = 0; i < read1.CanNum; i++){
		if(!read1.AlnReportArr[i].isPrime){
			read1.AlnReportArr[i].SamFlag |= 0x100;
		}
	}
	for(j = 0; j < read2.CanNum; j++){
		if(!read2.AlnReportArr[j].isPrime){
			read2.AlnReportArr[j].SamFlag |= 0x100;
		}
	}								


}

void OutputPairedAlignments(ReadItem_t& read1, ReadItem_t& read2, char* buffer, int& myUniqueMapping, int& myUnMapping,int& myDistance, int& myPaired, vector<string>& SamOutputVec ,vector<string>& SamOutputVec1)
{
	string rqual;
	char *seq, *rseq;
	int i, j, dist = 0;

	if (read1.score == 0){
		myUnMapping++;
		(void)sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read1.header, read1.AlnReportArr[0].SamFlag, read1.seq, (FastQFormat ? read1.qual : "*"));
		SamOutputVec.push_back(buffer);
	}else{
		if (read1.mapq == Max_MAPQ) myUniqueMapping++;
		seq = read1.seq; rseq = NULL;
		for (i = read1.iBestAlnCanIdx; i < read1.CanNum; i++){
			if (read1.AlnReportArr[i].AlnScore > 0){
				if (read1.AlnReportArr[i].coor.bDir == false && rseq == NULL){
					rseq = new char[read1.rlen + 1]; rseq[read1.rlen] = '\0'; GetComplementarySeq(read1.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read1.qual; 
						reverse(rqual.begin(), rqual.end());
					}
				}
				if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0){
					dist = (int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen));
					if (i == read1.iBestAlnCanIdx){
						myPaired += 2;
						if (abs(dist) < 10000){
							myDistance += abs(dist);
						} 
					}
					(void)sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (long long)read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);
				}else{
					(void)sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);
				}
				SamOutputVec.push_back(buffer);
			}
			if (!bMultiHit) break;
		}
		if (rseq != NULL){
			delete[] rseq;
			rseq = NULL;
		}
	}

	if (read2.score == 0){
		myUnMapping++;
		(void)sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read2.header, read2.AlnReportArr[0].SamFlag, read2.seq, (FastQFormat ? read2.qual : "*"));
		SamOutputVec1.push_back(buffer);
	}else{
		if (read2.mapq == Max_MAPQ) myUniqueMapping++;

		rseq = read2.seq; seq = NULL;
		for (j = read2.iBestAlnCanIdx; j < read2.CanNum; j++){
			if (read2.AlnReportArr[j].AlnScore > 0)
			{
				if (read2.AlnReportArr[j].coor.bDir == true && seq == NULL){
					seq = new char[read2.rlen + 1]; 
					seq[read2.rlen] = '\0'; 
					GetComplementarySeq(read2.rlen, rseq, seq);
					if (FastQFormat){
						rqual = read2.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0){
					dist = 0 - ((int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen)));
					(void)sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (long long)read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);
				}
				else (void)sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);
				SamOutputVec1.push_back(buffer);
			}
			if (!bMultiHit) break;
		}
		if (seq != NULL){
			delete[] seq;
			seq = NULL;
		}
	}
}

int roundup(int x){
	--(x);
	(x)|=(x)>>1;
	(x)|=(x)>>2;
	(x)|=(x)>>4;
	(x)|=(x)>>8;
	(x)|=(x)>>16;
	++(x);
	// (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x));
	return x;
}


void OutputPairedAlignments_block(ReadItem_t& read1, ReadItem_t& read2, char* buffer, SamOutputItem_t * & samItem){
	string rqual;
	char *seq, *rseq;
	int i, j, dist = 0;
	int line_len =0;
	size_t tmp_offset =samItem->offset ;

	char *xa =(char *)malloc(1024);
	
	int prime_has_output =0;

	if (read1.score == 0){
		samItem->myUnMapping++;	
		line_len=sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read1.header, read1.AlnReportArr[0].SamFlag, read1.seq, (FastQFormat ? read1.qual : "*"));
		if(samItem->offset +  line_len >samItem->sam_block_size){
			samItem->sam_block_size*=2;
			samItem->samOutputBuffer =(char *)realloc(samItem->samOutputBuffer,samItem->sam_block_size);

		}				
		memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
		samItem->offset+=line_len;
		// SamOutputVec.push_back(buffer);
	}else{
		if (read1.mapq == Max_MAPQ) samItem->myUniqueMapping++;
		seq = read1.seq; rseq = NULL;
		string XA_Prime = "XA:Z:";
		char * for_or_back="-+";
		for(int x =0;x<read1.CanNum;x++){
			if(read1.AlnReportArr[x].isPrime && x!=read1.iBestAlnCanIdx){
				int xa_len= sprintf(xa, "%s,%c%lld,%s,%d;", ChromosomeVec[read1.AlnReportArr[x].coor.ChromosomeIdx].name,"-+"[read1.AlnReportArr[x].coor.bDir],(long long)read1.AlnReportArr[x].coor.gPos,read1.AlnReportArr[x].coor.CIGAR.c_str(),read1.AlnReportArr[x].NM);
				// strcat(XA_Prime, xa);
				XA_Prime.append(xa);
			}
		}				
		// printf("%d,supplement:%s\n",XA_Prime.size(),XA_Prime.c_str());

		// printf("%s\n",XA_Prime);

		for (i = read1.iBestAlnCanIdx; i < read1.CanNum; i++){
			if(prime_has_output && read1.AlnReportArr[i].isPrime){
				continue;
			}
			if (read1.AlnReportArr[i].AlnScore > 0){


				if (read1.AlnReportArr[i].coor.bDir == false && rseq == NULL){
					rseq = new char[read1.rlen + 1]; rseq[read1.rlen] = '\0'; GetComplementarySeq(read1.rlen, seq, rseq);
					if (FastQFormat)
					{
						rqual = read1.qual; 
						reverse(rqual.begin(), rqual.end());
					}
				}
				if ((j = read1.AlnReportArr[i].PairedAlnCanIdx) != -1 && read2.AlnReportArr[j].AlnScore > 0){
					dist = (int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen));
					if (i == read1.iBestAlnCanIdx){
						samItem->myPaired += 2;
						if (abs(dist) < 10000){
							samItem->myDistance += abs(dist);
						} 
					}
					if(XA_Prime.size()>5){
						line_len= sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\t%s\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (long long)read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score,XA_Prime.c_str());
	
					}else{
						line_len= sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (long long)read2.AlnReportArr[j].coor.gPos, dist, (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);

					}
				}else{
					if(XA_Prime.size()>5){
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\t%s\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score,XA_Prime.c_str());

					}else{
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read1.header, read1.AlnReportArr[i].SamFlag, ChromosomeVec[read1.AlnReportArr[i].coor.ChromosomeIdx].name, (long long)read1.AlnReportArr[i].coor.gPos, read1.mapq, read1.AlnReportArr[i].coor.CIGAR.c_str(), (read1.AlnReportArr[i].coor.bDir ? seq : rseq), (FastQFormat ? (read1.AlnReportArr[i].coor.bDir ? read1.qual : rqual.c_str()) : "*"), read1.rlen - read1.score, read1.score, read1.sub_score);

					}
				}


				if(samItem->offset +  line_len >samItem->sam_block_size){
					samItem->sam_block_size*=2;
					samItem->samOutputBuffer =(char *)realloc(samItem->samOutputBuffer,samItem->sam_block_size);

				}				
				memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
				samItem->offset+=line_len;
			}
			if (!bMultiHit){
				break;
			} else{
				prime_has_output=1;
				XA_Prime = "XA:Z:";
			}
		}
		if (rseq != NULL){
			delete[] rseq;
			rseq = NULL;
		}
	}
	prime_has_output =0;
	if (read2.score == 0){
		samItem->myUnMapping++;
		line_len=sprintf(buffer, "%s\t%d\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tAS:i:0\tXS:i:0\n", read2.header, read2.AlnReportArr[0].SamFlag, read2.seq, (FastQFormat ? read2.qual : "*"));
		// SamOutputVec1.push_back(buffer);
		if(samItem->offset +  line_len >samItem->sam_block_size){
			// printf("before sam_block_size:%d %p\n",samItem->sam_block_size,samItem->samOutputBuffer);
			samItem->sam_block_size*=2;
			samItem->samOutputBuffer =(char *)realloc(samItem->samOutputBuffer,samItem->sam_block_size);

		}				
		memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
		samItem->offset+=line_len;

	}else{
		if (read2.mapq == Max_MAPQ) samItem->myUniqueMapping++;

		rseq = read2.seq; seq = NULL;
		string XA_Prime = "XA:Z:";
		for(int x =0;x<read2.CanNum;x++){
			if(read2.AlnReportArr[x].isPrime && x!=read2.iBestAlnCanIdx){
				int xa_len= sprintf(xa, "%s,%c%lld,%s,%d;", ChromosomeVec[read2.AlnReportArr[x].coor.ChromosomeIdx].name,"-+"[read2.AlnReportArr[x].coor.bDir],(long long)read2.AlnReportArr[x].coor.gPos,read2.AlnReportArr[x].coor.CIGAR.c_str(),read2.AlnReportArr[x].NM);
				
				XA_Prime.append(xa);
				// strcat(XA_Prime, xa);
			}
		}
		// printf("%d,supplement:%s\n",XA_Prime.size(),XA_Prime.c_str());
		for (j = read2.iBestAlnCanIdx; j < read2.CanNum; j++){
			if(prime_has_output && read2.AlnReportArr[j].isPrime){
				continue;
			}
			if (read2.AlnReportArr[j].AlnScore > 0){

				if (read2.AlnReportArr[j].coor.bDir == true && seq == NULL){
					seq = new char[read2.rlen + 1]; 
					seq[read2.rlen] = '\0'; 
					GetComplementarySeq(read2.rlen, rseq, seq);
					if (FastQFormat){
						rqual = read2.qual; reverse(rqual.begin(), rqual.end());
					}
				}
				if ((i = read2.AlnReportArr[j].PairedAlnCanIdx) != -1 && read1.AlnReportArr[i].AlnScore > 0){
					dist = 0 - ((int)(read2.AlnReportArr[j].coor.gPos - read1.AlnReportArr[i].coor.gPos + (read1.AlnReportArr[i].coor.bDir ? read2.rlen : 0 - read1.rlen)));
					if(XA_Prime.size()>5){
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\t%s\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (long long)read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score,XA_Prime.c_str());

					}else{
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t=\t%lld\t%d\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (long long)read1.AlnReportArr[i].coor.gPos, dist, (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);

					}
				}
				else{
					if(XA_Prime.size()>5){
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\t%s\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score,XA_Prime.c_str());

					}else{
						line_len=sprintf(buffer, "%s\t%d\t%s\t%lld\t%d\t%s\t*\t0\t0\t%s\t%s\tNM:i:%d\tAS:i:%d\tXS:i:%d\n", read2.header, read2.AlnReportArr[j].SamFlag, ChromosomeVec[read2.AlnReportArr[j].coor.ChromosomeIdx].name, (long long)read2.AlnReportArr[j].coor.gPos, read2.mapq, read2.AlnReportArr[j].coor.CIGAR.c_str(), (read2.AlnReportArr[j].coor.bDir ? seq : rseq), (FastQFormat ? (read2.AlnReportArr[j].coor.bDir ? rqual.c_str() : read2.qual) : "*"), read2.rlen - read2.score, read2.score, read2.sub_score);

					}
				} 
				if(samItem->offset +  line_len >samItem->sam_block_size){
					// printf("before sam_block_size:%d %p\n",samItem->sam_block_size,samItem->samOutputBuffer);
					samItem->sam_block_size*=2;
					samItem->samOutputBuffer =(char *)realloc(samItem->samOutputBuffer,samItem->sam_block_size);
				}				
				memcpy((samItem->samOutputBuffer)+samItem->offset, buffer, line_len);
				samItem->offset+=line_len;
			}
			if (!bMultiHit){

				break;
			} else{
				prime_has_output =1;
				XA_Prime = "XA:Z:";
			}
		}
		if (seq != NULL){
			delete[] seq;
			seq = NULL;
		}
	}

	samItem->readNum+=2;
	free(xa);
}





void  process_mapping( void * arg ) {
    char* buffer;
    QueueItem_t *queueItem =(QueueItem_t *)arg;

 
	// fprintf(stderr,"子线程接收\n");

    SamOutputItem_t *samItem ;
    while (!spsc_queue_output_addr.try_dequeue(samItem)){
    //    cout << " try直到获取到内存空间! " << endl;
    }
	samItem->index =queueItem->index;
	samItem->localIndex =0;

    int threadId =pthread_self() ;
	bool bReadPairing;   

	int i, j;


	QuerySet qs;

    if (bPacBioData) buffer = new char[1024000];
	else buffer = new char[10240];


	if (bPacBioData){
		for (i = 0; i != queueItem->readNum; i++){
			if (bDebugMode) printf("\n\n\nMapping pacbio read#%d %s (len=%d):\n", i + 1, queueItem->readArr[i].header, queueItem->readArr[i].rlen);
			
			queueItem->EncodeSeq[i] = new uint8_t[queueItem->readArr[i].rlen]; 
			EnCodeReadSeq(queueItem->readArr[i].rlen, queueItem->readArr[i].seq, queueItem->EncodeSeq[i],&queueItem->readArr[i].gotN);
			// queueItem->SeedPairVec1[i] = IdentifySeedPairs_SensitiveMode(queueItem->readArr[i].rlen, queueItem->EncodeSeq[i]); delete[] queueItem->EncodeSeq[i];

			// queueItem->AlignmentVec1[i] = GenerateAlignmentCandidateForPacBioSeq(queueItem->readArr[i].rlen, queueItem->SeedPairVec1[i]);
			//if (bDebugMode) ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
			RemoveRedundantCandidates(queueItem->AlignmentVec1[i]);
			if (bDebugMode) ShowAlignmentCandidateInfo(1, queueItem->readArr[i].header, queueItem->AlignmentVec1[i]);
			GenMappingReport(true, queueItem->readArr[i], queueItem->AlignmentVec1[i]);

			SetSingleAlignmentFlag(queueItem->readArr[i]); 
			// EvaluateMAPQ(queueItem->readArr[i]);
			//if (bDebugMode) printf("\nEnd of mapping for read#%s (len=%d)\n%s\n", ReadArr[i].header, ReadArr[i].rlen, string().assign(100, '=').c_str());
		    OutputSingledAlignments(queueItem->readArr[i], buffer, queueItem->myUniqueMapping, queueItem->myUnMapping, queueItem->SamOutputVec[i]);            
		
		}

	}else if (bPairEnd && queueItem->readNum % 2 == 0){


		// for (auto& q : qs) {
		// 	out << "*" << q->query << " " << q->matching.size() << std::endl;
		// 	for (const auto& n : q->matching) {
		// 		out << n->name() << std::endl;
		// 	}
		// }
		// print_query_results(qs, o);
		

        for(i=0;i<queueItem->readNum;i++){
            queueItem->EncodeSeq[i] = new uint8_t[queueItem->readArr[i].rlen];
            EnCodeReadSeq(queueItem->readArr[i].rlen, queueItem->readArr[i].seq, queueItem->EncodeSeq[i],&queueItem->readArr[i].gotN);
			qs.emplace_back(new QueryInfo(queueItem->readArr[i].seq));
        }
		//在这里植入布隆过滤器

		query_batch(root, qs);

		for (auto& q : qs) {
			fprintf(bffile, "%d\t%f\t",q->matching.size(),QUERY_THRESHOLD );
			for (const auto& n : q->matching) {
				fprintf(bffile, "%s\t",n->name().c_str() );
			}
			fprintf(bffile, "\n" );
		}


		for(i=0;i<queueItem->readNum;i++){
			IdentifySeedPairs_FastMode_getN(queueItem->readArr[i].rlen, queueItem->EncodeSeq[i],queueItem->SeedPairVec1[i],queueItem->readArr[i].gotN);
			GenerateAlignmentCandidateForIlluminaSeq_Recycle(queueItem->readArr[i].rlen, queueItem->SeedPairVec1[i],queueItem->AlignmentVec1[i]);
			
		}
		for(i=0;i<queueItem->readNum;i++){
			delete[] queueItem->EncodeSeq[i];
		}



		for (i = 0, j = 1; i != queueItem->readNum; i += 2, j += 2){

			bReadPairing = CheckPairedAlignmentCandidates(queueItem->EstDistance, queueItem->AlignmentVec1[i], queueItem->AlignmentVec1[j]);

			if (!bReadPairing){
				bReadPairing = RescueUnpairedAlignment(queueItem->EstDistance, queueItem->readArr[i], queueItem->readArr[j], queueItem->AlignmentVec1[i], queueItem->AlignmentVec1[j]);
			} 
			if (bReadPairing){
				RemoveUnMatedAlignmentCandidates(queueItem->AlignmentVec1[i], queueItem->AlignmentVec1[j]);
			} 

			RemoveRedundantCandidates(queueItem->AlignmentVec1[i]); 
			RemoveRedundantCandidates(queueItem->AlignmentVec1[j]);
			// if (bDebugMode)
			// {
			// 	ShowAlignmentCandidateInfo(1, ReadArr[i].header, AlignmentVec1);
			// 	ShowAlignmentCandidateInfo(0, ReadArr[j].header, AlignmentVec2);
			// }
			GenMappingReport(true,  queueItem->readArr[i], queueItem->AlignmentVec1[i]);
			GenMappingReport(false, queueItem->readArr[j], queueItem->AlignmentVec1[j]);

			CheckPairedFinalAlignments(queueItem->readArr[i], queueItem->readArr[j]);

			
		
			SetPairedAlignmentFlag(queueItem->readArr[i], queueItem->readArr[j]);



			EvaluateMAPQ(queueItem->readArr[i]); 
			EvaluateMAPQ(queueItem->readArr[j]);

			//if (bDebugMode) printf("\nEnd of mapping for read#%s\n%s\n", ReadArr[i].header, string().assign(100, '=').c_str());
		}
		
		for (i = 0; i != queueItem->readNum; i += 2) {
			// OutputPairedAlignments(queueItem->readArr[i], queueItem->readArr[i+1], buffer, samItem->myUniqueMapping, samItem->myUnMapping, samItem->myDistance, samItem->myPaired,samItem->SamOutputVec[i],samItem->SamOutputVec[i+1]);
			OutputPairedAlignments_block(queueItem->readArr[i], queueItem->readArr[i+1], buffer, samItem);
		}		

	}else {
        // fprintf(stderr, "子进程%d,接收到第%d份数据的地址为%ld queueItem->readNum:%d\n",threadId,queueItem->index,queueItem,queueItem->readNum);
        for(i=0;i<queueItem->readNum;i++){
            queueItem->EncodeSeq[i] = new uint8_t[queueItem->readArr[i].rlen];
            EnCodeReadSeq(queueItem->readArr[i].rlen, queueItem->readArr[i].seq, queueItem->EncodeSeq[i],&queueItem->readArr[i].gotN);
			// IdentifySeedPairs_FastMode_1(queueItem->readArr[i].rlen, queueItem->EncodeSeq[i] ,queueItem);
        }

		for(i=0;i<queueItem->readNum;i++){

			IdentifySeedPairs_FastMode_Recycle(queueItem->readArr[i].rlen, queueItem->EncodeSeq[i],queueItem->SeedPairVec1[i]);
            GenerateAlignmentCandidateForIlluminaSeq_Recycle(queueItem->readArr[i].rlen, queueItem->SeedPairVec1[i],queueItem->AlignmentVec1[i]);
			// queueItem->SeedPairVec1[i] = IdentifySeedPairs_FastMode(queueItem->readArr[i].rlen, queueItem->EncodeSeq[i]); 

            // queueItem->AlignmentVec1[i] = GenerateAlignmentCandidateForIlluminaSeq(queueItem->readArr[i].rlen, queueItem->SeedPairVec1[i]);
            // queueItem->AlignmentVec1[i] =GenerateAlignmentCandidateForIlluminaSeq_1(queueItem->readArr[i].rlen, queueItem);

            GenMappingReport(true, queueItem->readArr[i], queueItem->AlignmentVec1[i]);
            SetSingleAlignmentFlag(queueItem->readArr[i]);
            // EvaluateMAPQ(queueItem->readArr[i]);
            OutputSingledAlignments_block(queueItem->readArr[i], buffer, samItem);            

		}
		for(i=0;i<queueItem->readNum;i++){
			delete[] queueItem->EncodeSeq[i];
		}
		// memset(&queueItem->EncodeSeq,0,sizeof(uint8_t *)* ReadChunkSize);
		// fprintf(stderr, "子进程%d,处理完的第%d份数据,地址为%ld queueItem->readNum:%d\n",threadId,queueItem->index,queueItem,queueItem->readNum);

        // myUniqueMapping = myUnMapping = 0; 
    }
	for (auto & p : qs) {
		delete p;
	}

	for (i=0;i<queueItem->readNum;i++){
		vector<SeedPair_t>().swap(queueItem->SeedPairVec1[i]);
		vector<AlignmentCandidate_t>().swap(queueItem->AlignmentVec1[i]);
		if(queueItem->readArr[i].CanNum > 0) delete[] queueItem->readArr[i].AlnReportArr;
	}

	delete[] buffer;


	spsc_queue_input.enqueue(queueItem);
	spsc_queue_output.enqueue(samItem); 
	// process_num.fetch_add(-1);



}






void processInit(){


    OutputFileFormat =0;
	OutputByOrder =0;

	if(threadNum <=8){
		ReadChunkSize =2048;
		InputQueueSize=2048;		
	}else if(threadNum >=32){
		ReadChunkSize =8192;
		InputQueueSize=512;
	}else{
		ReadChunkSize =threadNum*256;
		InputQueueSize=256;

	}

	// if(threadNum>=16){
	// 	ReadChunkSize =4096;
	// 	InputQueueSize=1024;
	// }	

	int sam_fd=0;
    if (OutputFileFormat == 0){

		sam_out_fd = open(outputFile, O_WRONLY   | O_TRUNC , 0645);
		close(sam_out_fd);

		sam_out_fd = open(outputFile, O_WRONLY   | O_CREAT , 0645);

	} else{
		bam_out = sam_open_format(outputFile, "wb", NULL);
	} 

	if (sam_out_fd == NULL && bam_out == NULL){
		fprintf(stderr, "Error! Cannot open file [%s]\n", outputFile);
		exit(1);
	}

	writeHeader(sam_out);



	gzReadFileHandler1 = gzReadFileHandler2 = NULL; 
    ReadFileHandler1 = ReadFileHandler2 = NULL;

	if (readFile1.substr(readFile1.find_last_of('.') + 1) == "gz"){
		gzCompressed = true;
	} else{
		gzCompressed = false;
	} 

	FastQFormat = CheckReadFormat(readFile1.c_str());

	if (gzCompressed){
		gzReadFileHandler1 = gzopen(readFile1.c_str(), "rb");

	}else{
		ReadFileHandler1 = fopen(readFile1.c_str(), "r");
	} 


	if (strlen(readFile1.c_str()) == strlen(readFile2.c_str())){
		bSepLibrary = bPairEnd = true;
		if (FastQFormat == CheckReadFormat(readFile2.c_str())){
			if (gzCompressed){
				gzReadFileHandler2 = gzopen(readFile2.c_str(), "rb");
	
			}else{
				ReadFileHandler2 = fopen(readFile2.c_str(), "r");
			
			} 
		}else{
			fprintf(stdout, "Error! %s and %s are with different format...\n", (char*)readFile1.c_str(), (char*)readFile2.c_str());
		}

	}else {
		bSepLibrary = false;
	}


    init_resources(2, &input_info, &output_info);
    input_info.queue_size =InputQueueSize;
    output_info.queue_size =OutputQueueSize;
    output_info.ioFlag =0;
    output_info.runFlag =1;

    input_info.shutdown = output_info.shutdown = 0;

	tmp1 = new char[block_size*2];
	tmp2 = new char[block_size*2];
	tmp1_block_1 =new char[block_size*2];
	tmp1_block_2 =new char[block_size*2];
	// for (MinSeedLength = 13; MinSeedLength < 20; MinSeedLength++) {
	// 	if (TwoGenomeSize < pow(4, MinSeedLength)) {
	// 		break;
	// 	}
	// }

	MinSeedLength=20;
	actual_block_size = block_size - 1;
}

void   processFree(){

    free_resources(2, &input_info, &output_info);

}




void process(){
    int  consumer_count=0,i;


    processInit();


	QueueItem_t *queueItem =new QueueItem_t[InputQueueSize];
	for(int i=0;i<InputQueueSize;i++){
		int   chunkSize =ReadChunkSize;
		queueItem[i].readArr =new ReadItem_t[chunkSize];

		for(int j=0;j<chunkSize;j++){
			(queueItem[i].readArr+j)->header =new char[512];
			(queueItem[i].readArr+j)->seq =new char[512];
			(queueItem[i].readArr+j)->qual =new char[512];
			(queueItem[i].readArr+j)->gotN =0;
		}

		queueItem[i].SamOutputVec =new vector<string>[chunkSize];
		queueItem[i].AlignmentVec1 =new vector<AlignmentCandidate_t>[chunkSize];
		queueItem[i].SeedPairVec1 =new vector<SeedPair_t>[chunkSize];
		queueItem[i].EncodeSeq =new uint8_t*[chunkSize];
		// printf("%ld\n",&queueItem[i]);
		spsc_queue_input.enqueue(&queueItem[i]);
	}


	SamOutputItem_t *samOutputItem =new SamOutputItem_t[OutputQueueSize];
	for(int i=0;i<OutputQueueSize;i++){
		int chunkSize =ReadChunkSize;
		int offset =0;
		samOutputItem[i].sam_block_size =output_block_size;
		samOutputItem[i].samOutputBuffer =(char *)malloc(output_block_size);
		samOutputItem[i].offset=0;
		samOutputItem[i].index =0;
		samOutputItem[i].readNum =0;
		samOutputItem[i].myUniqueMapping=0;
		samOutputItem[i].myUnMapping=0;
		samOutputItem[i].myDistance=0;
		samOutputItem[i].myPaired=0;
		// samOutputItem[i].SamOutputVec =new vector<string>[chunkSize];
		spsc_queue_output_addr.enqueue(&samOutputItem[i]);
	}


	int len;
	char buffer[1024];
	kstring_t str = { 0, 0, NULL };


	StartProcessTime = time(NULL);


	// hashes = new HashPair;
	// int num_hashes =0;
	// hashes = get_hash_function("/home/b8402/22_liangjialang/tools/bloomtree/hashfile", num_hashes);


	//  bf =new UncompressedBF("/home/b8402/22_liangjialang/tools/bloomtree/chr5_3.bf.bv",*hashes,num_hashes);



	root = read_bloom_tree_hash("/home/b8402/22_liangjialang/dataset/chroms_sbt_2/chromSBT.bloomtree","/home/b8402/22_liangjialang/tools/bloomtree/hashfile");


	// bf->load();
	// printf("wdnmdddd:%d\n",bf->size());


	// fclose(rfile_text);
	// precision, recall
	bffile = fopen("sbt_2_40.csv", "w");



	//写线程
	std::thread threadOutput;

	if(OutputByOrder){
		threadOutput =std::thread(process_write);
	}else{
		threadOutput =std::thread(process_write_disorder);
	}
    
	int real_threadNum =threadNum >=2 ? threadNum-2 :1;
    auto *executor = new Executor<QueueItem_t *>([&](QueueItem_t * read_address) {process_mapping(read_address);}, real_threadNum);
	if(bBlockRead){

		process_read_block_aio(executor);
		// process_read_block(executor);		
	}else{
		process_read(executor);

	}

    output_info.shutdown =1;


	std::unique_lock<std::mutex> lock(run_mutex_t);
	// lock.lock();
	while(output_info.runFlag){
		run_cv_t.wait(lock);
	}
	// lock.unlock();

    executor->join();
    delete executor;
	threadOutput.join();
    fprintf(stderr, "线程结束:%d\n",consumer_count );



	for(int i=0;i<InputQueueSize;i++){

		delete[] queueItem[i].readArr;
		// delete[] queueItem[i].SamOutputVec ;
		delete[] queueItem[i].AlignmentVec1;
		delete[] queueItem[i].SeedPairVec1;
		delete[] queueItem[i].EncodeSeq ;
		// printf("%ld\n",&queueItem[i]);
	}

	for(int i=0;i<OutputQueueSize;i++){
		free(samOutputItem[i].samOutputBuffer)  ;

	}	
	delete[] queueItem;
	delete[] samOutputItem;

	// free(block_buff);
    processFree();

	// delete bf;
	fclose(bffile);


	fprintf(stderr, "\rAll the %lld %s reads have been processed in %lld seconds.\n", (long long)iTotalReadNum, (bPairEnd? "paired-end":"single-end"), (long long)(time(NULL) - StartProcessTime));

	if(iTotalReadNum > 0){
		if (bPairEnd) fprintf(stderr, "\t# of total mapped sequences = %lld (sensitivity = %.2f%%)\n\t# of paired sequences = %lld (%.2f%%), average insert size = %d\n", (long long)(iTotalReadNum - iUnMapping), (int)(10000 * (1.0*(iTotalReadNum - iUnMapping) / iTotalReadNum) + 0.5) / 100.0, (long long)iPaired, (int)(10000 * (1.0*iPaired / iTotalReadNum) + 0.5) / 100.0, (iPaired > 1 ? (int)(iDistance / (iPaired >> 1)) : 0));
		else fprintf(stderr, "\t# of total mapped sequences = %lld (sensitivity = %.2f%%)\n", (long long)(iTotalReadNum - iUnMapping), (int)(10000 * (1.0*(iTotalReadNum - iUnMapping) / iTotalReadNum) + 0.5) / 100.0);
		fprintf(stderr, "Alignment output: %s\n", outputFile);
	}

	// exit(0);
    if (OutputFileFormat == 0) {
		// fprintf(stderr, "线程:fclose\n");
		// fsync(fileno(sam_out));
		// fclose(sam_out);

		close(sam_out_fd);
	}else{
		sam_close(bam_out);
	} 

}


