/* The MIT License

   Copyright (c) 2008 Genome Research Ltd (GRL).

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Contact: Heng Li <lh3@sanger.ac.uk> */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "bntseq.h"
#include "bwt.h"
#include "utils.h"

int64_t bwa_seq_len(const char *fn_pac)
{
	FILE *fp;
	int64_t pac_len;
	ubyte_t c;
	fp = xopen(fn_pac, "rb");
	err_fseek(fp, -1, SEEK_END);
	pac_len = err_ftell(fp);
	err_fread_noeof(&c, 1, 1, fp);
	err_fclose(fp);
	return (pac_len - 1) * 4 + (int)c;
}

#define bwt_B00(b, k) ((b)->bwt[(k)>>4]>>((~(k)&0xf)<<1)&3)

void bwt_bwtupdate_core(bwt_t *bwt)
{
	bwtint_t i, k, c[4], n_occ;
	uint32_t *buf;

	n_occ = (bwt->seq_len + OCC_INTERVAL - 1) / OCC_INTERVAL + 1;
	bwt->bwt_size += n_occ * sizeof(bwtint_t); // the new size
	buf = (uint32_t*)calloc(bwt->bwt_size, 4); // will be the new bwt
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = 0; i < bwt->seq_len; ++i) {
		if (i % OCC_INTERVAL == 0) {
			memcpy(buf + k, c, sizeof(bwtint_t) * 4);
			k += sizeof(bwtint_t); // in fact: sizeof(bwtint_t)=4*(sizeof(bwtint_t)/4)
		}
		if (i % 16 == 0) buf[k++] = bwt->bwt[i/16]; // 16 == sizeof(uint32_t)/2
		++c[bwt_B00(bwt, i)];
	}
	// the last element
	memcpy(buf + k, c, sizeof(bwtint_t) * 4);
	xassert(k + sizeof(bwtint_t) == bwt->bwt_size, "inconsistent bwt_size");
	// update bwt
	free(bwt->bwt); bwt->bwt = buf;
}




void gen_tenlen_seq(bwt_t *bwt,const char *fn){
	char * base ="ACTG";
	char **prefix_list ={"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
	uint32_t prefix_no_list[ 16 ] ={0b0000,0b0001,0b0010,0b0011,0b0100,0b0101,0b0110,0b0111,0b1000,0b1001,0b1010,0b1011,0b1100,0b1101,0b1110,0b1111};
	// uint32_t *prefix_no_list ={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
	// uint32_t seq_no_list[1048576]={0};
	// bwtint_t x0[1048576]={0};
	// bwtint_t x1[1048576]={0};
	// bwtint_t x2[1048576]={0};

	uint32_t *seq_no_list =(uint32_t *)malloc(sizeof(uint32_t)*1048576);
	bwtint_t *x0=(bwtint_t *)malloc(sizeof(bwtint_t)*1048576);
	bwtint_t *x1=(bwtint_t *)malloc(sizeof(bwtint_t)*1048576);
	bwtint_t *x2=(bwtint_t *)malloc(sizeof(bwtint_t)*1048576);


	memset(x0,0,sizeof(bwtint_t)*1048576);
	memset(x1,0,sizeof(bwtint_t)*1048576);
	memset(x2,0,sizeof(bwtint_t)*1048576);

	printf("%s\n",fn);
	FILE *fp;
	fp = xopen(fn, "wb");

	// for(int x=0;x<1048576;x++){
	// 	// printf("i:%d\t",x);
	// 	// printf("empty:%d\n",seq_no_list[x]);
	// 	printf("i:%d k:%ld l:%ld s:%ld\n",seq_no_list[x],x0[x],x1[x],x2[x]);


	// }	

	uint32_t seq_no;

	for(int i=0;i<16;i++){
		// seq_no =(uint32_t)i;
		int b;

		bwtintv_t ik, ok[4];
		bwtint_t tk[4], tl[4];

		seq_no = prefix_no_list[i];
		
		uint32_t p =(prefix_no_list[i] >>2)  ;


		ik.x[0] = bwt->L2[p] + 1;
		ik.x[1] = bwt->L2[3 - p] + 1;
		//interval
		ik.x[2] = bwt->L2[p + 1] - bwt->L2[p];


		//正向匹配
		bwt_2occ4(bwt, ik.x[1] - 1, ik.x[1] - 1 + ik.x[2], tk, tl);
		for ( b = 0; b != 4; ++b) {
			ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
			ok[b].x[2] = tl[b] - tk[b];
		}
		ok[3].x[0] = ik.x[0] + (ik.x[1] <= bwt->primary && ik.x[1] + ik.x[2] - 1 >= bwt->primary);
		ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
		ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
		ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
		b = 3 -  (prefix_no_list[i] & 0b0011);

		// printf("%d %d ",prefix_no_list[i] >>2,prefix_no_list[i] & 0b0011);

		//interval
		if (ok[b].x[2] == 0){
			seq_no =(seq_no >>4);
			continue; // extension ends
		} else{
			ik = ok[b];
		} 



		for(int j=0;j<16;j++){
			seq_no =(seq_no <<4) +prefix_no_list[j];

			bwtintv_t ik_34 =ik;
			//正向匹配
			bwt_2occ4(bwt, ik_34.x[1] - 1, ik_34.x[1] - 1 + ik_34.x[2], tk, tl);
			for ( b = 0; b != 4; ++b) {
				ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
				ok[b].x[2] = tl[b] - tk[b];
			}
			ok[3].x[0] = ik_34.x[0] + (ik_34.x[1] <= bwt->primary && ik_34.x[1] + ik_34.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			b = 3 - (prefix_no_list[j] >>2 );
			//interval
			if (ok[b].x[2] == 0){
				seq_no =(seq_no >>4);
				continue; // extension ends
			} else{
				ik_34 = ok[b];
			} 

			bwt_2occ4(bwt, ik_34.x[1] - 1, ik_34.x[1] - 1 + ik_34.x[2], tk, tl);
			for ( b = 0; b != 4; ++b) {
				ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
				ok[b].x[2] = tl[b] - tk[b];
			}
			ok[3].x[0] = ik_34.x[0] + (ik_34.x[1] <= bwt->primary && ik_34.x[1] + ik_34.x[2] - 1 >= bwt->primary);
			ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
			ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
			ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
			b = 3 - (prefix_no_list[j] & 0b0011);
			//interval
			if (ok[b].x[2] == 0){
				seq_no =(seq_no >>4);
				continue; // extension ends
			} else{
				ik_34 = ok[b];
			} 
			// printf("%d %d ",prefix_no_list[j] >>2,(prefix_no_list[j] & 0b0011));

			// fprintf(stdout,"k:%ld l:%d s:%ld\n",ik_34.x[0],ik_34.x[1],ik_34.x[2]);
			
			for(int k=0;k<16;k++){
				seq_no =(seq_no <<4) +prefix_no_list[k];

				bwtintv_t ik_56 =ik_34;
				
				//正向匹配
				bwt_2occ4(bwt, ik_56.x[1] - 1, ik_56.x[1] - 1 + ik_56.x[2], tk, tl);
				for ( b = 0; b != 4; ++b) {
					ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
					ok[b].x[2] = tl[b] - tk[b];
				}
				ok[3].x[0] = ik_56.x[0] + (ik_56.x[1] <= bwt->primary && ik_56.x[1] + ik_56.x[2] - 1 >= bwt->primary);
				ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
				ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
				ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
				b = 3 - (prefix_no_list[k] >>2);
				//interval
				if (ok[b].x[2] == 0){
					seq_no =(seq_no >>4);
					continue; // extension ends
				} else{
					ik_56 = ok[b];
				} 

				bwt_2occ4(bwt, ik_56.x[1] - 1, ik_56.x[1] - 1 + ik_56.x[2], tk, tl);
				for ( b = 0; b != 4; ++b) {
					ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
					ok[b].x[2] = tl[b] - tk[b];
				}
				ok[3].x[0] = ik_56.x[0] + (ik_56.x[1] <= bwt->primary && ik_56.x[1] + ik_56.x[2] - 1 >= bwt->primary);
				ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
				ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
				ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
				b = 3 - (prefix_no_list[k] & 0b0011);
				//interval
				if (ok[b].x[2] == 0){
					seq_no =(seq_no >>4);
					continue; // extension ends
				} else{
					ik_56 = ok[b];
				} 

				// printf("%d %d ",prefix_no_list[i] >>2,prefix_no_list[i] & 0b0011);
				// printf("%d %d ",prefix_no_list[j] >>2,(prefix_no_list[j] & 0b0011));

				// printf("%d %d ",prefix_no_list[k] >>2,(prefix_no_list[k] & 0b0011));
				// printf("k:%ld l:%d s:%ld\n",ik_56.x[0],ik_56.x[1],ik_56.x[2]);
				for(int l=0;l<16;l++){
					seq_no =(seq_no <<4) +prefix_no_list[l];
					bwtintv_t ik_78 =ik_56;

					//正向匹配
					bwt_2occ4(bwt, ik_78.x[1] - 1, ik_78.x[1] - 1 + ik_78.x[2], tk, tl);
					for ( b = 0; b != 4; ++b) {
						ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
						ok[b].x[2] = tl[b] - tk[b];
					}
					ok[3].x[0] = ik_78.x[0] + (ik_78.x[1] <= bwt->primary && ik_78.x[1] + ik_78.x[2] - 1 >= bwt->primary);
					ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
					ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
					ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
					b = 3 - (prefix_no_list[l] >>2);
					//interval
					if (ok[b].x[2] == 0){
						seq_no =(seq_no >>4);
						continue; // extension ends
					} else{
						ik_78 = ok[b];
					} 

					bwt_2occ4(bwt, ik_78.x[1] - 1, ik_78.x[1] - 1 + ik_78.x[2], tk, tl);
					for ( b = 0; b != 4; ++b) {
						ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
						ok[b].x[2] = tl[b] - tk[b];
					}
					ok[3].x[0] = ik_78.x[0] + (ik_78.x[1] <= bwt->primary && ik_78.x[1] + ik_78.x[2] - 1 >= bwt->primary);
					ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
					ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
					ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
					b = 3 - (prefix_no_list[l] & 0b0011);
					//interval
					if (ok[b].x[2] == 0){
						seq_no =(seq_no >>4);
						continue; // extension ends
					} else{
						ik_78 = ok[b];
					} 
					// printf("%d %d ",prefix_no_list[i] >>2,prefix_no_list[i] & 0b0011);
					// printf("%d %d ",prefix_no_list[j] >>2,(prefix_no_list[j] & 0b0011));
					// printf("%d %d ",prefix_no_list[k] >>2,(prefix_no_list[k] & 0b0011));
					// printf("%d %d ",prefix_no_list[l] >>2,(prefix_no_list[l] & 0b0011));
					// printf("k:%ld l:%d s:%ld\n",ik_78.x[0],ik_78.x[1],ik_78.x[2]);


					for(int m=0;m<16;m++){
						seq_no =(seq_no <<4) +prefix_no_list[m];
						bwtintv_t ik_9x =ik_78;

						//正向匹配
						bwt_2occ4(bwt, ik_9x.x[1] - 1, ik_9x.x[1] - 1 + ik_9x.x[2], tk, tl);
						for ( b = 0; b != 4; ++b) {
							ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
							ok[b].x[2] = tl[b] - tk[b];
						}
						ok[3].x[0] = ik_9x.x[0] + (ik_9x.x[1] <= bwt->primary && ik_9x.x[1] + ik_9x.x[2] - 1 >= bwt->primary);
						ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
						ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
						ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
						b = 3 -(prefix_no_list[m] >>2) ;
						//interval
						if (ok[b].x[2] == 0){
							seq_no =(seq_no >>4);
							continue; // extension ends
						} else{
							ik_9x = ok[b];
						} 

						bwt_2occ4(bwt, ik_9x.x[1] - 1, ik_9x.x[1] - 1 + ik_9x.x[2], tk, tl);
						for ( b = 0; b != 4; ++b) {
							ok[b].x[1] = bwt->L2[b] + 1 + tk[b];
							ok[b].x[2] = tl[b] - tk[b];
						}
						ok[3].x[0] = ik_9x.x[0] + (ik_9x.x[1] <= bwt->primary && ik_9x.x[1] + ik_9x.x[2] - 1 >= bwt->primary);
						ok[2].x[0] = ok[3].x[0] + ok[3].x[2];
						ok[1].x[0] = ok[2].x[0] + ok[2].x[2];
						ok[0].x[0] = ok[1].x[0] + ok[1].x[2];
						b = 3 - (prefix_no_list[m] & 0b0011);
						//interval
						if (ok[b].x[2] == 0){
							seq_no =(seq_no >>4);
							continue; // extension ends
						} else{
							ik_9x = ok[b];
						} 

									// fprintf(stdout, "%d\n",seq_no);
						x0[seq_no]=ik_9x.x[0];
						x1[seq_no]=ik_9x.x[1];
						x2[seq_no]=ik_9x.x[2];
						seq_no_list[seq_no]=1;
						// printf("%d %d ",prefix_no_list[i] >>2,prefix_no_list[i] & 0b0011);
						// printf("%d %d ",prefix_no_list[j] >>2,(prefix_no_list[j] & 0b0011));
						// printf("%d %d ",prefix_no_list[k] >>2,(prefix_no_list[k] & 0b0011));
						// printf("%d %d ",prefix_no_list[l] >>2,(prefix_no_list[l] & 0b0011));
						// printf("%d %d ",prefix_no_list[m] >>2,(prefix_no_list[m] & 0b0011));
						// printf("i:%d k:%ld l:%ld s:%ld\n",seq_no,x0[seq_no],x1[seq_no],x2[seq_no]);
						seq_no =(seq_no >>4);
					}
					seq_no =(seq_no >>4);

				}
				seq_no =(seq_no >>4);

			}
			seq_no =(seq_no >>4);

		}
		// seq_no =(seq_no >>4);

	}


	char * line =(char *)malloc(100);

	sprintf(line, "%d\t%ld\n", 10,1048576);
	err_fwrite(line, sizeof(char), strlen(line), fp);
	// printf("%d %s",strlen(line), line);

	// // printf("wdnmddddd\n");
	for(int z=0;z<1048576;z++){
		
		// printf("i:%d\t",x);
		// printf("empty:%d\n",seq_no_list[x]);
		// err_fwrite(z, sizeof(int), 1, fp);
		// err_fwrite(bwt->sa + 1, sizeof(bwtint_t), bwt->n_sa - 1, fp);
		sprintf(line, "%d\t%ld\t%ld\t%ld\n", z,x0[z],x1[z],x2[z]);
		err_fwrite(line, sizeof(char), strlen(line), fp);
	
		// printf("%d %s",strlen(line), line);

		// memset(line,0,1024);
		// printf("z:%d i:%d k:%ld l:%ld s:%ld\n",z ,seq_no_list[z],x0[z],x1[z],x2[z]);

	}	

	err_fflush(fp);
	err_fclose(fp);
	free(x0);	
	free(x1);	
	free(x2);	

}



int bwa_idx_build(const char *fa, const char *prefix) //**//
{
	//int algo_type = 0, block_size = 10000000;
	extern void bwa_pac_rev_core(const char *fn, const char *fn_rev);

	char *str, *str2, *str3,*str4;
	clock_t t;
	int64_t l_pac;

	str  = (char*)calloc(strlen(prefix) + 10, 1);
	str2 = (char*)calloc(strlen(prefix) + 10, 1);
	str3 = (char*)calloc(strlen(prefix) + 10, 1);
	str4 = (char*)calloc(strlen(prefix) + 10, 1);

	// { // nucleotide indexing
	// 	gzFile fp = xzopen(fa, "r");
	// 	t = clock();
	// 	fprintf(stdout, "[bwt_index] Pack FASTA... ");
	// 	l_pac = bns_fasta2bntseq(fp, prefix, 0);
	// 	fprintf(stdout, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// 	err_gzclose(fp);
	// }
	// //if (algo_type == 0) algo_type = l_pac > 50000000? 2 : 3; // set the algorithm for generating BWT
	// {
	// 	//printf("algo_type = %d\n", algo_type);

	// 	strcpy(str, prefix); strcat(str, ".pac");
	// 	strcpy(str2, prefix); strcat(str2, ".bwt");
	// 	t = clock();
	// 	fprintf(stdout, "[bwt_index] Construct BWT for the packed sequence...\n");
	// 	//if (algo_type == 2) 
	// 		bwt_bwtgen2(str, str2, 10000000);
	// 	//else if (algo_type == 1 || algo_type == 3) {
	// 		//bwt_t *bwt;
	// 		//bwt = bwt_pac2bwt(str/*, algo_type == 3*/);
	// 		//bwt_dump_bwt(str2, bwt);
	// 		//bwt_destroy(bwt);
	// 	//}
	// 	fprintf(stdout, "[bwt_index] %.2f seconds elapse.\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// }
	// {
	// 	bwt_t *bwt;
	// 	strcpy(str, prefix); strcat(str, ".bwt");
	// 	t = clock();
	// 	fprintf(stdout, "[bwt_index] Update BWT... ");
	// 	bwt = bwt_restore_bwt(str);
	// 	bwt_bwtupdate_core(bwt);
	// 	bwt_dump_bwt(str, bwt);
	// 	bwt_destroy(bwt);
	// 	fprintf(stdout, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// }
	// {
	// 	gzFile fp = xzopen(fa, "r");
	// 	t = clock();
	// 	fprintf(stdout, "[bwt_index] Pack forward-only FASTA... ");
	// 	l_pac = bns_fasta2bntseq(fp, prefix, 1);
	// 	fprintf(stdout, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	// 	err_gzclose(fp);
	// }
	{
		bwt_t *bwt;
		strcpy(str, prefix); strcat(str, ".bwt");
		strcpy(str3, prefix); strcat(str3, ".sa");
		strcpy(str4, prefix); strcat(str4, ".hsh");

		t = clock();
		fprintf(stdout, "[bwt_index] Construct SA from BWT and Occ... ");
		bwt = bwt_restore_bwt(str);

		fprintf(stdout, "\n[bwt_index] Construct Hash from BWT and Occ... ");
		gen_tenlen_seq(bwt,str4);
		// bwt_cal_sa(bwt, 32);
		// bwt_dump_sa(str3, bwt);
		bwt_destroy(bwt);
		fprintf(stdout, "%.2f sec\n", (float)(clock() - t) / CLOCKS_PER_SEC);
	}
	free(str3); free(str2); free(str);
	return 0;
}
