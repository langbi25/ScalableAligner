#include "util.h"

bool CompByPosDiff(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.PosDiff == p2.PosDiff) return p1.rPos < p2.rPos;
	else return p1.PosDiff < p2.PosDiff;
}




vector<SeedPair_t> IdentifySeedPairs_FastMode(int rlen, uint8_t* EncodeSeq)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos ,seed_len;
	vector<SeedPair_t> SeedPairVec;
	vector<bwtint_t> il_list;
	vector<int> interval_list;
	vector<int> match_len_list;
	vector<int> match_beg_list;		
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; 
	
	// pos = 0, end_pos = rlen - MinSeedLength;
	seed_len =20;	
	pos = 0, end_pos = rlen - seed_len;

	int last_rightest =0;


	while (pos < end_pos){
		if (EncodeSeq[pos] > 3) {
			pos++;
		}else{
			bwtSearchResult = BWT_Search_Forward_3(EncodeSeq, pos, rlen,last_rightest,&il_list, &interval_list,&match_len_list,&match_beg_list);
			last_rightest =bwtSearchResult.rightset;
			if(last_rightest ==rlen){
				break;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}



	if(bwtSearchResult.rightset!=rlen){
		bwtSearchResult =BWT_Search_Backward(EncodeSeq, rlen, 0,&il_list, &interval_list,&match_len_list,&match_beg_list);
		// fprintf(stderr,"%d-%d interval:%d\n",rlen,bwtSearchResult.len ,bwtSearchResult.freq);
	}

	for(int i=0;i<il_list.size();i++){	
		// fprintf(stderr,"%d-%d-%d-%ld\n",match_beg_list[i],match_len_list[i]+match_beg_list[i],interval_list[i],il_list[i]);
		bwt_sa_batch(il_list[i],interval_list[i],match_len_list[i],match_beg_list[i],SeedPairVec);

	}

	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByPosDiff);
	return SeedPairVec;
}



void IdentifySeedPairs_FastMode_Recycle(int rlen, uint8_t* EncodeSeq ,vector<SeedPair_t>& SeedPairVec)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos ,seed_len;
	// vector<SeedPair_t> SeedPairVec;
	vector<bwtint_t> il_list;
	vector<int> interval_list;
	vector<int> match_len_list;
	vector<int> match_beg_list;		
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; 
	
	// pos = 0, end_pos = rlen - MinSeedLength;
	seed_len =20;	
	pos = 0, end_pos = rlen - seed_len;

	int last_rightest =0;
	// for(int i=0;i<rlen;i++){
	// 	printf("%d",EncodeSeq[i]);
	// }
	// printf("\n");

	while (pos < end_pos){
		if (EncodeSeq[pos] > 3) {
			pos++;
		}else{
			if(bHash){
				bwtSearchResult = BWT_Search_Forward_hash(EncodeSeq, pos, rlen,last_rightest,&il_list, &interval_list,&match_len_list,&match_beg_list);
			}else{
				bwtSearchResult = BWT_Search_Forward_3(EncodeSeq, pos, rlen,last_rightest,&il_list, &interval_list,&match_len_list,&match_beg_list);
			}
			
			last_rightest =bwtSearchResult.rightset;
			if(last_rightest ==rlen){
				break;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}



	if(bwtSearchResult.rightset!=rlen){
		if(bHash){
			bwtSearchResult =BWT_Search_Backward_hash(EncodeSeq, rlen, 0,&il_list, &interval_list,&match_len_list,&match_beg_list);
		}else{

			bwtSearchResult =BWT_Search_Backward(EncodeSeq, rlen, 0,&il_list, &interval_list,&match_len_list,&match_beg_list);
		}
		// fprintf(stderr,"%d-%d interval:%d\n",rlen,bwtSearchResult.len ,bwtSearchResult.freq);
	}

	for(int i=0;i<il_list.size();i++){	
		bwt_sa_batch(il_list[i],interval_list[i],match_len_list[i],match_beg_list[i],SeedPairVec);
		// fprintf(stderr,"%d-%d-%d-%ld\n",match_beg_list[i],match_len_list[i]+match_beg_list[i],interval_list[i],il_list[i]);

	}

	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByPosDiff);
	// return SeedPairVec;
}

void IdentifySeedPairs_FastMode_getN(int rlen, uint8_t* EncodeSeq ,vector<SeedPair_t>& SeedPairVec ,int gotN)
{
	SeedPair_t SeedPair;
	int i, pos, end_pos ,seed_len;
	// vector<SeedPair_t> SeedPairVec;
	vector<bwtint_t> il_list;
	vector<int> interval_list;
	vector<int> match_len_list;
	vector<int> match_beg_list;		
	bwtSearchResult_t bwtSearchResult;

	SeedPair.bSimple = true; 
	
	// pos = 0, end_pos = rlen - MinSeedLength;
	seed_len =20;	
	pos = 0, end_pos = rlen - seed_len;

	int last_rightest =0;
	// for(int i=0;i<rlen;i++){
	// 	printf("%d",EncodeSeq[i]);
	// }
	// printf("\n");
	// printf("gotN:%d\n",gotN);
	while (pos < end_pos){
		if (EncodeSeq[pos] > 3) {
			pos++;
		}else{
			if(bHash && !gotN){
				bwtSearchResult = BWT_Search_Forward_hash(EncodeSeq, pos, rlen,last_rightest,&il_list, &interval_list,&match_len_list,&match_beg_list);
			}else{
				bwtSearchResult = BWT_Search_Forward_3(EncodeSeq, pos, rlen,last_rightest,&il_list, &interval_list,&match_len_list,&match_beg_list);
			}
			
			last_rightest =bwtSearchResult.rightset;
			if(last_rightest ==rlen){
				break;
			}
			pos += (bwtSearchResult.len + 1);
		}
	}



	if(bwtSearchResult.rightset!=rlen){
		if(bHash && !gotN){
			bwtSearchResult =BWT_Search_Backward_hash(EncodeSeq, rlen, 0,&il_list, &interval_list,&match_len_list,&match_beg_list);
		}else{

			bwtSearchResult =BWT_Search_Backward(EncodeSeq, rlen, 0,&il_list, &interval_list,&match_len_list,&match_beg_list);
		}
		// fprintf(stderr,"%d-%d interval:%d\n",rlen,bwtSearchResult.len ,bwtSearchResult.freq);
	}

	for(int i=0;i<il_list.size();i++){	
		bwt_sa_batch(il_list[i],interval_list[i],match_len_list[i],match_beg_list[i],SeedPairVec);
		// fprintf(stderr,"%d-%d-%d-%ld\n",match_beg_list[i],match_len_list[i]+match_beg_list[i],interval_list[i],il_list[i]);

	}

	sort(SeedPairVec.begin(), SeedPairVec.end(), CompByPosDiff);
	// return SeedPairVec;
}




bool CheckCoordinateValidity(vector<SeedPair_t>& SeedVec)
{
	bool bValid = true;
	int64_t gPos1=0, gPos2=TwoGenomeSize;
	map<int64_t, int>::iterator iter1, iter2;

	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos1 = iter->gPos;
			break;
		}
	}
	for (vector<SeedPair_t>::reverse_iterator iter = SeedVec.rbegin(); iter != SeedVec.rend(); iter++)
	{
		if (iter->gLen > 0)
		{
			gPos2 = iter->gPos + iter->gLen - 1;
			break;
		}
	}	
	if ((gPos1 < GenomeSize && gPos2 >= GenomeSize) || (gPos1 >= GenomeSize && gPos2 < GenomeSize) || ((iter1 = ChrLocMap.lower_bound(gPos1)) == ChrLocMap.end() || (iter2 = ChrLocMap.lower_bound(gPos2)) == ChrLocMap.end() || iter1->second != iter2->second))
	{
		bValid = false;
		//if (bDebugMode) fprintf(stdout, "%lld and %lld are not in the same chromosome!\n", gPos1, gPos2);
	}
	return bValid;
}
bool CompByReadPos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	return p1.rPos < p2.rPos;
}
vector<AlignmentCandidate_t> GenerateAlignmentCandidateForIlluminaSeq(int rlen, vector<SeedPair_t> SeedPairVec)
{
	int64_t gPos_end;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;
	vector<AlignmentCandidate_t> AlignmentVec;

	vector<pair<int,int>> rpos_list;
	vector<pair<int64_t,int64_t>> gpos_list;
	int rpos_num=0;
	if ((thr = (int)(rlen*0.2)) > 20) thr = 20;

	//if (bDebugMode) printf("\n\nRaw seeds:\n"), ShowSeedInfo(SeedPairVec);
	AlignmentCandidate.PairedAlnCanIdx = -1; num = (int)SeedPairVec.size();

	i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
	for (; i < num;)
	{
		AlignmentCandidate.Score = SeedPairVec[i].rLen; 
		AlignmentCandidate.PosDiff =SeedPairVec[i].PosDiff;
		AlignmentCandidate.PosDiff_End =SeedPairVec[i].PosDiff +rlen;
		gPos_end = GetAlignmentBoundary(SeedPairVec[i].gPos); 
		AlignmentCandidate.SeedVec.clear();
		rpos_list.clear();
		gpos_list.clear();
		//AlignmentCandidate.SeedVec.resize(1); AlignmentCandidate.SeedVec[0] = SeedPairVec[i]; 
		//if (bDebugMode) 		printf("\nMaster seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen, SeedPairVec[i].PosDiff);
		rpos_list.push_back(make_pair(SeedPairVec[i].rPos ,SeedPairVec[i].rPos+SeedPairVec[i].rLen-1));
		gpos_list.push_back(make_pair(SeedPairVec[i].gPos ,SeedPairVec[i].gPos+SeedPairVec[i].gLen-1));
		for (j = i, k = i + 1; k < num; k++){
			if (SeedPairVec[k].gPos > gPos_end || (SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff) > MaxGaps) break;
			else
			{	
				rpos_num =gpos_list.size();
				for (int q = 0; q != rpos_num; q++){
					bool if_read_overlap =SeedPairVec[k].rPos >=rpos_list[q].first && SeedPairVec[k].rPos <=rpos_list[q].second;
					bool if_ref_overlap =SeedPairVec[k].gPos >=gpos_list[q].first && SeedPairVec[k].gPos <=gpos_list[q].second;
					if(if_read_overlap || if_ref_overlap){
						// printf("%d-%d no way add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", if_read_overlap,if_ref_overlap,SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
						int indel =(rpos_list[q].second)-SeedPairVec[k].rPos;
						SeedPairVec[k].rPos +=(indel+1);
						SeedPairVec[k].gPos +=(indel+1);
						SeedPairVec[k].gLen -=(indel+1);
						SeedPairVec[k].rLen -=(indel+1);
						break;
					}
				}
				//if (bDebugMode) 				printf("add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
				AlignmentCandidate.PosDiff_End =SeedPairVec[k].PosDiff +rlen;
				AlignmentCandidate.Score += SeedPairVec[k].rLen;
				rpos_list.push_back(make_pair(SeedPairVec[k].rPos ,SeedPairVec[k].rPos+SeedPairVec[k].rLen-1));
				gpos_list.push_back(make_pair(SeedPairVec[k].gPos ,SeedPairVec[k].gPos+SeedPairVec[k].gLen-1));
				//AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
				j = k;
			}
		}
		if (AlignmentCandidate.Score > thr){
			copy(SeedPairVec.begin() + i, SeedPairVec.begin() + k, back_inserter(AlignmentCandidate.SeedVec));
			if (AlignmentCandidate.Score - 50 > thr) thr = AlignmentCandidate.Score - 50;
			// AlignmentCandidate.PosDiff = AlignmentCandidate.SeedVec[0].PosDiff;
			if (AlignmentCandidate.PosDiff < 0) AlignmentCandidate.PosDiff = 0;
			sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByReadPos);

			AlignmentVec.push_back(AlignmentCandidate);
			//if (bDebugMode)
			//{
				// printf("Candidate score = %d\n", AlignmentCandidate.Score);
			//	ShowSeedLocationInfo(AlignmentCandidate.SeedVec[0]);
			//	ShowSeedInfo(AlignmentCandidate.SeedVec);
			//}
		}
		i = k;

	}

	int ali_thr;
	vector<AlignmentCandidate_t>::iterator iter;

	if (AlignmentVec.size() <= 1) return AlignmentVec;
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
		if (bPacBioData || score1 == score2 || score1 - score2 > 20) ali_thr = score1;
		else ali_thr = score2;

		//if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
		// for(int x=0;x<AlignmentVec.size();x++){
		// 	fprintf(stderr,"score:%d\n",AlignmentVec[x].Score);
		// }
		for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++){
			if (iter->Score < ali_thr){
				iter->Score = 0;
			} else{

				int i, j, rGaps, gGaps, num ,glen =-1;
				SeedPair_t SeedPair;
				if (iter->SeedVec.size() > 1)
				{	
					// printf("new seedvec");
					// for(int x=0;x <iter->SeedVec.size();x++){
					// 	printf("%d\t",iter->SeedVec[x].rPos);
					// }
					// printf("\n");


					SeedPair.bSimple = false;
					num = (int)iter->SeedVec.size();

					for (i = 0, j = 1; j < num; i++, j++){
						rGaps = iter->SeedVec[j].rPos - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen); 
						if (rGaps < 0) rGaps = 0;
						gGaps = iter->SeedVec[j].gPos - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen); 
						if (gGaps < 0) gGaps = 0;
						//printf("check %d and %d: rGaps=%d, gGaps=%d\n", i + 1, j + 1, rGaps, gGaps); fflush(stdout);
						if (rGaps > 0 || gGaps > 0){
							SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
							SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
							SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
							SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
							//printf("Add normal pair:\nR1[%d-%d] G1[%lld-%lld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%lld-%lld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
							iter->SeedVec.push_back(SeedPair);
						}
					}
					if ((int)iter->SeedVec.size() > num) inplace_merge(iter->SeedVec.begin(), iter->SeedVec.begin()+num, iter->SeedVec.end(), CompByGenomePos);
					//if(bDebugMode) printf("After filling gaps\n"), ShowSeedInfo(SeedVec);
				}
				// Check missing blocks at both ends
				if (iter->SeedVec.size() > 0){
					i = 0;
					//printf("Identify heading pair R[%d-%d]=%d, G[%lld-%lld]=%d\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].rLen, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].gLen); fflush(stdout);
					rGaps = iter->SeedVec[i].rPos > 0 ? iter->SeedVec[i].rPos : 0;
					gGaps = glen > 0 ? iter->SeedVec[i].gPos : rGaps;
					//printf("rGaps=%d, gGaps=%d\n", rGaps, gGaps); fflush(stdout);
					if (rGaps > 0 || gGaps > 0){
						SeedPair.rPos = 0;
						SeedPair.gPos = iter->SeedVec[i].gPos - gGaps;
						if (SeedPair.gPos < 0) { SeedPair.gPos = 0; gGaps += SeedPair.gPos; }
						SeedPair.PosDiff = SeedPair.gPos;
						SeedPair.bSimple = false;
						SeedPair.rLen = rGaps;
						SeedPair.gLen = gGaps;
						iter->SeedVec.insert(iter->SeedVec.begin(), SeedPair);
						//if (bDebugMode) printf("Add missing head: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
					}
					i = (int)iter->SeedVec.size() - 1;
					rGaps = rlen - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen);
					gGaps = glen > 0 ? glen - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen) : rGaps;

					if (rGaps > 0 || gGaps > 0){
						SeedPair.bSimple = false;
						SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
						SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
						SeedPair.rLen = rGaps;
						SeedPair.gLen = gGaps;
						iter->SeedVec.push_back(SeedPair);
						//if (bDebugMode) printf("Add missing tail: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
					}
				}

			}
		}
	}

	return AlignmentVec;
}




void GenerateAlignmentCandidateForIlluminaSeq_Recycle(int rlen, vector<SeedPair_t>& SeedPairVec,vector<AlignmentCandidate_t>& AlignmentVec)
{
	int64_t gPos_end;
	int i, j, k, thr, num;
	AlignmentCandidate_t AlignmentCandidate;


	vector<pair<int,int>> rpos_list;
	vector<pair<int64_t,int64_t>> gpos_list;
	int rpos_num=0;
	if ((thr = (int)(rlen*0.2)) > 20) thr = 20;

	//if (bDebugMode) printf("\n\nRaw seeds:\n"), ShowSeedInfo(SeedPairVec);
	AlignmentCandidate.PairedAlnCanIdx = -1; num = (int)SeedPairVec.size();

	i = 0; while (i < num && SeedPairVec[i].PosDiff < 0) i++;
	for (; i < num;){
		AlignmentCandidate.Score = SeedPairVec[i].rLen; 
		AlignmentCandidate.PosDiff =SeedPairVec[i].PosDiff;
		AlignmentCandidate.PosDiff_End =SeedPairVec[i].PosDiff +rlen;
		gPos_end = GetAlignmentBoundary(SeedPairVec[i].gPos); 
		AlignmentCandidate.SeedVec.clear();
		rpos_list.clear();
		gpos_list.clear();
		//AlignmentCandidate.SeedVec.resize(1); AlignmentCandidate.SeedVec[0] = SeedPairVec[i]; 
		//if (bDebugMode) 				printf("\nMaster seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[i].rPos, SeedPairVec[i].rPos + SeedPairVec[i].rLen - 1, SeedPairVec[i].gPos, SeedPairVec[i].gPos + SeedPairVec[i].gLen - 1, SeedPairVec[i].rLen, SeedPairVec[i].PosDiff);
		rpos_list.push_back(make_pair(SeedPairVec[i].rPos ,SeedPairVec[i].rPos+SeedPairVec[i].rLen-1));
		gpos_list.push_back(make_pair(SeedPairVec[i].gPos ,SeedPairVec[i].gPos+SeedPairVec[i].gLen-1));
		for (j = i, k = i + 1; k < num; k++){
			if (SeedPairVec[k].gPos > gPos_end || (SeedPairVec[k].PosDiff - SeedPairVec[j].PosDiff) > MaxGaps) break;
			else
			{	
				rpos_num =gpos_list.size();
				for (int q = 0; q != rpos_num; q++){
					bool if_read_overlap =SeedPairVec[k].rPos >=rpos_list[q].first && SeedPairVec[k].rPos <=rpos_list[q].second;
					bool if_ref_overlap =SeedPairVec[k].gPos >=gpos_list[q].first && SeedPairVec[k].gPos <=gpos_list[q].second;
					if(if_read_overlap || if_ref_overlap){
						// printf("%d-%d no way add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", if_read_overlap,if_ref_overlap,SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
						int indel =(rpos_list[q].second)-SeedPairVec[k].rPos;
						SeedPairVec[k].rPos +=(indel+1);
						SeedPairVec[k].gPos +=(indel+1);
						SeedPairVec[k].gLen -=(indel+1);
						SeedPairVec[k].rLen -=(indel+1);
						break;
					}
				}
				//if (bDebugMode) 				printf("add seed: r[%d-%d] g[%lld-%lld], len=%d, PosDiff=%lld\n", SeedPairVec[k].rPos, SeedPairVec[k].rPos + SeedPairVec[k].rLen - 1, SeedPairVec[k].gPos, SeedPairVec[k].gPos + SeedPairVec[k].gLen - 1, SeedPairVec[k].rLen, SeedPairVec[k].PosDiff);
				AlignmentCandidate.PosDiff_End =SeedPairVec[k].PosDiff +rlen;
				AlignmentCandidate.Score += SeedPairVec[k].rLen;
				rpos_list.push_back(make_pair(SeedPairVec[k].rPos ,SeedPairVec[k].rPos+SeedPairVec[k].rLen-1));
				gpos_list.push_back(make_pair(SeedPairVec[k].gPos ,SeedPairVec[k].gPos+SeedPairVec[k].gLen-1));
				//AlignmentCandidate.SeedVec.push_back(SeedPairVec[k]);
				j = k;
			}
		}
		if (AlignmentCandidate.Score > thr){
			copy(SeedPairVec.begin() + i, SeedPairVec.begin() + k, back_inserter(AlignmentCandidate.SeedVec));
			if (AlignmentCandidate.Score - 50 > thr) thr = AlignmentCandidate.Score - 50;
			// AlignmentCandidate.PosDiff = AlignmentCandidate.SeedVec[0].PosDiff;
			if (AlignmentCandidate.PosDiff < 0) AlignmentCandidate.PosDiff = 0;
			sort(AlignmentCandidate.SeedVec.begin(), AlignmentCandidate.SeedVec.end(), CompByReadPos);

			AlignmentVec.push_back(AlignmentCandidate);
			//if (bDebugMode)
			//{
				// printf("Candidate score = %d\n", AlignmentCandidate.Score);
			//	ShowSeedLocationInfo(AlignmentCandidate.SeedVec[0]);
			//	ShowSeedInfo(AlignmentCandidate.SeedVec);
			//}
		}
		i = k;

	}

	// int ali_thr;
	// vector<AlignmentCandidate_t>::iterator iter;


	// // for(int x =0;x<AlignmentVec.size();x++){
	// // 	printf("Score:%d SeedVec.size():%d\n",AlignmentVec[x].Score,AlignmentVec[x].SeedVec.size());
	// // 	for(int y=0;y<AlignmentVec[x].SeedVec.size();y++){
	// // 		printf("bSimple:%d rLen:%d\n",AlignmentVec[x].SeedVec[y].bSimple,AlignmentVec[x].SeedVec[y].rLen);
	// // 	}
	// // }



	// if (AlignmentVec.size() <= 1){
	// 	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++){


	// 		int i, j, rGaps, gGaps, num ,glen =-1;
	// 		SeedPair_t SeedPair;
	// 		if (iter->SeedVec.size() > 1){	
	// 			// printf("new seedvec");
	// 			// for(int x=0;x <iter->SeedVec.size();x++){
	// 			// 	printf("%d\t",iter->SeedVec[x].rPos);
	// 			// }
	// 			// printf("\n");


	// 			SeedPair.bSimple = false;
	// 			num = (int)iter->SeedVec.size();

	// 			for (i = 0, j = 1; j < num; i++, j++){
	// 				rGaps = iter->SeedVec[j].rPos - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen); 
	// 				if (rGaps < 0) rGaps = 0;
	// 				gGaps = iter->SeedVec[j].gPos - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen); 
	// 				if (gGaps < 0) gGaps = 0;
	// 				//printf("check %d and %d: rGaps=%d, gGaps=%d\n", i + 1, j + 1, rGaps, gGaps); fflush(stdout);
	// 				if (rGaps > 0 || gGaps > 0){
	// 					SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
	// 					SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
	// 					SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
	// 					SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
	// 					//printf("Add normal pair:\nR1[%d-%d] G1[%lld-%lld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%lld-%lld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
	// 					iter->SeedVec.push_back(SeedPair);
	// 				}
	// 			}
	// 			if ((int)iter->SeedVec.size() > num) inplace_merge(iter->SeedVec.begin(), iter->SeedVec.begin()+num, iter->SeedVec.end(), CompByGenomePos);
	// 			//if(bDebugMode) printf("After filling gaps\n"), ShowSeedInfo(SeedVec);
	// 		}
	// 		// Check missing blocks at both ends
	// 		if (iter->SeedVec.size() > 0){
	// 			i = 0;
	// 			//printf("Identify heading pair R[%d-%d]=%d, G[%lld-%lld]=%d\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].rLen, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].gLen); fflush(stdout);
	// 			rGaps = iter->SeedVec[i].rPos > 0 ? iter->SeedVec[i].rPos : 0;
	// 			gGaps = glen > 0 ? iter->SeedVec[i].gPos : rGaps;
	// 			//printf("rGaps=%d, gGaps=%d\n", rGaps, gGaps); fflush(stdout);
	// 			if (rGaps > 0 || gGaps > 0){
	// 				SeedPair.rPos = 0;
	// 				SeedPair.gPos = iter->SeedVec[i].gPos - gGaps;
	// 				if (SeedPair.gPos < 0) { SeedPair.gPos = 0; gGaps += SeedPair.gPos; }
	// 				SeedPair.PosDiff = SeedPair.gPos;
	// 				SeedPair.bSimple = false;
	// 				SeedPair.rLen = rGaps;
	// 				SeedPair.gLen = gGaps;
	// 				iter->SeedVec.insert(iter->SeedVec.begin(), SeedPair);
	// 				//if (bDebugMode) printf("Add missing head: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
	// 			}
	// 			i = (int)iter->SeedVec.size() - 1;
	// 			rGaps = rlen - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen);
	// 			gGaps = glen > 0 ? glen - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen) : rGaps;

	// 			if (rGaps > 0 || gGaps > 0){
	// 				SeedPair.bSimple = false;
	// 				SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
	// 				SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
	// 				SeedPair.rLen = rGaps;
	// 				SeedPair.gLen = gGaps;
	// 				iter->SeedVec.push_back(SeedPair);
	// 				//if (bDebugMode) printf("Add missing tail: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
	// 			}
	// 		}

			
	// 	}


	// 	return ;
	// } else{
	// 	int score1, score2;

	// 	score1 = score2 = 0;
	// 	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++)
	// 	{
	// 		if (iter->Score > score2)
	// 		{
	// 			if (iter->Score >= score1)
	// 			{
	// 				score2 = score1;
	// 				score1 = iter->Score;
	// 			}
	// 			else score2 = iter->Score;
	// 		}
	// 	}
	// 	if (bPacBioData || score1 == score2 || score1 - score2 > 20) ali_thr = score1;
	// 	else ali_thr = score2;

	// 	//if (bDebugMode) printf("Candidate score threshold = %d\n", thr);
	// 	// for(int x=0;x<AlignmentVec.size();x++){
	// 	// 	fprintf(stderr,"ali_thr:%d score:%d\n",ali_thr,AlignmentVec[x].Score);
	// 	// }
	// 	for (iter = AlignmentVec.begin(); iter != AlignmentVec.end(); iter++){
	// 		if (iter->Score < ali_thr){
	// 			iter->Score = 0;
	// 		} else{

	// 			int i, j, rGaps, gGaps, num ,glen =-1;
	// 			SeedPair_t SeedPair;
	// 			if (iter->SeedVec.size() > 1){	
	// 				// printf("new seedvec");
	// 				// for(int x=0;x <iter->SeedVec.size();x++){
	// 				// 	printf("%d\t",iter->SeedVec[x].rPos);
	// 				// }
	// 				// printf("\n");
					
	// 				// printf("Score:%d SeedVec.size():%d\n",iter->Score,iter->SeedVec.size());
	// 				// for(int y=0;y<iter->SeedVec.size();y++){
	// 				// 	printf("bSimple:%d %d-%d\n",iter->SeedVec[y].bSimple,iter->SeedVec[y].rPos,iter->SeedVec[y].rPos+iter->SeedVec[y].rLen);
	// 				// }
					

	// 				SeedPair.bSimple = false;
	// 				num = (int)iter->SeedVec.size();

	// 				for (i = 0, j = 1; j < num; i++, j++){
	// 					rGaps = iter->SeedVec[j].rPos - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen); 
	// 					if (rGaps < 0) rGaps = 0;
	// 					gGaps = iter->SeedVec[j].gPos - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen); 
	// 					if (gGaps < 0) gGaps = 0;
	// 					//printf("check %d and %d: rGaps=%d, gGaps=%d\n", i + 1, j + 1, rGaps, gGaps); fflush(stdout);
	// 					if (rGaps > 0 || gGaps > 0){
	// 						SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
	// 						SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
	// 						SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
	// 						SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
	// 						// printf("Add normal pair:\nR1[%d-%d] G1[%lld-%lld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%lld-%lld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
	// 						iter->SeedVec.push_back(SeedPair);
	// 					}
	// 				}
	// 				if ((int)iter->SeedVec.size() > num){
	// 					inplace_merge(iter->SeedVec.begin(), iter->SeedVec.begin()+num, iter->SeedVec.end(), CompByGenomePos);
	// 				} 
	// 				//if(bDebugMode) printf("After filling gaps\n"), ShowSeedInfo(SeedVec);
	// 			}
	// 			// Check missing blocks at both ends
	// 			if (iter->SeedVec.size() > 0){
	// 				i = 0;
	// 				//printf("Identify heading pair R[%d-%d]=%d, G[%lld-%lld]=%d\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].rLen, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].gLen); fflush(stdout);
	// 				rGaps = iter->SeedVec[i].rPos > 0 ? iter->SeedVec[i].rPos : 0;
	// 				gGaps = glen > 0 ? iter->SeedVec[i].gPos : rGaps;
	// 				//printf("rGaps=%d, gGaps=%d\n", rGaps, gGaps); fflush(stdout);
	// 				if (rGaps > 0 || gGaps > 0){
	// 					SeedPair.rPos = 0;
	// 					SeedPair.gPos = iter->SeedVec[i].gPos - gGaps;
	// 					if (SeedPair.gPos < 0) { SeedPair.gPos = 0; gGaps += SeedPair.gPos; }
	// 					SeedPair.PosDiff = SeedPair.gPos;
	// 					SeedPair.bSimple = false;
	// 					SeedPair.rLen = rGaps;
	// 					SeedPair.gLen = gGaps;
	// 					iter->SeedVec.insert(iter->SeedVec.begin(), SeedPair);
	// 					//if (bDebugMode) printf("Add missing head: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
	// 				}
	// 				i = (int)iter->SeedVec.size() - 1;
	// 				rGaps = rlen - (iter->SeedVec[i].rPos + iter->SeedVec[i].rLen);
	// 				gGaps = glen > 0 ? glen - (iter->SeedVec[i].gPos + iter->SeedVec[i].gLen) : rGaps;

	// 				if (rGaps > 0 || gGaps > 0){
	// 					SeedPair.bSimple = false;
	// 					SeedPair.rPos = iter->SeedVec[i].rPos + iter->SeedVec[i].rLen;
	// 					SeedPair.gPos = iter->SeedVec[i].gPos + iter->SeedVec[i].gLen;
	// 					SeedPair.rLen = rGaps;
	// 					SeedPair.gLen = gGaps;
	// 					iter->SeedVec.push_back(SeedPair);
	// 					//if (bDebugMode)printf("Add missing tail: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
	// 				}
	// 			}

	// 		}
	// 	}
	// }


	// 	printf("finally \n");
	// 	for(int x=0;x<AlignmentVec.size();x++){
	// 		if(AlignmentVec[x].Score >0){
	// 			for(int y=0;y<AlignmentVec[x].SeedVec.size();y++){
	// 				printf("bSimple:%d %d-%d\n",AlignmentVec[x].SeedVec[y].bSimple,AlignmentVec[x].SeedVec[y].rPos,AlignmentVec[x].SeedVec[y].rPos+AlignmentVec[x].SeedVec[y].rLen);
	// 			}
	// 			fprintf(stderr,"ali_thr:%d score:%d\n",ali_thr,AlignmentVec[x].Score);			
	// 		}
	// 	}
	// // return AlignmentVec;
}



int GapPenalty(vector<pair<int, char> >& cigar_vec)
{
	int GP = 0;
	for (vector<pair<int, char> >::iterator iter = cigar_vec.begin(); iter != cigar_vec.end(); iter++)
	{
		if (iter->second == 'I' || iter->second == 'D') GP += iter->first;
	}
	//if (bDebugMode) printf("GapPenalty=%d\n", GP), fflush(stdout);

	return GP;
}





bool CompByGenomePos(const SeedPair_t& p1, const SeedPair_t& p2)
{
	if (p1.gPos == p2.gPos) return p1.rPos < p2.rPos;
	else return p1.gPos < p2.gPos;
}


bool CompByFirstInt(const pair<int, int>& p1, const pair<int, int>& p2)
{
	return p1.first < p2.first;
}
 

void RemoveNullSeeds(vector<SeedPair_t>& SeedVec)
{
	for (vector<SeedPair_t>::iterator iter = SeedVec.begin(); iter != SeedVec.end();)
	{
		if (iter->rLen == 0) iter = SeedVec.erase(iter);
		else iter++;
	}
}
//rPos一样的统统删掉
void RemoveTandemRepeatSeeds(vector<SeedPair_t>& SeedVec)
{
	// printf("RemoveTandemRepeatSeeds");
	// for(int x=0;x <SeedVec.size();x++){
	// 	printf("%d\t",SeedVec[x].rPos);
	// }
	// printf("\n");
	int i, j, k, num;
	bool bTandemRepeat = false;
	vector<pair<int, int> > vec;
	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank
	// printf("SeedVec add:%ld num:%d \n",num,SeedVec);

	for (i = 0; i < num;) // identify all repetitive rPos
	{
		j = i + 1; while (j < num && vec[j].first == vec[i].first) j++;
		if (j - i > 1)
		{
			bTandemRepeat = true;
			// printf("Tandem repeat found: rPos = %d\n", vec[i].first);
			for (k = i; k < j; k++) SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
		}
		i = j;
	}
	//vec.clear(); vector<pair<int, int> >().swap(vec);
	if (bTandemRepeat) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After tandem repeat checking\n"), ShowSeedInfo(SeedVec);
}
int IdentifyTranslocationRange(int i, int num, vector<pair<int, int> >& vec, vector<SeedPair_t>& SeedVec)
{
	int j, max_idx = vec[i].second;

	for (j = i + 1; j <= max_idx; j++)
	{
		if (vec[j].second > max_idx) max_idx = vec[j].second;
	}
	return max_idx;
}

void RemoveTranslocatedSeeds(vector<SeedPair_t>& SeedVec)
{
	int i, j, k, s1, s2, num;
	bool bTranslocation = false;
	vector<pair<int, int> > vec;

	num = (int)SeedVec.size(); if (num < 2) return;
	vec.resize(num);
	for (i = 0; i < num; i++) vec[i] = make_pair(SeedVec[i].rPos, i);
	sort(vec.begin(), vec.end(), CompByFirstInt); //first:rPos, second: gPos rank

	// checking translocation
	for (i = 0; i < num; i++)
	{
		if (vec[i].first != SeedVec[i].rPos)
		{
			bTranslocation = true;
			j = IdentifyTranslocationRange(i, num, vec, SeedVec);
			s1 = 0; s2 = 0;
			for (k = i; k <= j; k++)
			{
				if (k < vec[k].second) s1 += SeedVec[vec[k].second].rLen;
				else s2 += SeedVec[vec[k].second].rLen;
			}
			if (s1 > s2)
			{
				for (k = i; k <= j; k++)
					if (k > vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			else
			{
				for (k = i; k <= j; k++)
					if (k < vec[k].second)
					{
						//printf("set #%d seed to zero length\n", vec[k].second + 1);
						SeedVec[vec[k].second].rLen = SeedVec[vec[k].second].gLen = 0;
					}
			}
			i = j;
		}
	}
	//vec.clear(); vector<pair<int, int> >().swap(vec);
	if (bTranslocation) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("After translocation checking\n"), ShowSeedInfo(SeedVec);
}

bool CheckSeedOverlapping(SeedPair_t& p1, SeedPair_t& p2)
{
	int iOverlap;
	bool bMaster = true;

	//printf("[1]: p1: r[%d-%d]=%d g[%lld-%lld]=%d vs p2: r[%d-%d]=%d g[%lld-%lld]=%d\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	if ((iOverlap = p1.rPos + p1.rLen - p2.rPos) > 0)
	{
		if (p1.rLen < p2.rLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	if ((p1.rLen > 0 && p2.rLen > 0) && (iOverlap = p1.gPos + p1.gLen - p2.gPos) > 0)
	{
		if (p1.gLen < p2.gLen)
		{
			bMaster = false;
			if (p1.rLen > iOverlap)
			{
				p1.gLen = (p1.rLen -= iOverlap);
			}
			else p1.rLen = p1.gLen = 0;
		}
		else
		{
			if (p2.rLen > iOverlap)
			{
				p2.rPos += iOverlap; p2.gPos += iOverlap;
				p2.gLen = (p2.rLen -= iOverlap);
			}
			else p2.rLen = p2.gLen = 0;
		}
	}
	//printf("[2]: p1: r[%d-%d]=%d g[%lld-%lld]=%d vs p2: r[%d-%d]=%d g[%lld-%lld]=%d\n\n", p1.rPos, p1.rPos + p1.rLen - 1, p1.rLen, p1.gPos, p1.gPos + p1.gLen - 1, p1.gLen, p2.rPos, p2.rPos + p2.rLen - 1, p2.rLen, p2.gPos, p2.gPos + p2.gLen - 1, p2.gLen);
	return bMaster;
}

int LocateThePreviousSeedIdx(int i, vector<SeedPair_t>& SeedVec)
{
	while (i > 0 && SeedVec[i].rLen == 0) i--;
	if (i < 0) return 0;
	else return i;
}

void CheckOverlappingSeeds(vector<SeedPair_t>& SeedVec)
{
	int64_t gEnd;
	int i, j, num, rEnd;
	bool bNullSeed = false;

	num = (int)SeedVec.size(); if (num < 2) return;
	for (i = 0; i < num;)
	{
		if (SeedVec[i].rLen > 0)
		{
			rEnd = SeedVec[i].rPos + SeedVec[i].rLen - 1;
			gEnd = SeedVec[i].gPos + SeedVec[i].gLen - 1;

			//printf("overlap check for seed#%d, rEnd=%d, gEned=%lld\n", i + 1, rEnd, gEnd);
			for (j = i + 1; j < num; j++)
			{
				if (SeedVec[j].rLen == 0) continue;
				if (rEnd < SeedVec[j].rPos && gEnd < SeedVec[j].gPos) break;
				//printf("\ttest seed#%d\n", j + 1);
				if (CheckSeedOverlapping(SeedVec[i], SeedVec[j]) == false) break;
			}
			if (SeedVec[i].rLen == 0)
			{
				bNullSeed = true;
				i = LocateThePreviousSeedIdx(i - 1, SeedVec);
			}
			else i++;
		}
		else
		{
			bNullSeed = true; i++;
		}
	}
	if (bNullSeed) RemoveNullSeeds(SeedVec);
	//if (bDebugMode) printf("after overlap checking\n"), ShowSeedInfo(SeedVec);
}

void IdentifyNormalPairs(int rlen, int glen, vector<SeedPair_t>& SeedVec)
{
	SeedPair_t SeedPair;
	int i, j, rGaps, gGaps, num;

	if (SeedVec.size() > 1)
	{	
		// printf("new seedvec");
		// for(int x=0;x <SeedVec.size();x++){
		// 	printf("%d\t",SeedVec[x].rPos);
		// }
		// printf("\n");

		RemoveTandemRepeatSeeds(SeedVec);
		RemoveTranslocatedSeeds(SeedVec);
		CheckOverlappingSeeds(SeedVec);

		SeedPair.bSimple = false;
		num = (int)SeedVec.size();

		for (i = 0, j = 1; j < num; i++, j++){
			rGaps = SeedVec[j].rPos - (SeedVec[i].rPos + SeedVec[i].rLen); 
			if (rGaps < 0) rGaps = 0;
			gGaps = SeedVec[j].gPos - (SeedVec[i].gPos + SeedVec[i].gLen); 
			if (gGaps < 0) gGaps = 0;
			//printf("check %d and %d: rGaps=%d, gGaps=%d\n", i + 1, j + 1, rGaps, gGaps); fflush(stdout);
			if (rGaps > 0 || gGaps > 0){
				SeedPair.rPos = SeedVec[i].rPos + SeedVec[i].rLen;
				SeedPair.gPos = SeedVec[i].gPos + SeedVec[i].gLen;
				SeedPair.PosDiff = SeedPair.gPos - SeedPair.rPos;
				SeedPair.rLen = rGaps; SeedPair.gLen = gGaps;
				//printf("Add normal pair:\nR1[%d-%d] G1[%lld-%lld]\nNR[%d-%d]=%d NG[%ld-%ld]=%d\nR2[%d-%d] G2[%lld-%lld]\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen, SeedVec[j].rPos, SeedVec[j].rPos + SeedVec[j].rLen - 1, SeedVec[j].gPos, SeedVec[j].gPos + SeedVec[j].gLen - 1); fflush(stdout);
				SeedVec.push_back(SeedPair);
			}
		}
		if ((int)SeedVec.size() > num) inplace_merge(SeedVec.begin(), SeedVec.begin()+num, SeedVec.end(), CompByGenomePos);
		//if(bDebugMode) printf("After filling gaps\n"), ShowSeedInfo(SeedVec);
	}
	// Check missing blocks at both ends
	if (SeedVec.size() > 0)
	{
		i = 0;
		//printf("Identify heading pair R[%d-%d]=%d, G[%lld-%lld]=%d\n", SeedVec[i].rPos, SeedVec[i].rPos + SeedVec[i].rLen - 1, SeedVec[i].rLen, SeedVec[i].gPos, SeedVec[i].gPos + SeedVec[i].gLen - 1, SeedVec[i].gLen); fflush(stdout);
		rGaps = SeedVec[i].rPos > 0 ? SeedVec[i].rPos : 0;
		gGaps = glen > 0 ? SeedVec[i].gPos : rGaps;
		//printf("rGaps=%d, gGaps=%d\n", rGaps, gGaps); fflush(stdout);
		if (rGaps > 0 || gGaps > 0)
		{
			SeedPair.rPos = 0;
			SeedPair.gPos = SeedVec[i].gPos - gGaps;
			if (SeedPair.gPos < 0) { SeedPair.gPos = 0; gGaps += SeedPair.gPos; }
			SeedPair.PosDiff = SeedPair.gPos;
			SeedPair.bSimple = false;
			SeedPair.rLen = rGaps;
			SeedPair.gLen = gGaps;
			SeedVec.insert(SeedVec.begin(), SeedPair);

			//if (bDebugMode) printf("Add missing head: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
		}
		i = (int)SeedVec.size() - 1;
		rGaps = rlen - (SeedVec[i].rPos + SeedVec[i].rLen);
		gGaps = glen > 0 ? glen - (SeedVec[i].gPos + SeedVec[i].gLen) : rGaps;

		if (rGaps > 0 || gGaps > 0)
		{
			SeedPair.bSimple = false;
			SeedPair.rPos = SeedVec[i].rPos + SeedVec[i].rLen;
			SeedPair.gPos = SeedVec[i].gPos + SeedVec[i].gLen;
			SeedPair.rLen = rGaps;
			SeedPair.gLen = gGaps;
			SeedVec.push_back(SeedPair);

			//if (bDebugMode) printf("Add missing tail: r[%d-%d]=%d, g[%lld-%lld]=%d\n", SeedPair.rPos, SeedPair.rPos + SeedPair.rLen - 1, SeedPair.rLen, SeedPair.gPos, SeedPair.gPos + SeedPair.gLen - 1, SeedPair.gLen);
		}
	}
	//if (bDebugMode) printf("after generating normal pairs\n"), ShowSeedInfo(SeedVec);
}



string GenerateCIGAR(vector<pair<int, char> >& cigar_vec)
{
	int i, num, c;
	char state, buf[8];
	string cigar_str;

	for (state = '\0', num = (int)cigar_vec.size(), c = 0, i = 0; i != num; i++)
	{
		if (cigar_vec[i].second != state)
		{
			if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

			c = cigar_vec[i].first; state = cigar_vec[i].second;
		}
		else c += cigar_vec[i].first;
	}
	if (c > 0) sprintf(buf, "%d%c", c, state), cigar_str += buf;

	if (bDebugMode) printf("CIGAR=%s\n\n\n", cigar_str.c_str());

	return cigar_str;
}
  
Coordinate_t GenCoordinateInfo(bool bFirstRead, int64_t gPos, int64_t end_gPos, vector<pair<int, char> >& cigar_vec)
{
	Coordinate_t coor;
	map<int64_t, int>::iterator iter;

	if (gPos < GenomeSize) // forward strand
	{
		if (bFirstRead) coor.bDir = true;
		else coor.bDir = false;

		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = gPos + 1;
		}
		else
		{
			iter = ChrLocMap.lower_bound(gPos);
			coor.ChromosomeIdx = iter->second;
			coor.gPos = gPos + 1 - ChromosomeVec[coor.ChromosomeIdx].FowardLocation;
			//if (bDebugMode) printf("matched chr=%s, loc=%lld, gPos: %lld -> %lld\n", ChromosomeVec[coor.ChromosomeIdx].name, ChromosomeVec[coor.ChromosomeIdx].FowardLocation, gPos, coor.gPos);
		}
	}
	else
	{
		if (bFirstRead) coor.bDir = false;
		else coor.bDir = true;

		reverse(cigar_vec.begin(), cigar_vec.end());

		if (iChromsomeNum == 1)
		{
			coor.ChromosomeIdx = 0;
			coor.gPos = TwoGenomeSize - end_gPos;
		}
		else
		{
			iter = ChrLocMap.lower_bound(gPos);
			coor.gPos = iter->first - end_gPos + 1; coor.ChromosomeIdx = iter->second;
			//if(bDebugMode) printf("matched chr=%s, loc=%lld, gPos: %lld -> %lld\n", ChromosomeVec[coor.ChromosomeIdx].name, ChromosomeVec[coor.ChromosomeIdx].ReverseLocation, gPos, coor.gPos);
		}
	}
	//if (bDebugMode) printf("gPos: %lld --> %lld %s\n", gPos, coor.gPos, (coor.bDir? "Forward":"Reverse"));

	coor.CIGAR = GenerateCIGAR(cigar_vec);
	// printf("coor.CIGAR:%s\n", coor.CIGAR.c_str());

	return coor;
}

void GenMappingReport(bool bFirstRead, ReadItem_t& read, vector<AlignmentCandidate_t>& AlignmentVec)
{
	int i, j, num, s;
	map<int, int>::iterator iter;
	vector<pair<int, char> > cigar_vec;
	map<int64_t, int>::iterator ChrIter;

	//if (bDebugMode) printf("\n\n%s\nGenerate alignment for read %s (%d cans)\n", string().assign(100, '=').c_str(), read.header, (int)AlignmentVec.size()), fflush(stdout);
	// fprintf(stderr,"%s-%d\n",read.header,bFirstRead);
	read.score = read.sub_score = read.iBestAlnCanIdx = 0;
	if ((read.CanNum = (int)AlignmentVec.size()) > 0){
		read.AlnReportArr = new AlignmentReport_t[read.CanNum];
		for (i = 0; i != read.CanNum; i++){
			int NM =0;
			read.AlnReportArr[i].isPrime=0;
			read.AlnReportArr[i].isSecondary=0;
			read.AlnReportArr[i].SamFlag =0;
			read.AlnReportArr[i].AlnScore = 0;
			read.AlnReportArr[i].PairedAlnCanIdx = AlignmentVec[i].PairedAlnCanIdx;

			if (AlignmentVec[i].Score == 0)  ;
			if (bPacBioData && read.score > 0)
			{
				read.sub_score = read.score;
				continue;
			}
			// fprintf(stderr,"before\n");
			// for(int x=0;x<AlignmentVec[i].SeedVec.size();x++){
			// 	fprintf(stderr,"simple:%d %d-%d\n",AlignmentVec[i].SeedVec[x].bSimple,AlignmentVec[i].SeedVec[x].rPos,AlignmentVec[i].SeedVec[x].rLen);
			// 	fprintf(stderr,"%ld-%d\n",AlignmentVec[i].SeedVec[x].gPos,AlignmentVec[i].SeedVec[x].gLen);

			// }
			IdentifyNormalPairs(read.rlen, -1, AlignmentVec[i].SeedVec); // fill missing framgment pairs (normal pairs) between simple pairs
			// fprintf(stderr,"after\n");
			// for(int x=0;x<AlignmentVec[i].SeedVec.size();x++){
			// 	fprintf(stderr,"simple:%d %d-%d\n",AlignmentVec[i].SeedVec[x].bSimple,AlignmentVec[i].SeedVec[x].rPos,AlignmentVec[i].SeedVec[x].rPos+AlignmentVec[i].SeedVec[x].rLen);
			// 	fprintf(stderr,"%ld-%ld\n",AlignmentVec[i].SeedVec[x].gPos,AlignmentVec[i].SeedVec[x].gLen+AlignmentVec[i].SeedVec[x].gPos);

			// }
			
			// if (bDebugMode)
			// {
			// 	printf("Process candidate#%d (Score = %d, SegmentPair#=%d): \n", i + 1, AlignmentVec[i].Score, (int)AlignmentVec[i].SeedVec.size());
			// 	ShowSeedInfo(AlignmentVec[i].SeedVec);
			// }
			if (CheckCoordinateValidity(AlignmentVec[i].SeedVec) == false) continue;

			cigar_vec.clear();
			for (num = (int)AlignmentVec[i].SeedVec.size(), j = 0; j != num; j++){
				if (AlignmentVec[i].SeedVec[j].rLen == 0 && AlignmentVec[i].SeedVec[j].gLen == 0) continue;
				else if (AlignmentVec[i].SeedVec[j].bSimple)
				{
					//if (bDebugMode) ShowFragmentPair(AlignmentVec[i].SeedVec[j]);
					cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'M'));
					read.AlnReportArr[i].AlnScore += AlignmentVec[i].SeedVec[j].rLen;
				}
				else
				{
					//if (bDebugMode) printf("Check normal pair#%d: R[%d-%d]=%d G[%lld-%lld]=%d\n", j + 1, AlignmentVec[i].SeedVec[j].rPos, AlignmentVec[i].SeedVec[j].rPos + AlignmentVec[i].SeedVec[j].rLen - 1, AlignmentVec[i].SeedVec[j].rLen, AlignmentVec[i].SeedVec[j].gPos, AlignmentVec[i].SeedVec[j].gPos + AlignmentVec[i].SeedVec[j].gLen - 1, AlignmentVec[i].SeedVec[j].gLen);
					if (j == 0)
					{
						if (AlignmentVec[i].SeedVec[0].rLen > 3000)
						{
							cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[0].rLen, 'S'));
							AlignmentVec[i].SeedVec[0].gPos = AlignmentVec[i].SeedVec[1].gPos;
							AlignmentVec[i].SeedVec[0].gLen = 0;
						}else{
							s = ProcessHeadSequencePair(read.seq, AlignmentVec[i].SeedVec[0], cigar_vec);
							read.AlnReportArr[i].AlnScore += s;
							if (s == 0){
								AlignmentVec[i].SeedVec[0].gPos = AlignmentVec[i].SeedVec[1].gPos;
								AlignmentVec[i].SeedVec[0].gLen = 0;
							}else{
								NM += (AlignmentVec[i].SeedVec[0].rLen-s);
							}
						}
					}
					else if (j == num - 1)
					{
						if (AlignmentVec[i].SeedVec[j].rLen > 3000)
						{
							
							cigar_vec.push_back(make_pair(AlignmentVec[i].SeedVec[j].rLen, 'S'));
							AlignmentVec[i].SeedVec[j].gPos = AlignmentVec[i].SeedVec[j-1].gPos + AlignmentVec[i].SeedVec[j - 1].gLen;
							AlignmentVec[i].SeedVec[j].gLen = 0;
						}
						else
						{
							s = ProcessTailSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
							read.AlnReportArr[i].AlnScore += s;
							if (s == 0)
							{
								AlignmentVec[i].SeedVec[j].gPos = AlignmentVec[i].SeedVec[j - 1].gPos + AlignmentVec[i].SeedVec[j - 1].gLen;
								AlignmentVec[i].SeedVec[j].gLen = 0;
							}else{
								NM += (AlignmentVec[i].SeedVec[j].rLen-s);

							}
						}
					}
					else
					{	s =ProcessNormalSequencePair(read.seq, AlignmentVec[i].SeedVec[j], cigar_vec);
						read.AlnReportArr[i].AlnScore += s;
						if(s>0){
							NM += (AlignmentVec[i].SeedVec[j].rLen-s);
						}
					}
				}
			}
			//if (bDebugMode) 			printf("Alignment score = %d (rlen=%d) NM=%d\n", read.AlnReportArr[i].AlnScore, read.rlen ,NM), fflush(stdout);



			read.AlnReportArr[i].NM =NM;
			if (!bPacBioData && cigar_vec.size() > 1)
			{
				read.AlnReportArr[i].AlnScore -= GapPenalty(cigar_vec);
				if (read.AlnReportArr[i].AlnScore <= 0)
				{
					read.AlnReportArr[i].AlnScore = 0;
					continue;
				}
			}
			if(cigar_vec.size() == 0 || (read.AlnReportArr[i].coor = GenCoordinateInfo(bFirstRead, AlignmentVec[i].SeedVec[0].gPos, (AlignmentVec[i].SeedVec[num - 1].gPos + AlignmentVec[i].SeedVec[num - 1].gLen - 1), cigar_vec)).gPos <= 0){
				read.AlnReportArr[i].AlnScore = 0;
			} 

			if (read.AlnReportArr[i].AlnScore > read.score){
				// read.AlnReportArr[read.iBestAlnCanIdx].isSecondary=1;
				// read.AlnReportArr[read.iBestAlnCanIdx].isPrime=0;

				read.iBestAlnCanIdx = i;
				read.sub_score = read.score;
				read.score = read.AlnReportArr[i].AlnScore;

				// read.AlnReportArr[i].isPrime =1;
				// read.AlnReportArr[i].isSecondary =0;

			}else if (read.AlnReportArr[i].AlnScore == read.score){
				// read.sub_score = read.score;
				// read.AlnReportArr[i].isPrime =1;
				// read.AlnReportArr[i].isSecondary =0;				
				if (!bMultiHit && ChromosomeVec[read.AlnReportArr[i].coor.ChromosomeIdx].len > ChromosomeVec[read.AlnReportArr[read.iBestAlnCanIdx].coor.ChromosomeIdx].len){
					read.iBestAlnCanIdx = i;
				} 
			}else if(read.AlnReportArr[i].AlnScore >= read.sub_score){
				read.sub_score = read.AlnReportArr[i].AlnScore;
			}
			
		}

		for(int x=0;x<read.CanNum;x++){
			if(read.AlnReportArr[x].AlnScore ==read.score){
				read.AlnReportArr[x].isPrime =1;
				read.AlnReportArr[x].isSecondary =0;	
			}else if(read.AlnReportArr[x].AlnScore ==read.sub_score){
				read.AlnReportArr[x].isPrime =0;
				read.AlnReportArr[x].isSecondary =1;
			}
			// printf("%d-%s-%lld-%s-%d-NM:%d-prime:%d-secondray:%d\n",read.AlnReportArr[x].AlnScore,ChromosomeVec[read.AlnReportArr[x].coor.ChromosomeIdx].name,read.AlnReportArr[x].coor.gPos,read.AlnReportArr[x].coor.CIGAR.c_str(),read.rlen -read.AlnReportArr[x].AlnScore,read.AlnReportArr[x].NM,read.AlnReportArr[x].isPrime,read.AlnReportArr[x].isSecondary);
		}


	}
	else
	{
		read.CanNum = 1; read.iBestAlnCanIdx = 0;
		read.AlnReportArr = new AlignmentReport_t[1];
		read.AlnReportArr[0].AlnScore = 0;
		read.AlnReportArr[0].PairedAlnCanIdx = -1;
	}


	// for(int x=0;x<read.CanNum;x++){
	// 	printf("%d-%s-%lld-%s-%d-NM:%d-prime:%d-secondray:%d\n",read.AlnReportArr[x].AlnScore,ChromosomeVec[read.AlnReportArr[x].coor.ChromosomeIdx].name,read.AlnReportArr[x].coor.gPos,read.AlnReportArr[x].coor.CIGAR.c_str(),read.rlen -read.AlnReportArr[x].AlnScore,read.AlnReportArr[x].NM,read.AlnReportArr[x].isPrime,read.AlnReportArr[x].isSecondary);
	// }

}



