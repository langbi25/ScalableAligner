



#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <string.h>
#define iShift 30

using namespace std;

bool bShowWrongCase=false;
bool bShowUnmapped=true;
bool bShowBadAlign=false;

int TotalQuery=0;

string indexSAMFile;

struct posInfo
{
    string TrueChName;	   
	unsigned int TruePos1;   
	unsigned int TruePos2;   
};

struct resultInfo
{
    int readMAPQ;	   
	int readLen;   
	vector< pair<string,unsigned int> > predictPos;   
};


map<string, posInfo> QueryMap;


map<string, resultInfo > ResultMap;

vector< pair<string,unsigned int> >  truePos;


bool CheckPosConsistency(int Length, int TrueLocation, int PredictedLocation){
	if(TrueLocation >= PredictedLocation && TrueLocation - PredictedLocation < (iShift+ Length)) return true;
	else if (TrueLocation < PredictedLocation && PredictedLocation - TrueLocation < (iShift + Length)) return true;
	else return false;
}






void getAnswer(string answerSamFileName){

printf("loading answer QueryMap ...\n");

	fstream file;
	stringstream ss;
	bool bCorLocation, bMapped;
	unsigned int pos, TrueLocation;
	long long iTotalBase, iTotalCorBase;
	int iFlag, MAPQ, iReadNum=0, iUnmapped=0, iBadMAPQ=0;
	string str, Prevheader, Readheader, PosInfoStr, ChName, tmp, seq, region, CIGAR;
	int PreFlag=0;	

	file.open(answerSamFileName.c_str(), ios_base::in);

	if(file.is_open()){
		while(!file.eof())
		{
			getline(file, str); if(str=="") break;
			if(str[0]=='@') continue;

			ss.clear(); ss.str(str);
			ss >> Readheader >> iFlag >> ChName >> PosInfoStr >> MAPQ >> tmp >> tmp >> tmp >> tmp >> seq;
			if((iFlag >>6 )&1){
				// printf("%d read1\n",iFlag );
				if(Readheader != Prevheader){
					Prevheader = Readheader;
					QueryMap[Readheader].TrueChName =ChName;
					if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str());
					QueryMap[Readheader].TruePos1 =pos;

					iReadNum++;
				}

			}else{
				if((PreFlag >>6 )&1  && Readheader == Prevheader){
					if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str());
					QueryMap[Readheader].TruePos2 =pos;
					iReadNum++;
				}
				// printf("%d read2\n",iFlag);
			}
			PreFlag =iFlag;

		}
		file.close();
	}
	printf("finishing loading answer QueryMap with %d pairs of reads\n",QueryMap.size());

}

bool CheckPosConsistency_2(int Length, string PredictedChName, int PredictedLocation,string Readheader){

	
	if(QueryMap.find(Readheader) != QueryMap.end()){
		if(QueryMap[Readheader].TrueChName == PredictedChName){
			if(QueryMap[Readheader].TruePos1 >= PredictedLocation && QueryMap[Readheader].TruePos1 - PredictedLocation < (iShift+ Length)) return true;
			else if (QueryMap[Readheader].TruePos1 < PredictedLocation && PredictedLocation - QueryMap[Readheader].TruePos1 < (iShift + Length)) return true;
			
			else if(QueryMap[Readheader].TruePos2 >= PredictedLocation && QueryMap[Readheader].TruePos2 - PredictedLocation < (iShift+ Length)) return true;
			else if (QueryMap[Readheader].TruePos2 < PredictedLocation && PredictedLocation - QueryMap[Readheader].TruePos2 < (iShift + Length)) return true;
			else return false;
		}else{
			return false;
		}
	}else{
		printf("something go wrong..\n");
		return false;
	}

	
	return false;

}


// bool CheckPosConsistency_3(int Length,string Readheader){

// 	for(int i=0;i< PredictPos.size();i++){
// 		//先判断这个header有没有被记录
// 		if(QueryMap.find(Readheader) != QueryMap.end()){

// 			if(QueryMap[Readheader].TrueChName == PredictPos[i].first){
// 				if(QueryMap[Readheader].TruePos1 >= PredictPos[i].second && QueryMap[Readheader].TruePos1 - PredictPos[i].second < (iShift+ Length)) return true;
// 				else if (QueryMap[Readheader].TruePos1 < PredictPos[i].second && PredictPos[i].second - QueryMap[Readheader].TruePos1 < (iShift + Length)) return true;
				
// 				else if(QueryMap[Readheader].TruePos2 >= PredictPos[i].second && QueryMap[Readheader].TruePos2 - PredictPos[i].second < (iShift+ Length)) return true;
// 				else if (QueryMap[Readheader].TruePos2 < PredictPos[i].second && PredictPos[i].second - QueryMap[Readheader].TruePos2 < (iShift + Length)) return true;
// 				else continue;
// 			}else{
// 				continue;
// 			}
// 		}else{
// 			printf("something go wrong..\n");
// 			return false;
// 		}


// 	}


	
// 	return false;

// }


bool CheckPosConsistency_4(int Length,string Readheader){

	string header =Readheader.substr(0,Readheader.length()-2);

	// printf("Readheader:%s\n",Readheader.c_str());


	vector< pair<string,unsigned int> >  PredictPos =ResultMap[Readheader].predictPos;


	



	for(int i=0;i< PredictPos.size();i++){
		//先判断这个header有没有被记录
		if(QueryMap.find(header) != QueryMap.end()){

			if(QueryMap[header].TrueChName == PredictPos[i].first){
				if(QueryMap[header].TruePos1 >= PredictPos[i].second && QueryMap[header].TruePos1 - PredictPos[i].second < (iShift+ Length)) return true;
				else if (QueryMap[header].TruePos1 < PredictPos[i].second && PredictPos[i].second - QueryMap[header].TruePos1 < (iShift + Length)) return true;
				
				else if(QueryMap[header].TruePos2 >= PredictPos[i].second && QueryMap[header].TruePos2 - PredictPos[i].second < (iShift+ Length)) return true;
				else if (QueryMap[header].TruePos2 < PredictPos[i].second && PredictPos[i].second - QueryMap[header].TruePos2 < (iShift + Length)) return true;
				else continue;
			}else{
				continue;
			}
		}else{
			printf("something go wrong..\n");
			return false;
		}


	}


	
	return false;

}



void Evaluation(string SamFileName ){
	fstream file;
	stringstream ss;
	bool bCorLocation, bMapped;
	unsigned int pos, TrueLocation;
	long long iTotalBase, iTotalCorBase;
	int iFlag, MAPQ, iReadNum=0, iUnmapped=0, iCorLocation=0, iBadMAPQ=0;
	string str, Prevheader, Readheader, PosInfoStr, ChName, tmp, seq, region, CIGAR,XZ;

	unsigned int realPos1,realPos2;

	int PreFlag=0,Primary,xz_pos;
	string realChName,tmpReadheader;

	iTotalBase = iTotalCorBase = 0;

	pos = 0; Prevheader = "";
	bMapped=bCorLocation=false;


	int flag_of_first;
	int second_read;

	// getAnswer(answerSamFileName);


	file.open(SamFileName.c_str(), ios_base::in);

	if(file.is_open()){
		while(!file.eof()){
			getline(file, str); if(str=="") break;
			if(str[0]=='@') continue;

			ss.clear(); ss.str(str);
			ss >> Readheader >> iFlag >> ChName >> PosInfoStr >> MAPQ >> tmp >> tmp >> tmp >> tmp >> seq >> tmp;


						// printf("wdnmd:%s -%s\n",Prevheader.c_str() ,Readheader.c_str());


			if(Readheader != Prevheader){
				iReadNum++;

				flag_of_first =(iFlag >>6 )&1;
				second_read =0;

				if(ChName=="*"){
					iUnmapped++;
					if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
					PreFlag =iFlag;
					Prevheader =Readheader;
					// iReadNum++;
					continue;
				}
				// resultInfo resInfo;
				// // vector< pair<string,unsigned int> >  PredictPos;
				// resInfo.readLen =(int)seq.length();
				// resInfo.readMAPQ = MAPQ;
				// string tmpHeader =Readheader;
				// string readNo =tmpHeader.append("_1");
				// if(ChName!="*" && MAPQ==0) {
				// 	// printf("wdnmd:%d\n",str.c_str());
				// 	iBadMAPQ++;
				// }
				// if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;
				// if(ChName=="*"){
				// 	iUnmapped++;
				// 	if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
				// }else{
				// 	resInfo.predictPos.push_back(pair<string,unsigned int>(ChName,pos));
				// }


				// xz_pos=str.find("XA:Z:");
				// if( xz_pos!= string::npos){
				// 	XZ =str.substr(xz_pos+5);
				// 	// printf("XZ:%s\n",XZ.c_str());
				// 	const char s[2] = ";";

				// 	char *token=strtok((char *)XZ.c_str(), s);
				// 	string ans=token;
				// 	int fisrt_douhao =ans.find(",");
				// 	string chr =ans.substr(0,fisrt_douhao);
				// 	ans =ans.substr(fisrt_douhao+1);
				// 	unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
				// 	pos1-=1;
				// 	resInfo.predictPos.push_back(pair<string,unsigned int>(chr,pos1));
				// 	/* 继续获取其他的子字符串 */
				// 	while( token != NULL ) {
				// 		// printf( "%s\n", token );
						
				// 		token = strtok(NULL, s);
				// 		if(token != NULL){
				// 			string ans=token;
				// 			int fisrt_douhao =ans.find(",");
				// 			string chr =ans.substr(0,fisrt_douhao);
				// 			ans =ans.substr(fisrt_douhao+1);
				// 			unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
				// 			resInfo.predictPos.push_back(pair<string,unsigned int>(chr,pos1));
				// 		}
				// 	}

				// }
				// ResultMap[readNo] =resInfo;				

			}else{
						// printf("flag_of_first:%d -%d\n",flag_of_first ,(iFlag >>6 )&1);
				
				//说明是mate pair的另外一半
				if( flag_of_first != (iFlag >>6 )&1){

					if(!second_read){
						iReadNum++;
						second_read =1;

						if(ChName=="*"){
							iUnmapped++;
							if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
							PreFlag =iFlag;
							Prevheader =Readheader;
							// iReadNum++;
							continue;
						}
						// resultInfo resInfo;
						// // vector< pair<string,unsigned int> >  PredictPos;
						// resInfo.readLen =(int)seq.length();
						// resInfo.readMAPQ = MAPQ;

						// 					string tmpHeader =Readheader;
						// string readNo =tmpHeader.append("_2");	
						// // string readNo =Readheader.append("_2");
						// if(ChName!="*" && MAPQ==0) {
						// 	// printf("wdnmd:%d\n",str.c_str());
						// 	iBadMAPQ++;
						// }
						// if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;
						// if(ChName=="*"){
						// 	iUnmapped++;
						// 	if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
						// }else{
						// 	resInfo.predictPos.push_back(pair<string,unsigned int>(ChName,pos));
						// }


						// xz_pos=str.find("XA:Z:");
						// if( xz_pos!= string::npos){
						// 	XZ =str.substr(xz_pos+5);
						// 	// printf("XZ:%s\n",XZ.c_str());
						// 	const char s[2] = ";";

						// 	char *token=strtok((char *)XZ.c_str(), s);
						// 	string ans=token;
						// 	int fisrt_douhao =ans.find(",");
						// 	string chr =ans.substr(0,fisrt_douhao);
						// 	ans =ans.substr(fisrt_douhao+1);
						// 	unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
						// 	pos1-=1;
						// 	resInfo.predictPos.push_back(pair<string,unsigned int>(chr,pos1));
						// 	/* 继续获取其他的子字符串 */
						// 	while( token != NULL ) {
						// 		// printf( "%s\n", token );
								
						// 		token = strtok(NULL, s);
						// 		if(token != NULL){
						// 			string ans=token;
						// 			int fisrt_douhao =ans.find(",");
						// 			string chr =ans.substr(0,fisrt_douhao);
						// 			ans =ans.substr(fisrt_douhao+1);
						// 			unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
						// 			resInfo.predictPos.push_back(pair<string,unsigned int>(chr,pos1));
						// 		}
						// 	}

						// }
						// // printf("%s\n",readNo);
						// ResultMap[readNo] =resInfo;	

					}else{
						// string tmpHeader =Readheader;
						// string readNo =tmpHeader.append("_2");
						// // string readNo =Readheader.append("_2");

						// if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;
						// if(ChName=="*"){
						// 	iUnmapped++;
						// 	if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
						// }else{
						// 	ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(ChName,pos));
						// }


						// xz_pos=str.find("XA:Z:");
						// if( xz_pos!= string::npos){
						// 	XZ =str.substr(xz_pos+5);
						// 	// printf("XZ:%s\n",XZ.c_str());
						// 	const char s[2] = ";";

						// 	char *token=strtok((char *)XZ.c_str(), s);
						// 	string ans=token;
						// 	int fisrt_douhao =ans.find(",");
						// 	string chr =ans.substr(0,fisrt_douhao);
						// 	ans =ans.substr(fisrt_douhao+1);
						// 	unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
						// 	pos1-=1;
						// 	ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(chr,pos1));
						// 	/* 继续获取其他的子字符串 */
						// 	while( token != NULL ) {
						// 		// printf( "%s\n", token );
								
						// 		token = strtok(NULL, s);
						// 		if(token != NULL){
						// 			string ans=token;
						// 			int fisrt_douhao =ans.find(",");
						// 			string chr =ans.substr(0,fisrt_douhao);
						// 			ans =ans.substr(fisrt_douhao+1);
						// 			unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
						// 			ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(chr,pos1));
						// 		}
						// 	}

						// }
					}


				}else {
					//说明是第一条的附属

					// flag_of_first =(iFlag >>6 )&1;
					// second_read =0;
					// string tmpHeader =Readheader;
					// string readNo =tmpHeader.append("_1");


					// if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;
					// if(ChName=="*"){
					// 	iUnmapped++;
					// 	if(bShowUnmapped && bShowWrongCase) cout << Readheader << endl;
					// }else{
					// 	ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(ChName,pos));
					// }


					// xz_pos=str.find("XA:Z:");
					// if( xz_pos!= string::npos){
					// 	XZ =str.substr(xz_pos+5);
					// 	// printf("XZ:%s\n",XZ.c_str());
					// 	const char s[2] = ";";

					// 	char *token=strtok((char *)XZ.c_str(), s);
					// 	string ans=token;
					// 	int fisrt_douhao =ans.find(",");
					// 	string chr =ans.substr(0,fisrt_douhao);
					// 	ans =ans.substr(fisrt_douhao+1);
					// 	unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
					// 	pos1-=1;
					// 	ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(chr,pos1));
					// 	/* 继续获取其他的子字符串 */
					// 	while( token != NULL ) {
					// 		// printf( "%s\n", token );
							
					// 		token = strtok(NULL, s);
					// 		if(token != NULL){
					// 			string ans=token;
					// 			int fisrt_douhao =ans.find(",");
					// 			string chr =ans.substr(0,fisrt_douhao);
					// 			ans =ans.substr(fisrt_douhao+1);
					// 			unsigned int pos1 =(unsigned int)atoi(ans.substr(1,ans.find(",")-1).c_str());
					// 			ResultMap[readNo].predictPos.push_back(pair<string,unsigned int>(chr,pos1));
					// 		}
					// 	}

					// }
				}
			}

			Prevheader = Readheader;
			PreFlag =iFlag;

			// bMapped=bCorLocation=false;
			// if(PosInfoStr!="*") pos = (unsigned int)atoi(PosInfoStr.c_str()); pos-=1;

		}
		file.close();
	}

	for(auto it : ResultMap){
		// string header =it.first.substr(0,it.first.length()-2);
				// printf("%s\n",it.first.c_str());

		// printf("%s\t%d\n",it.first.substr(0,it.first.length()-2).c_str(),it.second.readLen);

		bCorLocation=false;

		if(bCorLocation==false && CheckPosConsistency_4(it.second.readLen,it.first)){
			bCorLocation=true;
		}
			

		if(bCorLocation==false){
			// cout << Readheader << endl;
			if(MAPQ > 0 && bShowWrongCase){
				cout << Readheader << endl;
			} 
		}else{
			iCorLocation++;
		} 
	}


	printf("ResultMap:%d\n",ResultMap.size());



	if(TotalQuery==0) TotalQuery=iReadNum;
	else if(iReadNum > (int)(TotalQuery*0.995)) TotalQuery=iReadNum;


	cerr << endl << endl << "filename=" << SamFileName << endl;
	cerr << "# of reads= " << TotalQuery << endl;
	if(iReadNum>0) cerr << "# of mapped reads= " << (iReadNum-iUnmapped) << " (" << (int)(10000*(1.0*(iReadNum-iUnmapped)/TotalQuery)+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "# of mapq_0=" << iBadMAPQ << " (" << (int)(10000*(1.0*iBadMAPQ/iReadNum)+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "precision= " << iCorLocation << " (" << (int)(10000*(1.0*iCorLocation/(iReadNum-iUnmapped))+0.5)/100.0 << "%)" << endl;
	if(iReadNum>0) cerr << "recall= " << iCorLocation << " (" << (int)(10000*(1.0*iCorLocation/TotalQuery)+0.5)/100.0 << "%)" << endl;
	cerr << endl << endl;

	// precision, recall
	FILE *rfile = fopen("sen_exp.csv", "a");
	if (iReadNum > 0) {
		// fprintf(rfile, "%s: %.1f\t%.1f\n", SamFileName.c_str(), (int)(1000 * (1.0*iCorLocation / (iReadNum - iUnmapped)) + 0.5) / 10.0, (int)(1000 * (1.0*iCorLocation / TotalQuery) + 0.5) / 10.0);
		fprintf(rfile, "%s: %.2f\n", SamFileName.c_str(), (int)(10000*(1.0*(iReadNum-iUnmapped)/TotalQuery)+0.5)/100.0);
	
	}
	else fprintf(rfile, "%s: 0\t0\n", SamFileName.c_str());

	fclose(rfile);
}

int main(int argc, char* argv[])
{
	string str;

	if (argc == 1)
	{
		printf("usage: %s SamFile\n\n", argv[0]);
	}
	else
	{
		for (int i = 2; i < argc; i++)
		{
			if ((str = argv[i]) == "-d") bShowWrongCase = true;

			// if ((str = argv[i]) == "-i") indexSAMFile = argv[++i];
		}
	}			
	Evaluation(argv[1] );

	return 0;
}
