#include <sys/stat.h>
#include "util.h"


extern "C"
{
	int bwa_idx_build(const char *fa, const char *prefix);
}
const char* VersionStr = "1.0.0";
bwt_t *bwt;
bwaidx_t *bwtIdx;
string indexPath,readFile1 ,readFile2;
char *outputFile ,*indexFile ,*refSeq;
int threadNum, maxInsertSize, MaxGaps, MinSeedLength, OutputFileFormat,OutputByOrder;
// bool bPairEnd ,bPacBioData,FastQFormat;
// int iThreadNum, MaxInsertSize;
bool bDebugMode, bPairEnd, bPacBioData, bMultiHit, gzCompressed, FastQFormat, bSilent ,bHash,bBlockRead;

void Usage()
{
	fprintf(stdout, "kart v (Hsin-Nan Lin & Wen-Lian Hsu)\n\n");
	fprintf(stdout, "Usage:  -i Index_Prefix -f <ReadFile_A1 ReadFile_B1 ...> [-f2 <ReadFile_A2 ReadFile_B2 ...>] -o Output\n\n");
	fprintf(stdout, "Options: -t INT        number of threads [4]\n");
	fprintf(stdout, "         -fd            files with #1 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stdout, "         -fd2           files with #2 mates reads (format:fa, fq, fq.gz)\n");
	fprintf(stdout, "         -o            alignment filename in SAM format [output.sam]\n");
	fprintf(stdout, "         -p            pair-end mapping\n");
	fprintf(stdout, "\n");
}


bool CheckOutputFileName()
{
	bool bRet = true;
	if (strcmp(outputFile, "output.sam")!=0)
	{

		struct stat s;
		if (stat(outputFile, &s) == 0)
		{
			if (s.st_mode & S_IFDIR)
			{
				bRet = false;
				fprintf(stdout, "Warning: %s is a directory!\n", outputFile);
			}
		}
		int i, len = strlen(outputFile);
		for (i = 0; i < len; i++)
		{
			if (isalnum(outputFile[i]) || outputFile[i] == '/' || outputFile[i] == '.' || outputFile[i] == '-'|| outputFile[i] == '_');
			else
			{
				bRet = false;
				fprintf(stdout, "Warning: [%s] is not a valid filename!\n", outputFile);
				break;
			}
		}		
	}
		

	return bRet;
}

bool CheckInputFiles()
{
	struct stat s;
	bool bRet = true;
	
	
    if ((readFile1 != "") && stat(readFile1.c_str(), &s) == -1)
    {	
        bRet = false;
        fprintf(stdout, "Cannot access file1:[%s]\n", (char*)readFile1.c_str());
    }
    if ((readFile2 != "") &&stat(readFile2.c_str(), &s) == -1)
    {
        bRet = false;
        fprintf(stdout, "Cannot access file2:[%s]\n", (char*)readFile1.c_str());
    }
	return bRet;
}

int main(int argc,char* argv[]){

    int i;
    string parameter;
	threadNum = 1;
	outputFile = (char*)"output.sam";
	FastQFormat = true;
	bPacBioData = false;
		bMultiHit = false;

	bHash =true;
	bBlockRead=true;
	maxInsertSize = 1500;
	MinSeedLength = 0;
	MaxGaps = 5;
    if(argc ==1 || strcmp(argv[1],"-h") ==0){
        Usage();
    }else if (strcmp(argv[1], "index") == 0){
        if(argc ==4){
            bwa_idx_build(argv[2], argv[3]);
        }else{
			fprintf(stderr, "usage: %s index ref.fa prefix\n", argv[0]);
        }

    }
    
    else if (strcmp(argv[1], "align") == 0) {
        for(i=2;i<argc;i++){
            parameter =argv[i];
            if (parameter == "-i"){
                indexFile =argv[++i];
            }else if(parameter == "-fd"){
				
                readFile1 =argv[++i];

            }else if(parameter == "-fd2"){
                readFile2 =argv[++i];
            }else if(parameter == "-t" && i+1 <argc){
                threadNum = atoi(argv[++i]);
                if(threadNum <=0){
                    fprintf(stdout, "Warning! Thread number should be a positive number!\n");
					threadNum = 4;
                }
            }else if (parameter == "-m"){
				bMultiHit = true;
			} 
			//just for test
			else if (parameter == "-h"){
				bHash = false;
			}
			else if (parameter == "-b"){
				bBlockRead = false;
			}
			else if (parameter == "-o"){
				OutputFileFormat = 0;
				outputFile = argv[++i];
				
			}else if(parameter == "-pair" || parameter == "-p"){

				bPairEnd = true;

			}  else if (parameter == "-bo")
			{
				OutputFileFormat = 1;
				outputFile = argv[++i];
			}
			else{
                fprintf(stdout, "Error! Unknown parameter: %s\n", argv[i]);
				Usage();
				exit(1);
            }


        }
        if(readFile1==""){
            fprintf(stdout, "Error! read input can't be empty\n");
            Usage();
            exit(1);
        }
		if (CheckInputFiles() == false || CheckOutputFileName() == false){
            exit(0);
        }
        if (indexFile != NULL && CheckBWAIndexFiles(indexFile)) {
            bwtIdx = bwa_idx_load(indexFile);
        }else{
			fprintf(stdout, "Error! Please specify a valid reference index!\n");
			Usage();
			exit(1);
		}
		if (bwtIdx == 0){
			fprintf(stdout, "\n\nError! Index files are corrupt!\n");
			exit(1);
		}
		else
		{
			bwt = bwtIdx->bwt;
			RestoreReferenceInfo();
			process();	
			bwa_idx_destroy(bwtIdx);
			if (refSeq != NULL) {
               delete[] refSeq; 
            }

		}

    }else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
    return 0;
}