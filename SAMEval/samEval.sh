
g++ -g -o testSAM_mason  /home/b8402/22_liangjialang/mapping/efficient-iomapping/SAMEval/SamEvaluation_mason.cpp

# dataset=("50" "75" "100" "125" "150" "175" "200")
dataset=("200")

# algo=("eio" "kart" "kartM" "bwamem2" "bwamem2a" "hisat2" "bowtie2")

algo=("eiom" )


for data in "${dataset[@]}"
do

    for i in "${algo[@]}"
    do
        
        /home/b8402/22_liangjialang/mapping/efficient-iomapping/SAMEval/testSAM_mason  /home/b8402/22_liangjialang/mapping/mason_${i}_${data}.sam        -i /home/b8402/22_liangjialang/dataset/short/simulate/e_025/mason2/${data}.sam 


    done

done