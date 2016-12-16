# task: loop through random FASTQ files in the directory and run the script through the files

set -u  

for dir in real-2383_20161024 real-2384_20161024 real-2677_20161024 real-2678_20161024 real-2739_20161024 real-2740_20161024 real-2780_20161024 real-2781_20161024

do
    cd /lustre/jmlab/messme01/vMeyenn/$dir

    if [[ ! -d bam ]]
    then 
        mkdir bam
    fi

    if [[ ! -d log ]]
    then
        mkdir log
    fi

    for x in $(ls log); do echo $x; rm log/$x; done    #removes old log files if there are any


    alignname=Align$RANDOM   	# assign a random number to the alignment variable


    for fastq in $(ls fastq/ | egrep "_R1\.(fastq|fq)")     # pipe all files in fastq/ to egrep and look for all R1-files that contain fastq or fq
    do


        prefix=$(echo $fastq | sed -r "s/_R1\.(fq|fastq)(\.gz)?$//")

        #paired-end sequences:

        mate=$(echo $fastq | sed "s/_R1\.\(fq|fastq\)/_R2\.\1/")

        newbam=bam/${prefix}.bam
        newlog=log/${prefix}.err
        newout=log/${prefix}.out

        # alignment, to index, input fastq variable, use $fastq as first strand and mate as a second, phred +33, unique, RNA data, output in newbam
        cmd="set -e; subread-align -i /lustre/jmlab/resources/genomes/subread/hg38_ERCC -r fastq/$fastq -R fastq/$mate -P 3 -u -t 0 -o $newbam; touch log/$prefix.success"
    
    
        bsub -J $alignname  -R "rusage[mem=16000]" -n 1 -e $newlog -o $newout $cmd 
	
        # batch job -> finally puts the respective file in the lsf system
    done

    cd -

done
