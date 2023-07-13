
#!/bin/bash

#PBS -N CD34_CUX1_ATAC
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -o /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ATAC_seq/CD34_HSC_ATAC/logs/run_ATAC_wrapper.out
#PBS -e /gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ATAC_seq/CD34_HSC_ATAC/logs/run_ATAC_wrapper.err

date
module load gcc/6.2.0



#specify your project directory
project_dir=/gpfs/data/mcnerney-lab/NGS_analysis_tutorials/ATAC_seq/CD34_HSC_ATAC

#change directory to where the fastq files are
cd $project_dir/input

#this for loop will take the input fastq files and run the scripts for all of them one pair after another

for i in $(ls *R1*.gz)
do
otherfilename="${i/R1/R2}"
echo $i
echo $otherfilename


#here you need to specify whether to perform macs2 peak calling by include the -macs2 flag or not. If you include, you need to specify either -p or -q significance threshold followed by a number. Do not specify both p and q values

qsub -v project_path=$project_dir,fq_F=$i,fq_R=$otherfilename,-macs2,-p=0.1 $project_dir/scripts/run_job.sh 
      
done
