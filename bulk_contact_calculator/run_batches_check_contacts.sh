#!/bin/bash 

set -e;

check_pairs_dir=$1;
minnum=$2;
path_to_check_exe=$3;
path_to_rcsb=$4;
outdir=$5;

mkdir -p $outdir;

files_to_calc=$(mktemp);

wc -l $check_pairs_dir/*/contacts_needcalc.txt | awk -v min=$minnum '{if ($1>min) print $2}' | head -n -1 > $files_to_calc;

biginfile=$(mktemp)


for filename in $(cat $files_to_calc); 
do
    awk -F "_" '{print $1"\t"$2}' $filename >> $biginfile;
done


split $biginfile -l 500000 $outdir/split_;

slurm_job_ids=()
i=0;
for infile in $outdir/split_*;
do
    batchfile=$(mktemp);
    name=$(basename $infile);
    outfile="$outdir/out_$name";
    echo "#!/bin/bash" > $batchfile;
    echo "#SBATCH -A $slurm_project" >> $batchfile;
    echo "#SBATCH -p $slurm_queue" >> $batchfile;
    echo "#SBATCH -t 2-00:00:00" >> $batchfile;
    echo "#SBATCH --mem=4G" >> $batchfile;
    echo "#SBATCH -n 1" >> $batchfile;
    echo "#SBATCH -N 1" >> $batchfile;
    echo "$path_to_check_exe $infile $outfile 8 10 $path_to_rcsb" >> $batchfile;
    sbatch $batchfile;
    #slurm_result=$(sbatch $batchfile);
    #arr=($slurm_result);
    #job_id=${arr[3]};
    #slurm_job_ids=("${slurm_job_ids[@]}" "$job_id")
done


#ids_string=$( IFS=':'; echo "${slurm_job_ids[*]}" );

# wait for jobs
#sbatch --time=0:15:00 -A $slurm_project -p $slurm_queue --wait --dependency afterok:$ids_string /dev/stdin <<< '#!/bin/bash \n sleep 1'


