#!/bin/bash

# read arguments

if [[ "$#" -eq 0 ]]
then 
	echo ""
	echo "Usage: 1_Format_SummaryStatistics.sh [-q -l -c -h ] <Work.Dir> <Sum.Stats>"
	echo ""
	echo "<Work.Dir> = Parent directory of the project your summary statistic belongs to"
	echo "<Sum.Stats> = Filename of your summary statistics to be formatted"
	echo ""	
	echo "optional flags: "
	echo "-h = returns this usage"
	echo "-c = add coordinates to summary statistics"
	echo "-q = queue formatting job, more CPUs required than offered in command line"
	echo "-l = liftover from GRCh37 to GRCh38 required"
	echo "NOTE: if -l or -c are used, -q will be applied automatically"
	exit 1
fi


l_flag=0
c_flag=0
q_flag=0
h_flag=0


while getopts 'lcqh' flag; do
	case "${flag}" in
	l) l_flag=1;;
	c) c_flag=1;;
	q) q_flag=1;;
	h) h_flag=1;;
	esac
done

# usage - if asked for using the -h flag
if [[ $h_flag -eq 1 ]]
then 
	echo ""
	echo "Usage: 1_Format_SummaryStatistics.sh [-q -l -c -h ] <Work.Dir> <Sum.Stats>"
	echo ""
	echo "<Work.Dir> = Parent directory of the project your summary statistic belongs to"
	echo "<Sum.Stats> = Filename of your summary statistics to be formatted"
	echo ""	
	echo "optional flags: "
	echo "-h = returns this usage"
	echo "-c = add coordinates to summary statistics"
	echo "-q = queue formatting job, more CPUs required than offered in command line"
	echo "-l = liftover from GRCh37 to GRCh38 required"
	echo "NOTE: if -l or -c are used, -q will be applied automatically"
	exit 1
fi

LIFTOVER=l_flag
ADDCOORDS=c_flag

if [[ $l_flag -eq 1 ]]
then
	q_flag=1
fi

if [[ $c_flag -eq 1 ]]
then
	q_flag=1
fi


WORKDIR=${@:$OPTIND:1}
SS="${WORKDIR}data/summary_statistics/original/${@:$OPTIND+1:1}"
QUEUEJOB=$q_flag
ADDCOORDS=$c_flag
LIFTOVER=$l_flag

OUTDIR="${WORKDIR}data/summary_statistics/formatted/"
TRAIT=$(basename $SS | sed 's/\-.*//')
FILE=${OUTDIR}${TRAIT}.bed


if [[ $QUEUEJOB -eq 1 ]]
then
	echo "job being queued"
	
	echo -n "Give me you trait's abbreviation: "
	read ABBREV 
	
	# writing SHELL script
	echo "TRAIT=$TRAIT" > ${WORKDIR}jobs/Format_SS_$ABBREV.sh
	echo "FILE=$FILE" >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh
		
	if [[ $LIFTOVER -eq 1 ]]
	then
		echo "liftover to be performed (GRCh37 -> GRCh38) "
		echo "/apps/R/4.0.0/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/scripts/1_format_ss_GRCh38.R $SS T F $OUTDIR" >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh

	elif [[ $ADDCOORDS -eq 1 ]]
	then
		echo "genomic coordinates to be added"
		echo "/apps/R/4.0.0/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/scripts/1_format_ss_GRCh38.R $SS F T $OUTDIR" >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh

	else

		echo "/apps/R/4.0.0/INTEL/bin/Rscript /gpfs/scratch/bsc08/bsc08246/method/Gwas_TissueSetEnrich/scripts/1_format_ss_GRCh38.R $SS F F $OUTDIR" >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh

	fi
	
	#sorting 
	echo "sort -k1,1 -k2,2n  \$FILE | /apps/HTSLIB/1.8/INTEL/bin/bgzip -c > \$FILE.bgz " >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh
	#indexing
	echo "/apps/HTSLIB/1.8/INTEL/bin/tabix -p bed \$FILE.bgz" >> ${WORKDIR}jobs/Format_SS_$ABBREV.sh
	
	
	# writing GREASY file
	echo "#!/bin/bash" > ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --job-name='FormatSS_$ABBREV'" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --workdir=${WORKDIR}logs/" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --output=FormatSS_${ABBREV}_%j.out" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --error=FormatSS_${ABBREV}_%j.err" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --cpus-per-task=10" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --ntasks=1" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --qos=bsc_ls" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --time=01:00:00" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "#SBATCH --mail-type=all" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "module load R/4.0.0" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "module load htslib" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	echo "bash ${WORKDIR}jobs/Format_SS_${ABBREV}.sh" >> ${WORKDIR}jobs/Format_SS_$ABBREV.greasy

	sbatch ${WORKDIR}jobs/Format_SS_$ABBREV.greasy
	
else
	# no queuing job, no liftover, no adding coordinates
	echo "job running" 

	# bed 5 format
	/apps/R/4.0.0/INTEL/bin/Rscript 1b_format_ss.R $SS $OUTDIR


	#sort bed and bgzip
	echo "Sorting and Compressing"
	sort -k1,1 -k2,2n  $FILE | /apps/HTSLIB/1.8/INTEL/bin/bgzip -c > $FILE.bgz

	echo "Indexing"
	#index bed for faster reading
	/apps/HTSLIB/1.8/INTEL/bin/tabix -p bed $FILE.bgz	


fi	
