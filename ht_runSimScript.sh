#!/bin/sh -xv
# Determine the total number of jobs in the queue in  
function queue_total_jobs() {
         typeset -i X=`squeue | grep linly | wc -l`
         echo $X
}

# Set the *maximum* allowed time for a single job to run.  Generally speaking,
# you want this value to be no less than the actual time required for a job
# to run, but not much more than that.  Keeping the value low will help get 
# your jobs into the queue.
# arg 1: e.g., 3:00:00
WALLTIME_PERJOB=$1
# arg 2: e.g., 3000
MAX_TOTAL_JOBS_IN_QUEUE=$2

# Keep track of job ids to determine when we're finished
# $$ represents the process ID number (i.e., PID) of the current shell. For shell scripts, this is the process ID under which they are executing.
JOBIDS=/tmp/jobids_$$

# Move to the right directory
# arg 3: e.g., /projects/sequence_analysis/vol4/linly/workspace/simGWAS/
cd $3
DATA_DIR=`pwd`

# Indicate the common prefix of the batch of .sh scripts to be submitted
# arg 4: e.g., dose2geno_
PREFIX=$4

# Alternative 1
# for JOBFILE in `ls remainingJob[!1]*.sh`; do
#for i in {15..26..1}; do
    #JOBFILE=${DATA_DIR}/chatJob_${i}.sh 
for JOBFILE in `ls ${DATA_DIR}/${PREFIX}.sh`; do
    chmod +x $JOBFILE
    # submit current job	
    # sleep until queue has empty space
    while true; do
        JOBCT=`queue_total_jobs`
        if [ $JOBCT -ge $MAX_TOTAL_JOBS_IN_QUEUE ]; then
            echo "No enough space right now. Sleeping..."
            sleep 100
        else
            break
        fi
    done
    if [ "$JOBCT" -lt "$MAX_TOTAL_JOBS_IN_QUEUE" ]; then
        JOBID=`sbatch $JOBFILE | awk -F \. '{print $DATA_DIR}'`
        #JOBID=`sbatch -t $WALLTIME_PERJOB $JOBFILE | awk -F \. '{print $DATA_DIR}'`
 	  fi
    echo "Submitted job with id $JOBID...."
    echo $JOBID >> $JOBIDS
    sleep 8
done

# Alternative 2
# for i in {0..10..1}; do
#  chmod +x remainingJob${i}.sh
#    echo "Checking current queue count...  `queue_total_jobs` jobs in queue."
#    while [ "`queue_total_jobs`" -gt "$MAX_TOTAL_JOBS_IN_QUEUE" ]; do
#      echo "Sleeping...  `queue_total_jobs` jobs in queue."
#      sleep 2;
#    done;
#	  JOBID=qsub -q largemem ./remainingJob${i}.sh
#	  echo "Submitted job with id $JOBID...."
#	  echo $JOBID >> $JOBIDS
# done;

# let the queue fill
echo "Pausing to let queue fill..."

# sort the list of job ids
#sort -n $JOBIDS > $JOBIDS.sorted
#mv $JOBIDS.sorted $JOBIDS

# spin until all jobs complete
#TF=/tmp/running_$$
#while true;do
#  squeue | grep linly | awk '{print $DATA_DIR}' | sort -n > $TF
#	LEFT=`cat $JOBIDS $TF | sort -n | uniq -d | wc -l`
#	echo "$LEFT jobs left..."
#	if [ $LEFT -eq 0 ]; then
#		break;
#	fi
#	sleep 1000
#done;

#rm $TF
#rm $JOBIDS