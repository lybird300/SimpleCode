#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -o timeResult.txt

WORK_DIR=`pwd`
TIME_OUTFILE='ElapsedRealTime.txt'
NUM_REPLICA='50'

for RID in $(eval echo "{1..${NUM_REPLICA}..1}");do 
    { /usr/bin/time -f "%e" java -jar CoJava_noRecomb.jar -p ${WORK_DIR}/params_jeff -o ${WORK_DIR}/testRun_jeff/Rep${RID} >/dev/null ; } 2>> ${TIME_OUTFILE}
    #{ /usr/bin/time -f "%e" java -jar CoJava_noRecomb.jar -p ./params_jeff -o ./testRun_jeff/Rep4 ; } 2>> ElapsedRealTime.txt
done

AVG_TIME=`awk '{a+=$1} END{print a/NR}' ${TIME_OUTFILE}`
MAX_TIME=`sort -n ${TIME_OUTFILE} | head -n1`
MIN_TIME=`sort -n ${TIME_OUTFILE} | tail -n1`

echo "The execution time of ${NUM_REPLICA} replicates (seconds): AVG-${AVG_TIME}, MAX-${MAX_TIME}, MIN-${MIN_TIME}."
exit 0 

#Notes for myself
#(1) The time command gives you real, user, and system time. Use "/usr/bin/time" instead of just "time" allows you to get elapsed real time using the option "-f '%e' "
#(2) 1--stdout (standard output), 2--stderr (standard error)
#(3) The output of "time" command goes to stderr not stdout, but the output of your program (java -jar ... in this case) goes to stdout. To show the time output on the screen (stdout), use
#    { time "your program" ; } 2>&1 
#    To output the result to a file, use
#    { time "your program" ; } 2> output.txt 
#(4) To suppress stdout, direct it to /dev/null
#(5) To append rather than rewrite an existing file, use ">>" instead of ">"
#(6) for RID in $(eval echo "{1..${NUM_REPLICA}..1}") is used instead of "for RID in {1..${NUM_REPLICA}..1}". Otherwise ${RID} will not give the iterating number
#(7) The awk and sort statement lines are used to calculate the mean, max, and min values of a list of numbers, which are stored as consecutive lines of a file (ElapsedRealTime.txt)

#Useful references
#http://stackoverflow.com/questions/556405/what-do-real-user-and-sys-mean-in-the-output-of-time1
#http://www.cyberciti.biz/faq/unix-linux-iterate-over-a-variable-range-of-numbers-in-bash/
#http://www.thegeekstuff.com/2013/10/time-command-format/
#http://stackoverflow.com/questions/617182/with-bash-scripting-how-can-i-suppress-all-output-from-a-command
#http://stackoverflow.com/questions/13356628/is-there-a-way-to-redirect-time-output-to-file-in-linux
#http://stackoverflow.com/questions/9789806/command-line-utility-to-print-statistics-of-numbers-in-linux (the answer from "Skippy le Grand Gourou")