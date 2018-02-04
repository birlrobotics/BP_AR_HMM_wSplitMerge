#!/bin/bash 
#  SYNTAX:  runMocap6.sh <infName=SM> <initName=one> <timeLimit=# seconds>
#$ -S /bin/sh 
#$ -cwd 
# ------ attach job number
#$ -j n
# ------ send to particular queue
#$ -l long
# ------ send to particular machine
#$ -q '*@@ang'
# put stdout and stderr files in the right place for your system.
#   NOTE that $TASK_ID is the correct var here
#          but not in rest of script (where SGE_TASK_ID is correct)
#$ -o ../../logs/$JOB_ID.$TASK_ID.out
#$ -e ../../logs/$JOB_ID.$TASK_ID.err

# Edit this line!
PROJHOME=/home/vmrguser/matlab/probability/nonParametricMethods/foxUWash/NPBayesHMM/code/
echo "PROJECT_HOME_PATH" PROJHOME

# ======================================================= SANITIZE INPUT
EXPECTED_ARGS=3
E_BADARGS=65
if [ $# -lt $EXPECTED_ARGS ]  # use -lt instead of <
then
  echo "Usage: `basename $0` <infName=SM> <initName=one> <timeLimit=# seconds>"
  exit $E_BADARGS
fi
infName=$1
initName=$2
TimeLimit=$3

# ======================================================= SET JOB/TASK ID
if [[ ! $SGE_TASK_ID && ${SGE_TASK_ID-_} ]]
then
  # if not deployed on the grid, JOB_ID and SGE_TASK_ID are undefined
  # so we manually set both to make sure this works
  JOB_ID=1
  SGE_TASK_ID=1
fi

echo "JOB ID: " $JOB_ID
echo "SGE TASK ID: "$SGE_TASK_ID

# ======================================================= RUN GRID JOB
# matlabpath=/local/projects/matlab/current/bin/
matlabpath=/home/vmrguser/matlab/R2016b/bin/
export LD_LIBRARY_PATH=$matlabpath;
#MYPREFIX="dbstop if error; cd $PROJHOME; addpath( genpath( '.'));"
MYPREFIX="cd $PROJHOME; addpath( genpath( '.'));"
MYCMD="RunMocap6Experiment $JOB_ID $SGE_TASK_ID $infName $initName $TimeLimit"

echo "EXECUTING CMD:"
echo "              "$MYPREFIX
echo "              "$MYCMD
$matlabpath/matlab -nojvm -nodesktop -nodisplay -nosplash -r "$MYPREFIX; $MYCMD; exit;"
