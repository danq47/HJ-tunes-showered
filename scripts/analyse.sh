############################################
#                                          #
# Script to analyse the tarred event files #
#                                          #
############################################

do_lhef=1
do_py8=1
do_NLOplots=1
ANALYSIS=""

# Rename all scripts (just in case they aren't in ascending order)

if [ -e events ] ; then
    rm events
fi

if [ -e LHEF-output ] ; then
    rm -rf LHEF-output
fi

if [ -e PYTHIA8-output ] ; then
    rm -rf PYTHIA8-output
fi

if [ -e combined-output ] ; then
    rm -rf combined-output
fi

mkdir combined-output
mkdir LHEF-output
mkdir PYTHIA8-output

ls pwgevents*.gz >> events

numFiles=`ls pwgevents*.gz | wc -l` # this is how many files we have to analyse
numProcessors=$((numFiles/4)) # how many processors we need in total (will split this up into each script doing loops rather than submitting large numbers of scripts)
numScripts=30 # this is how many scripts we want to runnumScripts=$((numFiles/4)) # this is how many scripts we need (each script runs on 4 processors)
bigLoops=$((numProcessors/numScripts))
if [ ! $((numProcessors % numScripts)) -eq 0 ] ; then
    bigLoops=$((bigLoops+1)) # add an extra one to pick up files which don't divide evenly
fi

for i in `seq 1 $numScripts` ; do

    file=`basename $PWD`-analyse-$i.pbs
    cp ../scripts/analyse.pbs $file
    sed -i "s/nScripts=.*/nScripts=$numScripts/g" $file
    sed -i "s/scriptNumber=.*/scriptNumber=$i/g" $file
    sed -i "s/maxBigLoops=.*/maxBigLoops=$bigLoops/g" $file

    if [ $do_lhef -eq 1 ] ; then
        sed -i "s/do_lhef=.*/do_lhef=1/g" $file
    else
        sed -i "s/do_lhef=.*/do_lhef=0/g" $file
    fi

    if [ $do_py8 -eq 1 ] ; then
        sed -i "s/do_py8=.*/do_py8=1/g" $file
    else
        sed -i "s/do_py8=.*/do_py8=0/g" $file
    fi

    if [ $i -eq 1 ] ; then
        ANALYSIS=$(qsub $file)
    else
        ANALYSIS=$ANALYSIS:$(qsub $file)
    fi

done


