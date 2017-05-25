do_st1=1
do_st2=1
do_NLOplots=1
do_st3=1
do_st4=1
do_tunes_reweight=1
do_LHEF=1 # Can't do analysis on it's own yet, will add this at a later date
do_P8=1   # For now we can only do P8 and LHEF if we also do st4

startingSeed=1
numScripts=20
bigLoops=4
numst2=30

# Stage 1 - grids
if [ $do_st1 -eq 1 ] ; then
# Stage 1 - xgrid 1
    XG1=""
    for i in `seq 1 $numScripts` ; do
        cp ../scripts/st1-1.pbs `basename $PWD`-xg1-$i.pbs
        sed -i "s/nScripts=.*/nScripts=$numScripts/g" `basename $PWD`-xg1-$i.pbs
        sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-xg1-$i.pbs
	sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-xg1-$i.pbs
        if [ $i -eq 1 ] ; then
            XG1=$(qsub -q medium `basename $PWD`-xg1-$i.pbs)
        else
            XG1=$XG1:$(qsub -q medium `basename $PWD`-xg1-$i.pbs)
        fi
    done

# Stage 1.5 - filter out bad files

    XG1F=""
    cp ../scripts/filtering_scripts/filter_xg1.sh .
    cp ../scripts/filtering_scripts/filter_bad_xg1.py .
    XG1F=$(qsub -W depend=afterany:$XG1 filter_xg1.sh)

# Stage 1 - xgrid 2
    XG2=""
    for i in `seq 1 $numScripts` ; do
        cp ../scripts/st1-2.pbs `basename $PWD`-xg2-$i.pbs
        sed -i "s/nScripts=.*/nScripts=$numScripts/g" `basename $PWD`-xg2-$i.pbs
        sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-xg2-$i.pbs
	sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-xg2-$i.pbs
        if [ $i -eq 1 ] ; then
            XG2=$(qsub -q medium -W depend=afterany:$XG1F `basename $PWD`-xg2-$i.pbs)
        else
            XG2=$XG2:$(qsub -q medium -W depend=afterany:$XG1F `basename $PWD`-xg2-$i.pbs)
        fi
    done

# filter out bad xg2 files

    XG2F=""
    cp ../scripts/filtering_scripts/filter_xg2.sh .
    cp ../scripts/filtering_scripts/filter_bad_xg2.py .
    XG2F=$(qsub -W depend=afterany:$XG2 filter_xg2.sh)

fi

# Stage 2 - NLO and btilde upper bound
if [ $do_st2 -eq 1 ] ; then
    ST2=""
    for i in `seq 1 $numst2` ; do
        cp ../scripts/st2.pbs `basename $PWD`-st2-$i.pbs
# Check whether we should do NLO plots
        if [ $do_NLOplots -eq 1 ] ; then
            sed -i "s/do_NLOplots=.*/do_NLOplots=1/g" `basename $PWD`-st2-$i.pbs
        else
            sed -i "s/do_NLOplots=.*/do_NLOplots=0/g" `basename $PWD`-st2-$i.pbs
        fi

        sed -i "s/nScripts=.*/nScripts=$numst2/g" `basename $PWD`-st2-$i.pbs
        sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-st2-$i.pbs
	sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-st2-$i.pbs
# Check whether we depend on st1
        if [ $do_st1 -eq 1 ] ; then
            if [ $i -eq 1 ] ; then
                ST2=$(qsub -q long -W depend=afterany:$XG2F `basename $PWD`-st2-$i.pbs)
            else
                ST2=$ST2:$(qsub -q long -W depend=afterany:$XG2F `basename $PWD`-st2-$i.pbs)
            fi
        else
            if [ $i -eq 1 ] ; then
                ST2=$(qsub -q long `basename $PWD`-st2-$i.pbs)
            else
                ST2=$ST2:$(qsub -q long `basename $PWD`-st2-$i.pbs)
            fi
        fi
    done

# Stage 2.5 - filter out bad st2 files

    ST25=""
    cp ../scripts/filtering_scripts/filter_st2.sh .
    cp ../scripts/filtering_scripts/filter_bad_st2.py .
    ST25=$(qsub -W depend=afterany:$ST2 filter_st2.sh)

fi

# Stage 3 - upper bound for veto algorithm
if [ $do_st3 -eq 1 ] ; then
    ST3=""
    for i in `seq 1 $numScripts` ; do
        cp ../scripts/st3.pbs `basename $PWD`-st3-$i.pbs
        sed -i "s/nScripts=.*/nScripts=$numScripts/g" `basename $PWD`-st3-$i.pbs
        sed -i "s/scriptNumber=1/scriptNumber=$i/g" `basename $PWD`-st3-$i.pbs
	sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-st3-$i.pbs

# Check whether we are dependent on previous stages
        if [ $do_st2 -eq 1 ] ; then
            if [ $i -eq 1 ] ; then
                ST3=$(qsub -q medium -W depend=afterany:$ST25 `basename $PWD`-st3-$i.pbs)
            else
                ST3=$ST3:$(qsub -q medium -W depend=afterany:$ST25 `basename $PWD`-st3-$i.pbs)
            fi
        else
            if [ $i -eq 1 ] ; then
                ST3=$(qsub -q medium `basename $PWD`-st3-$i.pbs)
            else
                ST3=$ST3:$(qsub -q medium `basename $PWD`-st3-$i.pbs)
            fi
        fi
    done
fi

# Stage 4 - event generation
if [ $do_st4 -eq 1 ] ; then
    ST4=""
    for i in `seq 1 $numScripts` ; do
	if [ $do_tunes_reweight -eq 1 ] ; then
	    cp ../scripts/st4-tunes.pbs `basename $PWD`-st4-$i.pbs
	else
	    cp ../scripts/st4.pbs `basename $PWD`-st4-$i.pbs
	fi

# Check if we need to do LHEF/P8 showering
        if [ $do_LHEF -eq 1 ] ; then
            sed -i "s/do_lhef=.*/do_lhef=1/g" `basename $PWD`-st4-$i.pbs
        else
            sed -i "s/do_lhef=.*/do_lhef=0/g" `basename $PWD`-st4-$i.pbs
        fi

        if [ $do_P8 -eq 1 ] ; then
            sed -i "s/do_py8=.*/do_py8=1/g" `basename $PWD`-st4-$i.pbs
        else
            sed -i "s/do_py8=.*/do_py8=0/g" `basename $PWD`-st4-$i.pbs
        fi

        sed -i "s/scriptNumber=.*/scriptNumber=$i/g" `basename $PWD`-st4-$i.pbs
        sed -i "s/nScripts=.*/nScripts=$numScripts/g" `basename $PWD`-st4-$i.pbs
        sed -i "s/scriptNumber-1)\ )\ +\ 1/scriptNumber-1)\ )\ +\ $startingSeed/g" `basename $PWD`-st4-$i.pbs
	sed -i "s/maxBigLoops=.*/maxBigLoops=$bigLoops/g" `basename $PWD`-st4-$i.pbs

# Check whether we depend on previous stages
        if [ $do_st3 -eq 1 ] ; then
	    if [ $i -eq 1 ] ; then
		ST4=$(qsub -q long -W depend=afterany:$ST3 `basename $PWD`-st4-$i.pbs)
	    else
		ST4=$ST4:$(qsub -q long -W depend=afterany:$ST3 `basename $PWD`-st4-$i.pbs)
	    fi
        else
            if [ $i -eq 1 ] ; then
		ST4=$(qsub -q long `basename $PWD`-st4-$i.pbs)
	    else
		ST4=$ST4:$(qsub -q long `basename $PWD`-st4-$i.pbs)
	    fi
        fi
    done
fi


