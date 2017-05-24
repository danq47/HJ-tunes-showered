if [[ "$@" = "" ]] ; then
    echo "Error - need to specify the directory names you want to create"
    echo "e.g."
    echo "./setup_dirs_for_runs.sh dir1 dir2"
    echo "etc. for however many directories you want"
else
    make -j8 clean
    make -j8 veryclean
    make -j8 pwhg_main
    make -j8 lhef_analysis
    make -j8 main-PYTHIA8_31-lhef
    for DIR in "$@" ; do
	mkdir $DIR
	cd $DIR

	cp ../powheg.input-save .
	cp ../pwhg_main .
	cp ../pwgseeds.dat .
	cp ../lhef_analysis .
	cp ../main31.cmnd .
	cp ../main-PYTHIA8_31-lhef .
	cp ../scripts/event_generation_filtered.sh .
	
	cd ../
	
	echo ""
	echo "Completed " $DIR

    done

    echo ""
    echo "Now don't forget to modify the powheg.input-save files in each directory,"
    echo "plus any modifications you want to make to the shower (main31.cmnd)."
    echo "" 
    echo "Then, running ./event_generation_filtered.sh in the appropriate directory should do the rest."

fi

