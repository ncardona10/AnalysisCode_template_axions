#Directories
#/disco1/SIMULACIONES/w+jets/w+jets_460/Events/run_01/m_delphes_events.root
#--------------------------------------------------------------
ANALYZERFOLDER="/home/n.cardonac/AnalysisCode_Axions"
PHENOANALYZERFOLDER="/home/n.cardonac/AnalysisCode_Axions/PhenoAnalyzer/"
PHENOANALYZEREXE="PhenoAnalyzer"
INROOTFILE="Events/run_01/m_delphes_events.root"
TEMPORALFOLDER="/home/n.cardonac/AnalysisCode_Axions/RunCode/outPuts/temporal"
OUTPUTFOLDER="/home/n.cardonac/AnalysisCode_Axions/RunCode/outPuts/axions_v5"
#--------------------------------------------------------------
# Processes
#---------------------------------------------------------------
PROCESSFOLDER[1]="/disco4/SIMULACIONES/axions/"
PROCESSSSUBFOLDER[1]="axions"
RUNS[1]=100
TIMES[1]=1

PROCESSFOLDER[2]="/disco1/SIMULACIONES/www"
PROCESSSSUBFOLDER[2]="www"
RUNS[2]=200
TIMES[2]=1

PROCESSFOLDER[3]="/disco1/SIMULACIONES/zzw"
PROCESSSSUBFOLDER[3]="zzw"
RUNS[3]=50
TIMES[3]=1

PROCESSFOLDER[4]="/disco1/SIMULACIONES/zzz"
PROCESSSSUBFOLDER[4]="zzz"
RUNS[4]=200
TIMES[4]=1

PROCESSFOLDER[5]="/disco1/SIMULACIONES/wwz"
PROCESSSSUBFOLDER[5]="wwz"
RUNS[5]=50
TIMES[5]=1

PROCESSFOLDER[6]="/disco2/SIMULACIONES/ttbar"
PROCESSSSUBFOLDER[6]="ttbar"
RUNS[6]=500
TIMES[6]=1

PROCESSFOLDER[7]="/disco1/SIMULACIONES/z+jets"
PROCESSSSUBFOLDER[7]="z+jets"
RUNS[7]=450
TIMES[7]=1

PROCESSFOLDER[8]="/disco1/SIMULACIONES/w+jets"
PROCESSSSUBFOLDER[8]="w+jets"
RUNS[8]=550
TIMES[8]=1

PROCESSFOLDER[9]="/disco3/SIMULACIONES/ww"
PROCESSSSUBFOLDER[9]="ww"
RUNS[9]=250
TIMES[9]=1

PROCESSFOLDER[10]="/disco3/SIMULACIONES/zz"
PROCESSSSUBFOLDER[10]="zz"
RUNS[10]=200
TIMES[10]=1

PROCESSFOLDER[11]="/disco3/SIMULACIONES/wz"
PROCESSSSUBFOLDER[11]="wz"
RUNS[11]=200
TIMES[11]=1

#---------------------------------------------------------------

# index process to run
#-----------------------
INDEX[1]=1
INDEX[2]=1
#----------------------

# Cut value
#--------------
VARIABLE="nocut"
#--------------

typeset -i start_cut=1
typeset -i end_cut=1
typeset -i delta=1
#-------------------------
# loop over the cut values
#------------------------
for ((i = $start_cut; i <= $end_cut; i = i + $delta)); do
	parameters="sed s/mm/$i/ $PHENOANALYZERFOLDER/initial_parameters.in"
	$parameters >$PHENOANALYZERFOLDER/config.in
	cd $PHENOANALYZERFOLDER
	make compile_ROOT_Delphes
	cp *.h $ANALYZERFOLDER
	cp *.in $ANALYZERFOLDER
	cd -
	#-----------------------
	#loop over the processes
	#----------------------
	for j in $(seq ${INDEX[1]} ${INDEX[2]}); do
		typeset -i start_times=1
		typeset -i end_times=${TIMES[$j]} #Number of repetitions
		typeset -i start_runs=1
		typeset -i end_runs=${RUNS[$j]} #Number of runs
		#loop over repetitions
		for k in $(seq $start_times $end_times); do
			#loop over runs
			for l in $(seq $start_runs $end_runs); do
				$PHENOANALYZERFOLDER/$PHENOANALYZEREXE ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_$l.root &
				echo Running ${PROCESSSSUBFOLDER[$j]}_$l.root file
				echo ${PROCESSFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$l/$INROOTFILE >> allProcesses.dat
			done
			wait
			start_runs=start_runs+${RUNS[$j]}
			end_runs=end_runs+${RUNS[$j]}
		done
		hadd -f $OUTPUTFOLDER/${PROCESSSSUBFOLDER[$j]}/${PROCESSSSUBFOLDER[$j]}_$VARIABLE-$i.root $TEMPORALFOLDER/${PROCESSSSUBFOLDER[$j]}_*.root
		cd $TEMPORALFOLDER
		rm *
		cd -
		echo FINISH ${PROCESSSSUBFOLDER[$j]}
	done
done

