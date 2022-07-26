#! /bin/bash

while getopts 'haos' flag; do
    case "${flag}" in
        h) 
        echo "--------------------------------------------------------------"
        echo "./run_HeeP_Analysis.sh -{flags} {variable arguments, see help}"
        echo "--------------------------------------------------------------"
        echo
        echo "The following flags can be called for the heep analysis..."
        echo "    -h, help"
        echo "    -a, analyze"
	echo "        coin -> KIN=arg1"
	echo "        sing -> SPEC=arg1 KIN=arg2 (requires -s flag)"
	echo "    -s, single arm"
	echo "    -o, offset to replay applied"
        exit 0
        ;;
        a) a_flag='true' ;;
        o) o_flag='true' ;;
	s) s_flag='true' ;;
        *) print_usage
        exit 1 ;;
    esac
done

if [[ $a_flag = "true" || $o_flag = "true" || $s_flag = "true" ]]; then
    if [[ $s_flag = "true" ]]; then
	spec=$2
	SPEC=$(echo "$spec" | tr '[:lower:]' '[:upper:]')
	KIN=$3
    else
	KIN=$2
    fi
else
    KIN=$1
fi

# Runs script in the ltsep python package that grabs current path enviroment
if [[ ${HOSTNAME} = *"cdaq"* ]]; then
    PATHFILE_INFO=`python3 /home/cdaq/pionLT-2021/hallc_replay_lt/UTIL_PION/bin/python/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
elif [[ ${HOSTNAME} = *"farm"* ]]; then
    PATHFILE_INFO=`python3 /u/home/${USER}/.local/lib/python3.4/site-packages/ltsep/scripts/getPathDict.py $PWD` # The output of this python script is just a comma separated string
fi

# Split the string we get to individual variables, easier for printing and use later
VOLATILEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f1` # Cut the string on , delimitter, select field (f) 1, set variable to output of command
ANALYSISPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f2`
HCANAPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f3`
REPLAYPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f4`
UTILPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f5`
PACKAGEPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f6`
OUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f7`
ROOTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f8`
REPORTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f9`
CUTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f10`
PARAMPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f11`
SCRIPTPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f12`
ANATYPE=`echo ${PATHFILE_INFO} | cut -d ','  -f13`
USER=`echo ${PATHFILE_INFO} | cut -d ','  -f14`
HOST=`echo ${PATHFILE_INFO} | cut -d ','  -f15`
SIMCPATH=`echo ${PATHFILE_INFO} | cut -d ','  -f16`

if [[ $s_flag = "true" ]]; then
    InDATAFilename="Raw_Data_${SPEC}_${KIN}.root"
    InDUMMYFilename="Raw_DummyData_${SPEC}_${KIN}.root"
    InSIMCFilename="Heep_Coin_${SPEC}_${KIN}.root"
    OutDATAFilename="Analysed_Data_${SPEC}_${KIN}"
    OutDUMMYFilename="Analysed_DummyData_${SPEC}_${KIN}"
    if [[ $o_flag = "true" ]]; then
	OutFullAnalysisFilename="FullAnalysis_Offset_${SPEC}_${KIN}"
    else
	OutFullAnalysisFilename="FullAnalysis_${SPEC}_${KIN}"
    fi
else
    InDATAFilename="Raw_Data_${KIN}.root"
    InDUMMYFilename="Raw_DummyData_${KIN}.root"
    InSIMCFilename="Heep_Coin_${KIN}.root"
    OutDATAFilename="Analysed_Data_${KIN}"
    OutDUMMYFilename="Analysed_DummyData_${KIN}"
    if [[ $o_flag = "true" ]]; then
	OutFullAnalysisFilename="FullAnalysis_Offset_${KIN}"
    else
	OutFullAnalysisFilename="FullAnalysis_${KIN}"
    fi
fi

if [[ $a_flag = "true" ]]; then
    if [[ $s_flag = "true" ]]; then
	cd "${SIMCPATH}/scripts/SING"
	echo
	echo "Analysing ${SPEC} data..."
	echo

	for i in "${data[@]}"
	do
	    echo
	    echo "-----------------------------"
	    echo "Analysing data run $i..."
	    echo "-----------------------------"
	    echo
	    python3 Analysed_SING.py "$i" ${SPEC}
	    #root -l <<EOF 
	    #.x $SIMCPATH/Analysed_SING.C("$InDATAFilename","$OutDATAFilename")
	    #EOF
	done
	cd "${SIMCPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."  
	hadd -f ${OutDATAFilename}.root *_-1_Raw_Data.root
	rm -f *_-1_Raw_Data.root

	cd "${SIMCPATH}/scripts/SING"    
	echo
	echo "Analysing ${SPEC} dummy data..."
	echo

	for i in "${dummydata[@]}"
	do
	    echo
	    echo "-----------------------------------"
	    echo "Analysing dummy data run $i..."
	    echo "-----------------------------------"
	    echo
	    python3 Analysed_SING.py "$i" ${SPEC}
	    #root -l <<EOF 
	    #.x $SIMCPATH/Analysed_SING.C("$InDUMMYFilename","$OutDUMMYFilename")
	    #EOF
	done
	cd "${SIMCPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."
	hadd -f ${OutDUMMYFilename}.root *_-1_Raw_Data.root
	rm -f *_-1_Raw_Data.root	
    else
	cd "${SIMCPATH}/scripts/COIN"
	echo
	echo "Analysing data..."
	echo

	for i in "${data[@]}"
	do
	    echo
	    echo "-----------------------------"
	    echo "Analysing data run $i..."
	    echo "-----------------------------"
	    echo
	    python3 Analysed_COIN.py "$i"
	    #root -l <<EOF 
	    #.x $SIMCPATH/Analysed_COIN.C("$InDATAFilename","$OutDATAFilename")
	    #EOF
	done
	cd "${SIMCPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."  
	hadd -f ${OutDATAFilename}.root *_-1_Raw_Data.root
	rm -f *_-1_Raw_Data.root

	cd "${SIMCPATH}/scripts/COIN"    
	echo
	echo "Analysing dummy data..."
	echo

	for i in "${dummydata[@]}"
	do
	    echo
	    echo "-----------------------------------"
	    echo "Analysing dummy data run $i..."
	    echo "-----------------------------------"
	    echo
	    python3 Analysed_COIN.py "$i"
	    #root -l <<EOF 
	    #.x $SIMCPATH/Analysed_COIN.C("$InDUMMYFilename","$OutDUMMYFilename")
	    #EOF
	done
	cd "${SIMCPATH}/OUTPUT/Analysis/HeeP"
	echo
	echo "Combining root files..."
	hadd -f ${OutDUMMYFilename}.root *_-1_Raw_Data.root
	rm -f *_-1_Raw_Data.root
fi

cd "${SIMCPATH}/scripts"

DataChargeVal=()
DataEffVal=()
echo
echo "Calculating data total charge..."
for i in "${data[@]}"
do
    DataChargeVal+=($(python3 findcharge.py replay_coin_heep "$i" -1))
    DataEffVal+=($(python3 calculate_efficiency.py "$i"))
    #echo "${DataChargeVal[@]} mC"
done
DataChargeSum=$(IFS=+; echo "$((${DataChargeVal[*]}))") # Only works for integers
echo "${DataChargeSum} uC"

DummyChargeVal=()
DummyEffVal=()
echo
echo "Calculating dummy total charge..."
for i in "${dummydata[@]}"
do
    DummyChargeVal+=($(python3 findcharge.py replay_coin_heep "$i" -1))
    DummyEffVal+=($(python3 calculate_efficiency.py "$i"))
    #echo "${DummyChargeVal[@]} mC"
done
DummyChargeSum=$(IFS=+; echo "$((${DummyChargeVal[*]}))") # Only works for integers
echo "${DummyChargeSum} uC"

if [[ $s_flag = "true" ]]; then
    cd "${SIMCPATH}/scripts/SING"
    python3 HeepSing.py ${KIN} "${OutDATAFilename}.root" $DataChargeSum "${DataEffVal[*]}" "${OutDUMMYFilename}.root" $DummyChargeSum "${DummyEffVal[*]}" ${InSIMCFilename} ${OutFullAnalysisFilename} ${SPEC}
else
    cd "${SIMCPATH}/scripts/COIN"
    python3 HeepCoin.py ${KIN} "${OutDATAFilename}.root" $DataChargeSum "${DataEffVal[*]}" "${OutDUMMYFilename}.root" $DummyChargeSum "${DummyEffVal[*]}" ${InSIMCFilename} ${OutFullAnalysisFilename}
fi

cd "${SIMCPATH}"
evince "OUTPUT/Analysis/HeeP/${OutFullAnalysisFilename}.pdf"
