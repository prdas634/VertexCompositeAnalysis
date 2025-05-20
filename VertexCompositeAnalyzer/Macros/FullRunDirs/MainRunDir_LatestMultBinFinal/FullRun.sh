ThisDir=$PWD
for iDir in {0..16}
do
    rm -rf $ThisDir/$iDir
    mkdir -p $ThisDir/$iDir
    cp ../NewCodes/GetObs_v2.C $ThisDir/$iDir/.
    cp ../NewCodes/MyConfigAll_v2.h $ThisDir/$iDir/.
    touch $ThisDir/$iDir/RunScript.sh
    #echo "root -l -q 'GetObs_v2.C("$iDir",0,1)'" > $ThisDir/$iDir/RunScript.sh
    echo "root -l -q 'GetObs_v2.C("$iDir",1,1)'" >> $ThisDir/$iDir/RunScript.sh
    #touch $ThisDir/$iDir/RunScript2.sh
    echo "root -l -q 'GetObs_v2.C("$iDir",2,1)'" >> $ThisDir/$iDir/RunScript.sh
    #echo "root -l -q 'GetObs_v2.C("$iDir",3,1)'" >> $ThisDir/$iDir/RunScript.sh
    #echo "root -l -q 'GetObs_v2.C("$iDir",4,1)'" >> $ThisDir/$iDir/RunScript.sh
    #echo "root -l -q 'GetObs_v2.C("$iDir",5,1)'" >> $ThisDir/$iDir/RunScript.sh
    #echo "root -l -q 'GetObs_v2.C("$iDir",6,1)'" >> $ThisDir/$iDir/RunScript.sh
    cd $ThisDir/$iDir
    source RunScript.sh > logRunScript.out &

    cd $ThisDir
done
