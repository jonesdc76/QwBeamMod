#!/bin/bash
#Loops over input range of slugs creating files of runlet averages
#tree_type can be "reduced", "corrected", or "set#_reg_reduced"
#tree_name is name of tree in slugfile
#stem is additional file identifier such as "_new_monitors" or "_compton_bpm"
start=$1
end=$2
tree_type=$3
tree_name=$4
stem=$5
dir="/net/data3/paschkelab4/pass5c_reduced_slugfiles"
reduced_dir="/net/data3/paschkelab4/pass5c_reduced_slugfiles"
#for (( islug=${start};islug<=${end};islug+=1 ))

cd /net/data3/paschkelab4/pass5c_reduced_slugfiles/run1
n=0
for islug in $( ls red*| cut -d'g' -f 2 | cut -d'.' -f 1); 
  do
    cd ${BMOD_SRC}/macros
    if [ $islug -ge $start ] && [ $islug -le $end ];then
	echo "$islug"
	((++n))
	nice root -b -q MakeRunletAverages.C+\(${islug},\"${tree_type}\",\"${tree_name}\",\"${stem}\",\"${dir}\",\"${reduced_dir}\"\)
    fi
done

cd /net/data3/paschkelab4/pass5c_reduced_slugfiles/run2
for islug in $( ls red*| cut -d'g' -f 2 | cut -d'.' -f 1); 
  do
    cd ${BMOD_SRC}/macros
    if [ $islug -ge $start ] && [ $islug -le $end ];then
	echo "$islug"
	((++n))
	nice root -b -q MakeRunletAverages.C+\(${islug},\"${tree_type}\",\"${tree_name}\",\"${stem}\",\"${dir}\",\"${reduced_dir}\"\)
    fi
done
echo "$n files"
