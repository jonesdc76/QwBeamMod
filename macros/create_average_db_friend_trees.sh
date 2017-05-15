#!/bin/bash
for run in 1 2
do
    for type in 0
    do

	filename="hydrogen_cell_reduced_tree.root"
	treename="reduced_tree"
	nice root -b -q .x MakeAverageDBFriendTree.C+\($run,$type,\"$filename\",\"$treename\"\)
	filename="HYDROGEN-CELL_off_tree.root"
	treename="tree"
	nice root -b -q .x MakeAverageDBFriendTree.C+\($run,$type,\"$filename\",\"$treename\"\)

	filename="hydrogen_cell_corrected_tree.root"
	treename="corrected_tree"
	nice root -b -q .x MakeAverageDBFriendTree.C+\($run,$type,\"$filename\",\"$treename\"\)

	filename="hydrogen_cell_corrected_tree_compton_bpm.root"
	treename="corrected_tree"
	nice root -b -q .x MakeAverageDBFriendTree.C+\($run,$type,\"$filename\",\"$treename\"\)

	filename="hydrogen_cell_corrected_tree_alternate_monitors.root"
	treename="corrected_tree"
	nice root -b -q .x MakeAverageDBFriendTree.C+\($run,$type,\"$filename\",\"$treename\"\)

#	stem="_sine_only_nm"
#	treename="slopes"
#	nice root -b -q .x CreateAverageSlopesTree.C+\($run,\"$stem\",\"$treename\"\)


    done
done
	
