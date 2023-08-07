proc read_config {} {
tk_messageBox -message "Sorry, I can't read config files yet!" -type ok	
}

proc write_fasta {} {
	tk_messageBox -message "Writing FASTA sequence to \n$QWIKFOLD::fasta" -type ok 
	set outfile [ open $QWIKFOLD::fasta w ]
	set sequence [ $QWIKFOLD::topGui.nb.alphafold.fasta.sequence get 1.0 end ]
	puts $outfile $sequence
	close $outfile
}

proc load_fasta {} {
	set dir [tk_getOpenFile -initialdir [pwd] -title "FASTA file"]
	if {$dir != ""} {
		set QWIKFOLD::fasta $dir}
	tk_messageBox -message "Loading FASTA file:\n$QWIKFOLD::fasta" -type ok
	
	#  Slurp up the data file
	set fp [open $QWIKFOLD::fasta r]
	set file_data [read $fp]
	close $fp
	$QWIKFOLD::topGui.nb.alphafold.fasta.sequence insert 1.0 $file_data
}




proc load_models {} {
  set PDB_results [lsort [ glob -tails -directory  $QWIKFOLD::output_path ranked*.pdb ] ]
  set i 0
  foreach e $PDB_results {
    mol new $QWIKFOLD::output_path/$e
    mol modcolor 0 $i ColorID $i
    mol modstyle 0 $i NewCartoon 0.300000 10.000000 4.100000 0
    incr i 1
	}
    display resetview
}


proc  align_models {} {
  set PDB_results [lsort [ glob -tails -directory  $QWIKFOLD::output_path ranked*.pdb ] ]

  set i 0
  set sel0 [atomselect 0 all]

  foreach e $PDB_results {
 	set sel1 [atomselect $i all]	 
	set M [measure fit $sel1 $sel0]	 
	$sel1 move $M
	incr i 1
	}
}


proc write_alphafold {} {

tk_messageBox -message "Writing qwikfold.bash" -type ok 

set outfile [ open "qwikfold.bash" w ]
puts $outfile "\
#!/bin/bash 
eval \"\$(conda shell.bash hook)\"
conda activate af2

python3 run_alphafold_quickfold.py \
	--fasta_paths=$QWIKFOLD::fasta \
	--output_dir=$QWIKFOLD::output_path \
	--model_names='model_1,model_2,model_3,model_4,model_5' \
	--data_dir=$QWIKFOLD::alphafold_data/ \
	--uniref90_database_path=$QWIKFOLD::alphafold_data/uniref90/uniref90.fasta \
	--mgnify_database_path=$QWIKFOLD::alphafold_data/mgnify/mgy_clusters_2018_12.fa \
	--pdb70_database_path=$QWIKFOLD::alphafold_data/pdb70/pdb70 \
	--template_mmcif_dir=$QWIKFOLD::alphafold_data/pdb_mmcif/mmcif_files/ \
	--obsolete_pdbs_path=$QWIKFOLD::alphafold_data/pdb_mmcif/obsolete.dat \
	--preset=reduced_dbs \
	--small_bfd_database_path=$QWIKFOLD::alphafold_data/small_bfd/bfd-first_non_consensus_sequences.fasta \
	--max_template_date=2020-05-14 \
	--jackhmmer_binary_path=\$(which jackhmmer) \
	--hhsearch_binary_path=\$(which hhsearch) \
	--hhblits_binary_path=\$(which hhblits) \
	--kalign_binary_path=\$(which kalign)
    "
close $outfile

}

proc run_alphafold {} {
	set answer [tk_messageBox -message "Running AlphaFold will block VMD until completion\n\nMay I proceed?" -type yesno -icon question]
 switch -- $answer {
   yes {tk_messageBox -message "Running AlphaFold, sit back and relax.\n\[NOTE\] VMD will be unresponsive untill completion." -type ok
	   exec bash qwikfold.bash > $QWIKFOLD::output_path/$QWIKFOLD::job_id.log }
   no {tk_messageBox -message "Ok, we'll fold that latter." -type ok}
 }
	
	}
