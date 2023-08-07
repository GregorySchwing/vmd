proc read_config {} {
tk_messageBox -parent .qwikfold -message "Sorry, I can't read config files yet!" -type ok -icon info
}


proc check_output_folder {} {
  	
	if { ! [info exists QWIKFOLD::output_path] } { 
		tk_messageBox -parent .qwikfold -message "Set the results folder first" -type ok -icon info
		return 1
	} else 	 {
    	if { ! [file isdirectory $QWIKFOLD::output_path]} {
			tk_messageBox -parent .qwikfold -message "Results folder does not exist" -type ok -icon warning
			return 1
		}
	}
}

proc check_job_id {} { 
	if { ! [info exists QWIKFOLD::job_id] } {
		tk_messageBox -parent .qwikfold -message "Job name is not set" -type ok -icon error
		return 1
	} 
}


proc check_msa_folder {  } {
	set msa_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" ]
	if { ! [ file isdirectory $msa_folder ] } {
		tk_messageBox -parent .qwikfold -message "No sequence alignments to load" -type ok -icon info
	return 1
	} 
}



proc open_results {} {
	variable main_win
    global tcl_platform env

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	set QWIKFOLD::results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	if { [ file isdirectory $QWIKFOLD::results_folder ] } {
		switch $tcl_platform(os) {
			"Darwin" {
				eval exec open $QWIKFOLD::results_folder
			}
			"Linux" {
				# xdg-open is supposed to be "distribution agnostic"
				catch [eval exec xdg-open $QWIKFOLD::results_folder] errmsg
			}
			"Windows" {
				eval exec explorer $QWIKFOLD::results_folder
			}
			default {
				tk_messageBox -parent .qwikfold -message "Default file manager unreacheable\Please open $QWIKFOLD::results_folder manually" -type ok -icon error
				return
			}
		}
	} else {
		tk_messageBox -parent .qwikfold -message "Results folder does not exist" -type ok -icon error
		return
	}
}


proc write_fasta {} {
	variable main_win

	set QWIKFOLD::fasta_input [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]
	append QWIKFOLD::fasta_input ".fasta"
#	tk_messageBox -parent .qwikfold -message "Writing FASTA sequence to \n$QWIKFOLD::fasta_input" -type ok -icon info
	set outfile [ open $QWIKFOLD::fasta_input w ]
	set sequence [ $main_win.fasta.sequence get 1.0 end ]
	puts $outfile $sequence
	close $outfile
}

proc load_fasta {} {
	variable main_win

	set dir [tk_getOpenFile -initialdir [pwd] -title "FASTA file" -parent $main_win -filetypes [list { {.FASTA files} {.fasta .FASTA} } ] ]
	if {$dir != ""} {
		set QWIKFOLD::fasta_file $dir}
#	tk_messageBox -parent .qwikfold -message "Loading FASTA file:\n$QWIKFOLD::fasta_file" -type ok -icon info
	
	#  Slurp up the data file
	set fp [open $QWIKFOLD::fasta_file r]
	#set file_data [read $fp]
	set QWIKFOLD::fasta_sequence [read $fp]
	close $fp
	
	# Cleanup and load file contents
	$main_win.fasta.sequence delete 1.0 end
	$main_win.fasta.sequence insert 1.0 $QWIKFOLD::fasta_sequence
}

proc load_project {} {
	variable main_win
	puts "Not implemented"
}


proc load_coverage {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }
	if { [ check_msa_folder ]    == 1 } { return }

	# If already initialized, just turn on
	if { [winfo exists $main_win.coverage] } {
	wm deiconify $coverage_win
	return
	}

	set coverage_win [toplevel "$main_win.coverage"]
	wm title $coverage_win "Alignment" 
	wm resizable $coverage_win 0 0 
	wm transient $coverage_win $main_win
	raise $coverage_win

	# Set alignment folder
	set msa_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" ]

	# Add .sto alignments
	set alignment_list [lsort [ glob -tails -directory  $msa_folder {*.sto} ] ]

	# Appen d .a3m alignments
	lappend alignment_list [lsort [ glob -tails -directory  $msa_folder {*.a3m} ] ]

	grid [ ttk::combobox $coverage_win.mycombo -values $alignment_list ] -padx 5 -pady 5
	# Show 1st in combobox
	#$coverage_win.mycombo set [lindex $alignment_list 0]
	$coverage_win.mycombo configure -state readonly

	# Bind selection. 
	bind $coverage_win.mycombo <<ComboboxSelected>> { set alignment_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" [ %W get ] ] }

	# Plot.
	grid [ttk::button $coverage_win.button -text "Plot alignment\n  coverage" -command { 
		# Set 1st as "alignment file"
		if { ! [info exists alignment_file ] } { 
			tk_messageBox -parent .qwikfold -message "Choose an alignment" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/coverage.py $alignment_file 
		}
		} ] -row 1 -sticky news -padx 5 -pady 5
}

proc load_coverage_multimer {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }
	if { [ check_msa_folder ]    == 1 } { return }

	# If already initialized, just turn on
	if { [winfo exists $main_win.coverage] } {
	wm deiconify $coverage_win
	return
	}

	set coverage_win [toplevel "$main_win.coverage"]
	wm title $coverage_win "Alignment" 
	wm resizable $coverage_win 0 0 
	wm transient $coverage_win $main_win
	raise $coverage_win

	# Set alignment folder
	set msa_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" ]

	# Add .sto alignments
	set sto_list [lsort [ glob -tails -directory  $msa_folder {*/*.sto} ] ]

	# Add .a3m alignments
	set a3m_list [lsort [ glob -tails -directory  $msa_folder {*/*.a3m} ]]

	set alignment_list "$sto_list $a3m_list"

	# pdb_hits.sto not working with python/coverage.py  :(
    set alignment_list [ lsearch -all -inline -not $alignment_list *pdb_hits.sto ]

	grid [ ttk::combobox $coverage_win.mycombo -values $alignment_list ] -padx 5 -pady 5
	# Show 1st in combobox
	#$coverage_win.mycombo set [lindex $alignment_list 0]
	$coverage_win.mycombo configure -state readonly

	# Bind selection. 
	bind $coverage_win.mycombo <<ComboboxSelected>> { set alignment_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" [ %W get ] ] }

	# Plot.
	grid [ttk::button $coverage_win.button -text "Plot alignment\n  coverage" -command { 
		# Set 1st as "alignment file"
		if { ! [info exists alignment_file ] } { 
			tk_messageBox -parent .qwikfold -message "Choose an alignment" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/coverage.py $alignment_file 
		}
		} ] -row 1 -sticky news -padx 5 -pady 5
}



proc load_msa {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	# If already initialized, just turn on
	if { [winfo exists $main_win.msa] } {
	wm deiconify $msa_win
	return
	}

	set msa_win [toplevel "$main_win.msa"]
	wm title $msa_win "Alignment" 
	wm resizable $msa_win 0 0 
	wm transient $msa_win $main_win
	raise $msa_win

	# Set alignment folder
	set msa_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" ]

	# Add .sto alignments
	set alignment_list [lsort [ glob -tails -directory  $msa_folder {*.sto} ] ]

	# Appen d .a3m alignments
	lappend alignment_list [lsort [ glob -tails -directory  $msa_folder {*.a3m} ] ]

	grid [ ttk::combobox $msa_win.mycombo -values $alignment_list ] -padx 5 -pady 5
	# Show 1st in combobox
	#$msa_win.mycombo set [lindex $alignment_list 0]
	$msa_win.mycombo configure -state readonly

	# Bind selection. 
	bind $msa_win.mycombo <<ComboboxSelected>> { set alignment_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id "msas" [ %W get ] ] }

	# Plot.
	grid [ttk::button $msa_win.button -text "Plot alignment\n  coverage" -command { 
		# Set 1st as "alignment file"
		if { ! [info exists alignment_file ] } { 
			tk_messageBox -parent .qwikfold -message "Choose an alignment" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/convert.py $alignment_file
		}
		} ] -row 1 -sticky news -padx 5 -pady 5

}



proc load_distogram {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	# Set results folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	if { [ catch { set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ] } ] } {
		tk_messageBox -parent .qwikfold -message "No models found at $results_folder" -type ok
		return 1
		}

	# If already initialized, just turn on
	if { [winfo exists $main_win.distogram] } {
	wm deiconify $distogram_win
	return
	}

	set distogram_win [toplevel "$main_win.distogram"]
	wm title $distogram_win "Pick a model" 
	wm resizable $distogram_win 0 0 
	wm transient $distogram_win $main_win
	raise $distogram_win

	grid [ ttk::combobox $distogram_win.mycombo -values $model_list ] -padx 5 -pady 5
	#set model_file [lindex $model_list 0]  			; # Set 1st as "model file" # Not working here :(
	#$distogram_win.mycombo set [lindex $model_list 0] 	; # Show 1st in combobox
	$distogram_win.mycombo configure -state readonly

	# Bind selection. 
	bind $distogram_win.mycombo <<ComboboxSelected>> { set model_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id [ %W get ] ] }

	# Plot.
	grid [ttk::button $distogram_win.button -text "Plot distogram" -command { 
	if { ! [info exists model_file ]} { 
			tk_messageBox -parent .qwikfold -message "Pick a model first" -type ok	-icon info 
		} else {
	 
			exec python $env(QWIKFOLDDIR)/python/distogram.py  $model_file } }
	] -row 1 -sticky news -padx 5 -pady 5

}


proc load_contacts {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	# Set results folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	if { [ catch { set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ] } ] } {
		tk_messageBox -parent .qwikfold -message "No models found at $results_folder" -type ok
		return 1
		}


	# If Window already initialized, just turn on
	if { [winfo exists $main_win.contact] } {
		wm deiconify $contact_win
		return
	}

	set contact_win [toplevel "$main_win.contact"]
	wm title $contact_win "Pick a model" 
	wm resizable $contact_win 0 0 
	wm transient $contact_win $main_win
	raise $contact_win


	# Selections
	grid [ ttk::combobox $contact_win.mycombo -values $model_list ] -padx 5 -pady 5
	$contact_win.mycombo configure -state readonly

	# Bind selection. 
	bind $contact_win.mycombo <<ComboboxSelected>> { 
		set model_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id [ %W get ] ]
		}

	# Plot.
	grid [ttk::button $contact_win.button -text "Plot contact" -command { 
		if { ! [info exists model_file ]} { 
			tk_messageBox -parent .qwikfold -message "Pick a model first" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/contact.py  $model_file } }
		] -row 1 -sticky news -padx 5 -pady 5

}



proc load_lddt {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	# Set results folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	if { [ catch { set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ] } ] } {
		tk_messageBox -parent .qwikfold -message "No models found at $results_folder" -type ok
		return 1
		}


	# If already initialized, just turn on
	if { [winfo exists $main_win.lddt] } {
	wm deiconify $lddt_win
	return
	}

	set lddt_win [toplevel "$main_win.lddt"]
	wm title $lddt_win "Pick a model" 
	wm resizable $lddt_win 0 0 
	wm transient $lddt_win $main_win
	raise $lddt_win

	# Set alignment folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ]

	grid [ ttk::combobox $lddt_win.mycombo -values $model_list ] -padx 5 -pady 5
	#set model_file [lindex $model_list 0] 			; # Set 1st as "model file" (not actually working)
	#$lddt_win.mycombo set [lindex $model_list 0] 	; # Show 1st in combobox
	$lddt_win.mycombo configure -state readonly

	# Bind selection. 
	bind $lddt_win.mycombo <<ComboboxSelected>> { 
		set model_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id [ %W get ] ]
		}

	# Plot.
	grid [ttk::button $lddt_win.button -text "Plot lDDT" -command { 
		if { ! [info exists model_file ]} { 
			tk_messageBox -parent .qwikfold -message "Pick a model first" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/lddt.py  $model_file ] } } 
	] -row 1 -sticky news -padx 5 -pady 5

}




proc load_pae {} {
	variable main_win

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return }
	if { [ check_job_id ]        == 1 } { return }

	# Set results folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	if { [ catch { set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ] } ] } {
		tk_messageBox -parent .qwikfold -message "No models found at $results_folder" -type ok
		return 1
		}


	# If already initialized, just turn on
	if { [winfo exists $main_win.pae] } {
	wm deiconify $pae_win
	return
	}

	set pae_win [toplevel "$main_win.pae"]
	wm title $pae_win "Pick a model" 
	wm resizable $pae_win 0 0 
	wm transient $pae_win $main_win
	raise $pae_win

	# Set alignment folder
	set results_folder [file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]

	# List of model .pkl 
	set model_list [lsort [ glob -tails -directory  $results_folder {result*.pkl} ] ]

	grid [ ttk::combobox $pae_win.mycombo -values $model_list ] -padx 5 -pady 5
	#set model_file [lindex $model_list 0] 			; # Set 1st as "model file" (not actually working)
	#$pae_win.mycombo set [lindex $model_list 0] 	; # Show 1st in combobox
	$pae_win.mycombo configure -state readonly

	# Bind selection. 
	bind $pae_win.mycombo <<ComboboxSelected>> { 
		set model_file [file join $QWIKFOLD::output_path $QWIKFOLD::job_id [ %W get ] ]
		}

	# Plot.
	grid [ttk::button $pae_win.button -text "Plot Predicted\nAligned Error (PAE)" -command { 
		if { ! [info exists model_file ]} { 
			tk_messageBox -parent .qwikfold -message "Pick a model first" -type ok	
		} else {
			exec python $env(QWIKFOLDDIR)/python/pae.py  $model_file ] } } 
	] -row 1 -sticky news -padx 5 -pady 5

}



proc load_models {} {

	# Cleanup
	if { [ molinfo list ] > 0 } {
		foreach item [molinfo list] { mol delete $item }
	}

	# Sanity checks 1st
	if { [ check_output_folder ] == 1 } { return 1 }
	if { [ check_job_id ]        == 1 } { return 1 }
	
	set results_folder [ file join $QWIKFOLD::output_path $QWIKFOLD::job_id ]
	puts $results_folder
	puts $QWIKFOLD::job_id

	if { [ file isdirectory $results_folder ] } {
		set PDB_results [lsort [ glob -tails -directory  $results_folder ranked*.pdb ] ]
		set color_id 0
		
		foreach model $PDB_results {
			mol new $results_folder/$model
			set id [molinfo top]
			#mol modcolor 0 $id ColorID $color_id
			mol modcolor 0 $id Beta
			mol modstyle 0 $id NewCartoon 0.300000 10.000000 4.100000 0
			incr color_id 1
			}
		display resetview
	} else {
		tk_messageBox -parent .qwikfold -message "No models at found $results_folder" -type ok
		return
	}

	tk_messageBox -parent .qwikfold -message "Models colored by prediction quality" -type ok
	return
}

proc align_models {} {
	
	if { [ molinfo num ] == 0 } {
		tk_messageBox -parent .qwikfold -message "Nothing to alingn" -type ok
		return
	}

	# Get molecule list
	set mol_list [molinfo list]

	# Set first of list as TOP
    mol top [ lindex $mol_list 0 ]

	# Select all atoms from TOP structure
	set sel0 [atomselect top all]	

	foreach i $mol_list {
		set sel1 [atomselect $i all]
		set M [measure fit $sel1 $sel0]	 
		$sel1 move $M
	}

} 

proc run_alphafold {} {
	if { $QWIKFOLD::af_mode == "reduced_dbs" } { 
		set alphafold_cmd "python3 [ file join $QWIKFOLD::alphafold_path run_alphafold.py ] \
		--fasta_paths=$QWIKFOLD::fasta_input \
		--output_dir=$QWIKFOLD::output_path \
		--data_dir=$QWIKFOLD::alphafold_data/ \
		--uniref90_database_path=$QWIKFOLD::data_uniref90 \
		--mgnify_database_path=$QWIKFOLD::data_mgnify \
		--template_mmcif_dir=$QWIKFOLD::data_pdb_mmcif \
		--obsolete_pdbs_path=$QWIKFOLD::data_obsolete \
		--model_preset=$QWIKFOLD::model_preset \
		--db_preset=$QWIKFOLD::af_mode \
		--small_bfd_database_path=$QWIKFOLD::data_small_bfd \
		--jackhmmer_binary_path=jackhmmer \
		--hhsearch_binary_path=hhsearch \
		--kalign_binary_path=kalign \
		--max_template_date=$QWIKFOLD::max_template_date"

	} else {
		set alphafold_cmd "python3 [ file join $QWIKFOLD::alphafold_path run_alphafold.py ] \
		--fasta_paths=$QWIKFOLD::fasta_input \
		--output_dir=$QWIKFOLD::output_path \
		--data_dir=$QWIKFOLD::alphafold_data/ \
		--uniref90_database_path=$QWIKFOLD::data_uniref90 \
		--uniclust30_database_path=$QWIKFOLD::data_uniclust30 \
		--mgnify_database_path=$QWIKFOLD::data_mgnify \
		--template_mmcif_dir=$QWIKFOLD::data_pdb_mmcif \
		--obsolete_pdbs_path=$QWIKFOLD::data_obsolete \
		--model_preset=$QWIKFOLD::model_preset \
		--db_preset=$QWIKFOLD::af_mode \
		--bfd_database_path=$QWIKFOLD::data_bfd \
		--jackhmmer_binary_path=jackhmmer \
		--hhsearch_binary_path=hhsearch \
		--hhblits_binary_path=hhblits \
		--kalign_binary_path=kalign \
		--max_template_date=$QWIKFOLD::max_template_date"
	}

	if { $QWIKFOLD::model_preset == "multimer" } {
		append alphafold_cmd " \
		--pdb_seqres_database_path=$QWIKFOLD::data_pdb_seqres \
        --uniprot_database_path=$QWIKFOLD::data_uniprot"

	} else {
		append alphafold_cmd " \
		--pdb70_database_path=$QWIKFOLD::data_pdb70"
	}

	if { $QWIKFOLD::use_msa == "yes" } {
		append alphafold_cmd "\
		--use_precomputed_msas"
	}
	
	puts $alphafold_cmd

	# TEST
	#return

	if [ catch { eval exec $alphafold_cmd } alphafold_log ] {
		puts $alphafold_log
	}

}

proc submit {} {
	
	# Double check if alphafold is intalled/loadable
	check_environment 

	# Check required input fields
	if { [ catch validate_input_fields code ] } { return 1 }

	# Print summary
	print_summary

	# Write the FASTA file to "output" dir
	write_fasta

	# Actually run alphafold
	run_alphafold

	return
}
