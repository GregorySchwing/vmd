########################################################################
# LabelFrame for alphafold configuration ( .cf )
########################################################################
grid [ ttk::labelframe $main_win.cf -text "Configure" -relief groove ] -row 0 -columnspan 2 -padx 5 -pady 5 -sticky news

	########################################################################
	# Job ID
	########################################################################
	grid [ttk::label $main_win.cf.id_label -text "Job name" ] -column 0 -row 0
	grid [ttk::entry $main_win.cf.id_entry  -width 35 -textvariable QWIKFOLD::job_id -validate focus -validatecommand {
			if {[%W get] == "myjob"} {
				%W delete 0 end
			} elseif {[%W get] == ""} {
				set QWIKFOLD::job_id "myjob"
			}
			return 1
			}] -column 1 -columnspan 2 -row 0 -sticky news


	########################################################################
	# Path to OUTPUT files
	########################################################################
	grid [ttk::label $main_win.cf.output_label -text "Results" ] -column 0 -row 1
	 grid [ttk::entry $main_win.cf.output_entry -state readonly -width 35 -textvariable QWIKFOLD::output_path ] -column 1 -row 1
	# grid [ttk::entry $main_win.cf.output_entry -state readonly -width 35 -textvariable QWIKFOLD::output_path -validate focus -validatecommand {
	# 		if {[%W get] == "Output folder"} {
	# 			%W delete 0 end
	# 		} elseif {[%W get] == ""} {
	# 			set QWIKFOLD::output_path "Output folder"
	# 		}
	# 		return 1
	# 		}] -column 1 -row 1
			
	grid [ttk::button $main_win.cf.output_button -text "Browse" -command {
		set dir [tk_chooseDirectory -parent .qwikfold -initialdir [pwd] -title "Output folder"]
		if {$dir != ""} {
			set QWIKFOLD::output_path $dir}
			}] -row 1 -column 2 -sticky news -padx 5



########################################################################
# LabelFrame for run mode configuration ( .model_preset )
########################################################################
grid [ ttk::labelframe $main_win.model_preset -text "Model Preset" -relief groove ] -row 1 -column 0 -padx 5 -pady 5 -sticky news
	set model_preset $main_win.model_preset
	
	# grid [ttk::radiobutton $model_preset.monomer   -text "Monomer"        \
	# 	-variable QWIKFOLD::model_preset -value "monomer" ]        -column 0 -row 0
	
	# grid [ttk::radiobutton $model_preset.monomer_ptm  -text "Monomer_ptm" \
	# 	-variable QWIKFOLD::model_preset -value "monomer_ptm" ]    -column 1 -row 0
	
	# grid [ttk::radiobutton $model_preset.monomer_casp14  -text "Monomer_casp14"0\
	# 	-variable QWIKFOLD::model_preset -value "monomer_casp14" ] -column 0 -row 1

	# grid [ttk::radiobutton $model_preset.multimer  -text "Multimer" \
	# 	-variable QWIKFOLD::model_preset -value "multimer" ]       -column 1 -row 1
	set opt_list [list "monomer" "monomer_ptm" "monomer_casp14" "multimer" ] 
	grid [ ttk::combobox $model_preset.opt -textvariable QWIKFOLD::model_preset -values $opt_list ] -padx 5 -pady 5 -sticky news


########################################################################
# LabelFrame for run mode configuration ( .af_mode )
########################################################################
grid [ ttk::labelframe $main_win.af_mode -text "Database Preset" -relief groove ]  -row 1 -column 1 -padx 5 -pady 5 -sticky news
	set af_mode $main_win.af_mode
	grid [ttk::radiobutton $af_mode.reduced  -text "Reduced" -variable QWIKFOLD::af_mode -value "reduced_dbs" ] -column 0 -row 0 -padx 2 -pady 2 -sticky nsew
	grid [ttk::radiobutton $af_mode.complete -text "Full"    -variable QWIKFOLD::af_mode -value "full_dbs" ]    -column 1 -row 0 -padx 2 -pady 2 -sticky nsew
	#grid [ttk::radiobutton $af_mode.casp14   -text "CASP14"  -variable QWIKFOLD::af_mode -value "casp14" ]  -column 2 -row 0 -padx 2 -pady 2 -sticky nsew

########################################################################
# LabelFrame for FASTA ( .fasta )
########################################################################
grid [ ttk::labelframe $main_win.fasta -text "FASTA sequence" -relief groove ] -row 2 -column 0 -columnspan 2 -padx 5 -pady 5 -sticky news
	# Text field to input FASTA sequence
	grid [ text  $main_win.fasta.sequence -width 60 -height 7 -borderwidth 2 -relief sunken -setgrid true ] -row 1



########################################################################
# Submit AlphaFold
########################################################################
grid [ ttk::frame $main_win.run ] -row 0 -rowspan 3 -column 4 -columnspan 2 -padx 5 -pady 5 -sticky news

	grid [ ttk::labelframe $main_win.run.lf -text "Submit" -relief groove ] -row 0 -columnspan 2 -padx 5 -pady 5 -sticky news
		
		grid [ttk::button $main_win.run.lf.submit_button -text "Run AlphaFold" \
			-command {QWIKFOLD::submit}  ]         -row 1 -column 0 -padx 5 -pady 5 -sticky news

		grid [ttk::checkbutton $main_win.run.lf.check_button -text "Precomputed MSA?" -onvalue "yes" -offvalue "no" \
			-variable QWIKFOLD::use_msa ]  -row 1 -column 1  -padx 5 -pady 5 -sticky news
		
	grid [ ttk::labelframe $main_win.run.results -text "Results" -relief groove ] -row 1 -padx 5 -pady 5 -sticky news

		grid [ttk::button $main_win.run.results.open_folder   -text "Open Folder" \
			-command {QWIKFOLD::open_results}  ]     -row 0 -column 0 -padx 5 -pady 3 -sticky news

		grid [ttk::button $main_win.run.results.load_models   -text "Load Models" \
			-command {QWIKFOLD::load_models}  ]     -row 1 -column 0 -padx 5 -pady 3 -sticky news
		
		grid [ttk::button $main_win.run.results.align_models  -text " Align Models" \
			-command {QWIKFOLD::align_models}  ]    -row 2 -column 0 -padx 5 -pady 3 -sticky news

		grid [ttk::button $main_win.run.results.load_coverage   -text "Alignment Coverage" \
			-command {
				if { $QWIKFOLD::model_preset == "multimer" } {
					QWIKFOLD::load_coverage_multimer
					} else {
						QWIKFOLD::load_coverage
					}
			}  ]    -row 3 -column 0 -padx 5 -pady 3 -sticky news

		grid [ttk::button $main_win.run.results.load_pae   -text "Prediction Error" \
			-command {QWIKFOLD::load_pae}  ]  -row 0 -column 1 -padx 5 -pady 3 -sticky news

		grid [ttk::button $main_win.run.results.load_distogram   -text "Distogram" \
			-command {QWIKFOLD::load_distogram}  ]  -row 1 -column 1 -padx 5 -pady 3 -sticky news

		grid [ttk::button $main_win.run.results.load_contactmap   -text "Contact Map" \
			-command {QWIKFOLD::load_contacts}  ]  -row 2 -column 1 -padx 5 -pady 3 -sticky news	

		grid [ttk::button $main_win.run.results.lddt   -text "lDDT plot" \
			-command {QWIKFOLD::load_lddt}  ]  -row 3 -column 1 -padx 5 -pady 3 -sticky news

	#grid [ ttk::labelframe $main_win.run.view -text "Model(s)" -relief groove ] -row 2 -padx 5 -pady 5 -sticky news





# MSA using MultiSeq is impractical
# Simple molecules such as Ubiquitin render about 30K sequences.
#		grid [ttk::button $main_win.run.an.load_msa   -text "MSA" \
#			-command {QWIKFOLD::load_msa}  ]       -row 1 -column 0 -padx 5 -pady 5 -sticky news
	
#	grid [ ttk::labelframe $main_win.run.ft -text "Features" -relief groove ] -row 0 -column 2 -padx 5 -pady 5 -sticky news
