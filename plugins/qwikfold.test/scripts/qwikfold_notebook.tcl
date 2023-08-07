########################################################################
# # QwikFold Notebook including tabs (TTK Notebook)
########################################################################
ttk::style layout TNotebook.Tab {Notebook.tab -sticky nswe -children {
        Notebook.padding -expand 1 -sticky nswe -children {Notebook.label
        -expand 1 -sticky nesw }}}

grid [ttk::notebook $QWIKFOLD::topGui.nb ] -row 0 -column 0 -sticky news -pady 5

ttk::frame $QWIKFOLD::topGui.nb.alphafold
ttk::frame $QWIKFOLD::topGui.nb.modeller
$QWIKFOLD::topGui.nb add $QWIKFOLD::topGui.nb.alphafold -text "AlphaFold"
$QWIKFOLD::topGui.nb add $QWIKFOLD::topGui.nb.modeller  -text "Modeller"

########################################################################
# # LabelFrame for alphafold configuration ( .cf )
########################################################################

grid [ ttk::labelframe $QWIKFOLD::topGui.nb.alphafold.cf -text "Configure" -relief groove ] -row 0 -sticky news -pady 5

	########################################################################
	# Job ID
	########################################################################
	grid [ttk::label $QWIKFOLD::topGui.nb.alphafold.cf.id_label -text "Job ID" ] -column 0 -row 0
	grid [ttk::entry $QWIKFOLD::topGui.nb.alphafold.cf.id_entry  -width 40 -textvariable QWIKFOLD::job_id -validate focus -validatecommand {
			if {[%W get] == "myjob"} {
				%W delete 0 end
			} elseif {[%W get] == ""} {
				set QWIKFOLD::job_id "myjob"
			}
			return 1
			}] -column 1 -row 0

	########################################################################
	# Path to alphaFold2 "cloned" github 
	########################################################################
	grid [ttk::label $QWIKFOLD::topGui.nb.alphafold.cf.path_label -text "Path" ] -column 0 -row 1
	grid [ttk::entry $QWIKFOLD::topGui.nb.alphafold.cf.path_entry  -state readonly -width 40 -textvariable QWIKFOLD::alphafold_path -validate focus -validatecommand {
			if {[%W get] == "AlphaFold path"} {
				%W delete 0 end
			} elseif {[%W get] == ""} {
				set QWIKFOLD::alphafold_path "AlphaFold github path"
			}
			return 1
			}] -column 1 -row 1
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.cf.path_button -text "Browse" -command {
		set dir [tk_chooseDirectory -initialdir ~ -title "AlphaFold github path"]
		if {$dir != ""} {
			set QWIKFOLD::alphafold_path $dir}
			}] -row 1 -column 2 -sticky e -padx 5
			
	########################################################################
	# Path to alphaFold2 databases 
	########################################################################
	grid [ttk::label $QWIKFOLD::topGui.nb.alphafold.cf.data_label -text "Databases" ] -column 0 -row 2
	grid [ttk::entry $QWIKFOLD::topGui.nb.alphafold.cf.data_entry -state readonly -width 40 -textvariable QWIKFOLD::alphafold_data -validate focus -validatecommand {
			if {[%W get] == "AlphaFold databases"} {
				%W delete 0 end
			} elseif {[%W get] == ""} {
				set QWIKFOLD::alphafold_data "AlphaFold databases"
			}
			return 1
			}] -column 1 -row 2
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.cf.data_button -text "Browse" -command {
		set dir [tk_chooseDirectory -initialdir ~ -title "AlphaFold databases path"]
		if {$dir != ""} {
			set QWIKFOLD::alphafold_data $dir}
			}] -row 2 -column 2 -sticky e -padx 5


########################################################################
# Path to OUTPUT files
########################################################################
grid [ttk::label $QWIKFOLD::topGui.nb.alphafold.cf.output_label -text "Output folder" ] -column 0 -row 5
grid [ttk::entry $QWIKFOLD::topGui.nb.alphafold.cf.output_entry -state readonly -width 40 -textvariable QWIKFOLD::output_path -validate focus -validatecommand {
		if {[%W get] == "Output folder"} {
			%W delete 0 end
		} elseif {[%W get] == ""} {
			set QWIKFOLD::output_path "Output folder"
		}
		return 1
		}] -column 1 -row 5
grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.cf.output_button -text "Browse" -command {
	set dir [tk_chooseDirectory -initialdir [pwd] -title "Output folder"]
	if {$dir != ""} {
		set QWIKFOLD::output_path $dir}
		}] -row 5 -column 2 -sticky e -padx 5



########################################################################
# LabelFrame for FASTA ( .fasta )
########################################################################

grid [ ttk::labelframe $QWIKFOLD::topGui.nb.alphafold.fasta -text "FASTA sequece" -relief groove ] -column 0 -row 1 -sticky news -pady 5

	########################################################################
	# Text field to input FASTA sequence
	########################################################################
	# Create new frame for Fasta sequence
	grid [ text  $QWIKFOLD::topGui.nb.alphafold.fasta.sequence -width 60 -height 10 -borderwidth 2 -relief sunken -setgrid true ] -column 0 -columnspan 3 -row 1
	#$QWIKFOLD::topGui.nb.alphafold.fasta.sequence insert 1.0 $QWIKFOLD::fasta_sequence

	########################################################################
	# Path to FASTA input file Input file
	########################################################################
	grid [ttk::label $QWIKFOLD::topGui.nb.alphafold.fasta.label -text "FASTA file" ] -column 0 -row 2
	grid [ttk::entry $QWIKFOLD::topGui.nb.alphafold.fasta.entry -width 40 -textvariable QWIKFOLD::fasta -validate focus -validatecommand {
			if {[%W get] == "FASTA file"} {
				%W delete 0 end
			} elseif {[%W get] == ""} {
				set QWIKFOLD::fasta "FASTA file"
			}
			return 1
			}] -column 0 -row 3

# TODO : Prevent overwrite.			
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.fasta.write -text "Write" -command {QWIKFOLD::write_fasta }
		] -row 3 -column 1 -sticky e -padx 5
			
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.fasta.load -text "Load" -command { QWIKFOLD::load_fasta} ] -row 3 -column 2 -sticky e -padx 5




########################################################################
# Submit AlphaFold
########################################################################
grid [ ttk::labelframe $QWIKFOLD::topGui.nb.alphafold.lf -text "Submit & Analyze" -relief groove ] -row 0 -sticky news -pady 5 -row 6

	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.lf.read_button  -text "Read config"  \
		-command {QWIKFOLD::read_config} ]   -row 1 -column 0 -padx 2
	
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.lf.write_button  -text "Write config"  \
		-command {QWIKFOLD::write_alphafold} ] -row 1 -column 1 -padx 2
	
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.lf.submit_button -text "FOLD !" \
		-command {QWIKFOLD::run_alphafold}  ]         -row 1 -column 2 -padx 2
	
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.lf.load_button   -text "Load Models" \
		-command {QWIKFOLD::load_models}  ]  -row 1 -column 3 -padx 2
	
	grid [ttk::button $QWIKFOLD::topGui.nb.alphafold.lf.align_button  -text "Align Models" \
		-command {QWIKFOLD::align_models}  ] -row 1 -column 4 -padx 2



