#
# $Id: qwikfold.tcl,v 1.1 2021/12/03 23:44:39 johns Exp $
#
#==============================================================================
# QwikFold
#
# Authors:
#   Diego E. B. Gomes
#     Auburn University
#     dgomes@auburn.edu
#
#   Rafael C. Bernardi
#     Beckman Institute for Advanced Science and Technology
#       University of Illinois, Urbana-Champaign
#     Auburn University
#     rcbernardi@ks.uiuc.edu
#     http://www.ks.uiuc.edu/~rcbernardi/
#
# Usage:
#   QwikFold was designed to be used exclusively through its GUI,
#   launched from the "Extensions->Simulation" menu.
#
#   Also see http://www.ks.uiuc.edu/Research/vmd/plugins/qwikfold/ for the
#   accompanying documentation.
#
#=============================================================================

package provide qwikfold 0.7


# REMOVE THIS BEFORE RELEASE ##################################
set env(QWIKFOLDDIR) "/home/dgomes/github/QwikFold/qwikfold0.7"
################################################################

namespace eval ::QWIKFOLD:: {
    namespace export qwikfold
#	variable topGui ".qwikfold"

# 	Window handles
    variable main_win      			;	# handle to main window
    variable settings_win 			;	# handle to settings window
	variable review_win    			;	# handle to review window

#   Job ID will be used as output.
	variable job_id					;	# Job id
	variable af_mode  "reduced_dbs"	;   # AlphaFold Database mode
	variable model_preset "monomer" ;   # monomer, monomer_ptm, monomer_casp14, multimer.
	variable max_template_date "2021-11-01"
#   Where to hold the results
	variable results_folder			;   # Yet to be propagated in more functions

# 	Run options
	variable run_mode "local"		;   # AlphaFold run mode
#	variable alphafold_server		;   # DNS for our server # 	Reserved future use.
	variable use_msa "no"

# 	Original path	
	variable alphafold_path			; 	# Path to alphafold installation
	
# 	FASTA variables
	variable fasta_sequence    		;	# Contents of FASTA sequence
	variable fasta_file      		;	# Path to FASTA sequence file to READ   
	variable fasta_input      		;	# Path to FASTA sequence file to AlphaFold

# 	AlphaFold variables
	variable data_params
	variable data_bfd
	variable data_small_bfd
	variable data_mgnify
	variable data_pdb70
	variable data_pdb_mmcif
	variable data_obsolete
	variable data_uniclust30
	variable data_uniref90
	variable data_pdb_seqres
	variable data_uniprot
	
#   Model parameters ( no longer easly available with alphafold 2.01 )
#   Customization requires modifying the dictionary at "alphafold/models/config.py"
   	#variable m1 1
   	#variable m2 0
	#variable m3 0
   	#variable m4 0
	#variable m5 0
	#variable use_ptm 0	
	#variable gpu_relax 0

	#   Dictionary with model list
    	#variable model_list

# Dictionary
   variable data_paths
}


proc QWIKFOLD::qwikfold {} {
    global env
	variable main_win 

	# AlphaFold must be properlly installed and user must load its conda environment before to launching VMD.
	catch {set e [exec python -c "import alphafold"]} result
	if { $result != "" } {
		tk_messageBox -message "Could not load QwikFold.\n\nPlease load an AlphaFold conda environment before launching VMD" -icon error -type ok
		return
	}
	
	# Main window
	set           main_win [ toplevel .qwikfold ]
	wm title     $main_win "QwikFold 0.7" 
	wm resizable $main_win 0 0     ; #Not resizable

	if {[winfo exists $main_win] != 1} {
			raise $main_win

	} else {
			wm deiconify $main_win
	}

    # Source routines
	source $env(QWIKFOLDDIR)/qwikfold_menubar.tcl   ; # Menu Bar - File
	source $env(QWIKFOLDDIR)/qwikfold_settings.tcl  ; # Menu Bar - Edit
	source $env(QWIKFOLDDIR)/qwikfold_notebook.tcl  ; # Main notebook
	source $env(QWIKFOLDDIR)/qwikfold_functions.tcl	; # Functions
	source $env(QWIKFOLDDIR)/qwikfold_sanity_checks.tcl	; # Sanity checks

}

# Launch the main window
#QWIKFOLD::qwikfold
proc qwikfold {} { return [eval QWIKFOLD::qwikfold]}
