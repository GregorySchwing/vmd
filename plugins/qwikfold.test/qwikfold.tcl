#
# $Id: qwikfoldtcl,v 0.5b 2Mon Sep 20 14:40:36 CDT 2021 dgomes $
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

package provide qwikfold 0.5

namespace eval ::QWIKFOLD:: {
	variable topGui ".qwikfold"
	variable bindTop 0

	variable job_id				;	# Job id
	variable alphafold_path		; 	# Path to alphafold installation
	variable conda_path 		;	# Path to af2 environment at your 
								; 	# miniconda3 installation
	variable fasta_sequence    	;	# Fasta sequence  
	variable fasta_path      	;	# Fasta sequence  

#	variable FoldingEngine  ; # Which folding engine to use
#	variable Sequence_ID
#	variable MODELS
}


proc QWIKFOLD::qwikfold {} {
        global env

	set QWIKFOLD::topGui [ toplevel .qwikfold ]
	wm title 		$QWIKFOLD::topGui "QwikFold 0.5b"
	wm resizable 	$QWIKFOLD::topGui 0 0 						; 	#Not resizable

	if {[winfo exists $QWIKFOLD::topGui] != 1} {
			raise $QWIKFOLD::topGui

	} else {
			wm deiconify $QWIKFOLD::topGui
	}


# Source routines
	source $env(QWIKFOLDDIR)/scripts/qwikfold_defaults.tcl		; # Default paths
	source $env(QWIKFOLDDIR)/scripts/qwikfold_notebook.tcl  		; # Main notebook
	source $env(QWIKFOLDDIR)/scripts/qwikfold_functions.tcl		; # Functions

}

#QWIKFOLD::qwikfold

proc qwikfold {} { return [eval QWIKFOLD::qwikfold]}

