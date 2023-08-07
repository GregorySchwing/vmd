# Sanity checks


proc check_environment {} {
  # AlphaFold must be properlly installed
  # AND user must load its conda environment before to launching VMD.
  catch {set e [exec python -c "import alphafold"]} result
  if { $result != "" } {
      tk_messageBox -message "Could not load AlphaFold module.\n\nLoad AlphaFold conda environment\nbefore to launching VMD" -icon error -type ok
      break
  }
  
}

proc validate_input_fields {} {
  variable main_win
  # Jobname must not be empty
  if { [info exists QWIKFOLD::job_id]} {
    if { $QWIKFOLD::job_id == "" } {
      tk_messageBox -parent .qwikfold -message "Job name must not be empty" -icon error -type ok
      break
    }
  } else  {
      tk_messageBox -parent .qwikfold -message "Job name is required" -icon error -type ok
    break
  }


# TODO: Make "mkdir" fancier, in case it doesn't work.
  if {[info exists QWIKFOLD::output_path]} {
    if { ! [file isdirectory $QWIKFOLD::output_path]} {
      set answer [tk_messageBox -parent .qwikfold -message "Output folder does not exist.\nMay I create it?" -type yesno -icon question]
      switch -- $answer {
        yes {file mkdir $QWIKFOLD::output_path }
        no break
      }
    }

  } else {
      tk_messageBox -parent .qwikfold -message "Please set the output folder" -icon error -type ok
      break
  }

  # FASTA sequence must not be empty
  set QWIKFOLD::fasta_sequence [ string trim [ $main_win.fasta.sequence get 1.0 end ] ]
  if { $QWIKFOLD::fasta_sequence == "" } {
    tk_messageBox -parent .qwikfold -message "Type or load a FASTA sequence" -icon error -type ok
    break
  }
  
  # Settings for databases
  if { [info exists QWIKFOLD::alphafold_path ] } {
    if { $QWIKFOLD::alphafold_path == "" } {
        tk_messageBox -parent .qwikfold -message "Path to alphafold must be set.\nUse Edit->Settings to assign it." -icon error -type ok
        break
      }
    } else {
        tk_messageBox -parent .qwikfold -message "Path to alphafold must be set.\nUse Edit->Settings to assign it." -icon error -type ok
        break
  }

  if { [info exists QWIKFOLD::alphafold_data ] } {
    if { $QWIKFOLD::alphafold_data == "" } {
      tk_messageBox -parent .qwikfold -message "Path to alphafold databases must be set.\nUse Edit->Settings to assign it." -icon error -type ok
      break
    }
  } else {
      tk_messageBox -parent .qwikfold -message "Path to alphafold databases must be set.\nUse Edit->Settings to assign it." -icon error -type ok
      break
  }

  # Either BFD or Small_bfd must be present.
  if {  $QWIKFOLD::af_mode == "full" || $QWIKFOLD::af_mode == "casp14" } {
    if { $QWIKFOLD::data_bfd == "" } {
    tk_messageBox -parent .qwikfold -message "Path to bfd database must be set.\nUse Edit->Settings to assign it." -icon error -type ok
    break
    }
  }

  if { $QWIKFOLD::af_mode == "reduced" && $QWIKFOLD::data_small_bfd == "" } {
    tk_messageBox -parent .qwikfold -message "Path to small_bfd database must be set.\nUse Edit->Settings to assign it." -icon error -type ok
    break
  }

  # Other parameters are setup based on QWIKFOLD::alphafold_data, enough checking for now.

# Build up model names - No longer available after alphafold 2.01.
  # set model_list {}

  # if {$QWIKFOLD::use_ptm == 1 } {
  #   set ptm "_ptm"
  # } else {
  #   set ptm ""
  # }

  # if { $QWIKFOLD::m1 == 1 } {
  #   lappend model_list "model_1$ptm"
  # }

  # if { $QWIKFOLD::m2 == 1 } {
  #   lappend model_list "model_2$ptm"
  # }
  # if { $QWIKFOLD::m3 == 1 } {
  #   lappend model_list "model_3$ptm"
  # }
  # if { $QWIKFOLD::m4 == 1 } {
  #   lappend model_list "model_4$ptm"
  # }
  # if { $QWIKFOLD::m5 == 1 } {
  #   lappend model_list "model_5$ptm"
  # }

  # if { [llength $model_list] == 0 } {
  #   tk_messageBox -parent .qwikfold -message "Model Parameters missing\nSelect at least one \"model\"" -icon error -type ok
  #   break
  # }
  
  # Save to global variable as comma-separated-values
  #set QWIKFOLD::model_list [join $model_list ","]
puts "to aki"

}


proc print_summary {} {
puts "

# Config
-----------------------------------------------------------
 Job name: $QWIKFOLD::job_id
   Output: $QWIKFOLD::output_path

 Run Mode: $QWIKFOLD::run_mode
Databases: $QWIKFOLD::af_mode
-----------------------------------------------------------

# FASTA Sequence
-----------------------------------------------------------
$QWIKFOLD::fasta_sequence

-----------------------------------------------------------

# Database PATHS
-----------------------------------------------------------
 Alphafold: $QWIKFOLD::alphafold_path
    params: $QWIKFOLD::data_params
       bfd: $QWIKFOLD::data_bfd
 small_bfd: $QWIKFOLD::data_small_bfd
    mgnify: $QWIKFOLD::data_mgnify
     pdb70: $QWIKFOLD::data_pdb70
  obsolete: $QWIKFOLD::data_obsolete
uniclust30: $QWIKFOLD::data_uniclust30
  uniref90: $QWIKFOLD::data_uniref90
-----------------------------------------------------------
"
}



proc check_python_module { module } {
  # Checks if a python module is loadable

  #set module alphafold
  catch {set e [exec python -c "import $module"]} result

  if { $result != "" } {
    puts "[ Error ] $module not found
    conda install $module
    
    Review qwickfold install instructions at 
    http://www.ks.uiuc.edu/Research/vmd/plugins/qwikfold/"

  } else {
      puts "$module found :)"
  }
}
