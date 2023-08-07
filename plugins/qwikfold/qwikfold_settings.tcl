# 
# Open a window for changing QwikFold settings
# 
proc ::QWIKFOLD::settings {} {
  variable main_win
  
  # If already initialized, just turn on
  if { [winfo exists $main_win.settings] } {
    wm deiconify $settings_win
    return
  }

  set settings_win [toplevel "$main_win.settings"]
  wm title $settings_win "Settings" 
  wm resizable $settings_win 0 0 
  wm transient $settings_win $main_win
  raise $settings_win

  ########################################################################
  # Path to alphaFold2 "cloned" github 
  ########################################################################
  grid [ ttk::labelframe $settings_win.af -text "AlphaFold path" -relief groove ] -row 0 -padx 5 -pady 5 -sticky nsew

    grid [ttk::label     $settings_win.af.data_label -text "Path" ] -column 0 -row 0 -padx 5 -pady 5 -sticky nsew

    grid [ttk::entry     $settings_win.af.path_entry -state readonly -width 50 -textvariable QWIKFOLD::alphafold_path -validate focus -validatecommand {
            if {[%W get] == "Path to AlphaFold cloned from github"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::alphafold_path "Path to AlphaFold cloned from github"
            }
            return 1
            }] -column 1 -row 0 -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.af.path_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir [pwd] -title "AlphaFold github path"]
        if {$dir != ""} {
            set QWIKFOLD::alphafold_path $dir}
            }] -column 2 -row 0 -padx 5 -pady 5 -sticky e


  ########################################################################
  # Base path to alphaFold2 databases 
  ########################################################################
  grid [ ttk::labelframe $settings_win.dbs -text "AlphaFold Genetic Databases" -relief groove ] -row 1 -padx 5 -pady 5 -sticky nsew

    grid [ttk::label $settings_win.dbs.data_label -text "Base folder" ] -column 0 -row 0 -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.dbs.data_entry -state readonly -width 50 -textvariable QWIKFOLD::alphafold_data -validate focus -validatecommand {
            if {[%W get] == "AlphaFold databases"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::alphafold_data "AlphaFold databases"
            }
            return 1
            }] -column 1 -row 0 -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.dbs.data_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "AlphaFold databases path" ]
        if {$dir != ""} {
            set QWIKFOLD::alphafold_data $dir
            ::QWIKFOLD::set_dbs }
            }] -row 0 -column 2 -padx 5 -pady 5 -sticky nsew


  ########################################################################
  # Non-standard paths to alphaFold2 databases 
  ########################################################################

  grid [ ttk::labelframe $settings_win.full_dbs -text "Full Path to AlphaFold Genetic Databases" -relief groove ] -row 2 -padx 5 -pady 5 -sticky nsew

set dataset_row 0
# BFD
    grid [ttk::label $settings_win.full_dbs.bfd -text "bfd" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.bfd_entry -state readonly -width 50 -textvariable QWIKFOLD::data_bfd -validate focus -validatecommand {
            if {[%W get] == "bfd"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_bfd "bfd database"
            }
            return 1
            }] -column 1 -row 0 -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.bfd_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "bfd database path"]
        if {$dir != ""} {
            set QWIKFOLD::data_bfd $dir}
            }] -row 0 -column 2 -padx 5 -pady 5 -sticky nsew



incr dataset_row 1
# SMALL_BFD
    grid [ttk::label $settings_win.full_dbs.small_bfd -text "small_bfd" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.small_bfd_entry -state readonly -width 50 -textvariable QWIKFOLD::data_small_bfd -validate focus -validatecommand {
            if {[%W get] == "small_bfd"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_small_bfd "small_bfd database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.small_bfd_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "small bfd database path"]
        if {$dir != ""} {
            set QWIKFOLD::data_small_bfd $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# mgnify
    grid [ttk::label $settings_win.full_dbs.mgnify -text "mgnify" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.mgnify_entry -state readonly -width 50 -textvariable QWIKFOLD::data_mgnify -validate focus -validatecommand {
            if {[%W get] == "mgnify"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_mgnify "mgnify database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.mgnify_button -text "Browse" -command {
        set dir [tk_getOpenFile -parent .qwikfold -initialdir ~ -title "mgnify database path" -filetypes [list { {.FASTA files} {.fasta .FASTA} } ] ] 
        if {$dir != ""} {
            set QWIKFOLD::data_mgnify $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# pdb70
    grid [ttk::label $settings_win.full_dbs.pdb70 -text "pdb70" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.pdb70_entry -state readonly -width 50 -textvariable QWIKFOLD::data_pdb70 -validate focus -validatecommand {
            if {[%W get] == "pdb70"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_pdb70 "pdb70 database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.pdb70_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "pdb70 database path"]
        if {$dir != ""} {
            set QWIKFOLD::data_pdb70 $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# pdb_mmcif
    grid [ttk::label $settings_win.full_dbs.pdb_mmcif -text "pdb_mmcif" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.pdb_mmcif_entry -state readonly -width 50 -textvariable QWIKFOLD::data_pdb_mmcif -validate focus -validatecommand {
            if {[%W get] == "pdb_mmcif"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_pdb_mmcif "pdb_mmcif database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.pdb_mmcif_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "pdb_mmcif database path"]
        if {$dir != ""} {
            set QWIKFOLD::data_pdb_mmcif $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew

incr dataset_row 1
# obsolete
    grid [ttk::label $settings_win.full_dbs.obsolete -text "obsolete" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.obsolete_entry -state readonly -width 50 -textvariable QWIKFOLD::data_obsolete -validate focus -validatecommand {
            if {[%W get] == "obsolete"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_obsolete "obsolete database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.obsolete_button -text "Browse" -command {
        set dir [tk_getOpenFile -parent .qwikfold -initialdir ~ -title "obsolete PDB list" -filetypes [list { {.dat files} {.dat .txt} } ] ] 

        if {$dir != ""} {
            set QWIKFOLD::data_obsolete $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# uniclust30
    grid [ttk::label $settings_win.full_dbs.uniclust30 -text "uniclust30" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.uniclust30_entry -state readonly -width 50 -textvariable QWIKFOLD::data_uniclust30 -validate focus -validatecommand {
            if {[%W get] == "uniclust30"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_uniclust30 "uniclust30 database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.uniclust30_button -text "Browse" -command {
        set dir [tk_chooseDirectory -parent .qwikfold -initialdir ~ -title "uniclust30 database path"]
        if {$dir != ""} {
            set QWIKFOLD::data_uniclust30 $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# uniref90
    grid [ttk::label $settings_win.full_dbs.uniref90 -text "uniref90" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.uniref90_entry -state readonly -width 50 -textvariable QWIKFOLD::data_uniref90 -validate focus -validatecommand {
            if {[%W get] == "uniref90"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_uniref90 "uniref90 database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.uniref90_button -text "Browse" -command {
        set dir [tk_getOpenFile -parent .qwikfold -initialdir ~ -title "uniref90.fasta file" -filetypes [list { {.fasta files} {.fasta} } ] ] 

        if {$dir != ""} {
            set QWIKFOLD::data_uniref90 $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# uniprot
    grid [ttk::label $settings_win.full_dbs.uniprot -text "uniprot" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.uniprot_entry -state readonly -width 50 -textvariable QWIKFOLD::data_uniprot -validate focus -validatecommand {
            if {[%W get] == "uniprot"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_uniprot "uniprot database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.uniprot_button -text "Browse" -command {
        set dir [tk_getOpenFile -parent .qwikfold -initialdir ~ -title "uniprot.fasta file" -filetypes [list { {.fasta files} {.fasta} } ] ] 
        if {$dir != ""} {
            set QWIKFOLD::data_uniprot $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


incr dataset_row 1
# pdb_seqres
grid [ttk::label $settings_win.full_dbs.pdb_seqres -text "pdb_seqres" ] -column 0 -row $dataset_row -padx 5 -pady 5 -sticky e

    grid [ttk::entry $settings_win.full_dbs.pdb_seqres_entry -state readonly -width 50 -textvariable QWIKFOLD::data_pdb_seqres -validate focus -validatecommand {
            if {[%W get] == "pdb_seqres"} {
                %W delete 0 end
            } elseif {[%W get] == ""} {
                set QWIKFOLD::data_pdb_seqres "pdb_seqres database"
            }
            return 1
            }] -column 1 -row $dataset_row -padx 5 -pady 5 -sticky nsew

    grid [ttk::button $settings_win.full_dbs.pdb_seqres_button -text "Browse" -command {
        set dir [tk_getOpenFile -parent .qwikfold -initialdir ~ -title "pdb_seqres.txt file" -filetypes [list { {.txt files} {.txt} } ] ] 

        if {$dir != ""} {
            set QWIKFOLD::data_pdb_seqres $dir}
            }] -row $dataset_row -column 2 -padx 5 -pady 5 -sticky nsew


}


proc ::QWIKFOLD::set_dbs {} {
  set QWIKFOLD::data_params     $QWIKFOLD::alphafold_data/params
  set QWIKFOLD::data_bfd        $QWIKFOLD::alphafold_data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
  set QWIKFOLD::data_small_bfd  $QWIKFOLD::alphafold_data/small_bfd/bfd-first_non_consensus_sequences.fasta
  set QWIKFOLD::data_mgnify     $QWIKFOLD::alphafold_data/mgnify/mgy_clusters_2018_12.fa
  set QWIKFOLD::data_pdb70      $QWIKFOLD::alphafold_data/pdb70/pdb70
  set QWIKFOLD::data_pdb_mmcif  $QWIKFOLD::alphafold_data/pdb_mmcif/mmcif_files
  set QWIKFOLD::data_obsolete   $QWIKFOLD::alphafold_data/pdb_mmcif/obsolete.dat
  set QWIKFOLD::data_uniclust30 $QWIKFOLD::alphafold_data/uniclust30/uniclust30_2018_08/uniclust30_2018_08
  set QWIKFOLD::data_uniref90   $QWIKFOLD::alphafold_data/uniref90/uniref90.fasta
  set QWIKFOLD::data_pdb_seqres $QWIKFOLD::alphafold_data/pdb_seqres/pdb_seqres.txt
  set QWIKFOLD::data_uniprot    $QWIKFOLD::alphafold_data/uniprot/uniprot.fasta
  
# In the future we may migrate to a DICTIONARY
#   dict set QWIKFOLD::data_paths data       $QWIKFOLD::alphafold_data
#   dict set QWIKFOLD::data_paths params     $QWIKFOLD::alphafold_data/params
#   dict set QWIKFOLD::data_paths bfd        $QWIKFOLD::alphafold_data/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
#   dict set QWIKFOLD::data_paths small_bfd  $QWIKFOLD::alphafold_data/small_bfd/bfd-first_non_consensus_sequences.fasta
#   dict set QWIKFOLD::data_paths mgnify     $QWIKFOLD::alphafold_data/mgnify/mgy_clusters.fa
#   dict set QWIKFOLD::data_paths pdb70      $QWIKFOLD::alphafold_data/pdb70/pdb70
#   dict set QWIKFOLD::data_paths pdb_mmcif  $QWIKFOLD::alphafold_data/pdb_mmcif/mmcif_files
#   dict set QWIKFOLD::data_paths obsolete   $QWIKFOLD::alphafold_data/pdb_mmcif/obsolete.dat
#   dict set QWIKFOLD::data_paths uniclust30 $QWIKFOLD::alphafold_data/uniclust30/uniclust30_2018_08/uniclust30_2018_08
#   dict set QWIKFOLD::data_paths uniref90   $QWIKFOLD::alphafold_data/uniref90/uniref90.fasta
}