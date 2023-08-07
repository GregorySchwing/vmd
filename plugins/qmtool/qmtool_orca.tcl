
# Build route for ORCA
proc ::QMtool::build_route_ORCA { args } {
  variable w
  variable basename
  variable checkfile
  variable fromcheckfile
  variable memory
  variable nproc
  variable wavefunction
  variable method
  variable basisset
  variable simtype
  variable otherkey
  variable geometry
  variable guess
  variable autotitle {}
  variable extratitle {}
  variable coordtype
  variable route "!"
  set route2 ""
  variable totalcharge
  variable multiplicity
  variable ncoords
  variable natoms
  variable calcesp
  variable calcnpa
  variable calcnbo
  variable calcnboread
  variable havelewis
  set iops {}
  set popkeys {}
  
  if {[regexp {Closed} $wavefunction]} {
    append route " R"
  } elseif {[regexp {open} $wavefunction]} {
    append route " U"
  } elseif {[regexp {Restricted} $wavefunction]} {
    append route " RO"
  }
  
  if {[regexp {DFT} $method]} {
    append route "KS"
    append route " [lindex $method 1]"
    append route " $basisset"
  } elseif {[regexp {SE} $method]} {
    # The basisset is not needed in semiempirical methods, its always STO-3G.
    append route "HF"
    append route " [lindex $method 1]"    
  } elseif {[regexp {MP2} $method]} {
    append route "HF"
    append route " $method"
    append route " $basisset"
  } elseif {[regexp {HF} $method]} {
    append route "HF"
    append route " $basisset"
  }

  if {$simtype=="Geometry optimization" || $simtype=="Relaxed potential scan" || $simtype=="Rigid potential scan"} {
     variable optmaxcycles
     append route " Opt"
  } elseif {$simtype=="Frequency"} {
     append route " Freq"
     append autotitle "<qmtool> simtype=\"Frequency analysis\" </qmtool>"
  } elseif {$simtype=="Coordinate transformation"} {
  } else {
     append autotitle "<qmtool> simtype=\"Single point calculation\" </qmtool>"
  }

  variable PCMmethod
  variable solvent
  if {$solvent!="None"} {
     set auxsolv $solvent
     if {$solvent=="DiMethylSulfoxide"} {
       set auxsolv "DMSO"
     } elseif {$solvent=="CarbonTetrachloride"} {
       set auxsolv "CCL4"
     } elseif {$solvent=="Dichloromethane"} {   
       set auxsolv "CH2CL2"
     } elseif {$solvent=="TetraHydroFuran"} {   
       set auxsolv "THF"
     }
     append route " CPCM=($auxsolv)"
     if {$PCMmethod=="SMD"} {
       append route2 "\%cpcm smd true \n"
       append route2 "       SMDsolvent \"$auxsolv\" \n"
       append route2 "end"       
     }
  }

  # MARIANO: resume from here
  if {$calcesp} {
     variable molid
     lappend popkeys "ESP" 
     set sel [atomselect $molid "atomicnumber>17"]
     if {[$sel num]} { lappend popkeys "ReadRadii" }
     $sel delete

     # Print the fitting points
     # For paranoid quality also use 6/41=10, 6/42=15" (10 layers, desity 17=^2500 points/atom; 10=^1000points/atom)
     lappend iops    "6/33=2"; 
     append route " NoSymm"
  }
  if {$calcnbo} {
     lappend route "NBO"
  }
  if {$calcnboread} {
     lappend popkeys "NBORead"
     if {![llength $havelewis]} {
  tk_messageBox -icon error -type ok -title Message -parent .qm_setup \
     -message "NBO analysis: No Lewis structure defined!\nUse Molefacture to define it." 
  return 0
     }
  }

  set iops    [join [lsort -unique -ascii $iops] ","]
  set popkeys [join [lsort -unique -ascii $popkeys] ","]
  if {[llength $popkeys]} {
     append route " Pop=($popkeys)"
  }
  if {[llength $iops]} {
     append route " IOp($iops)"
  }
  append route " $otherkey"
  if {$PCMmethod=="SMD"} {
    append route "\n$route2"
  }
}

####################################
# write input for ORCA software    #
####################################
proc ::QMtool::write_input_ORCA { file } {
  variable route
  variable title
  variable autotitle
  variable extratitle
  variable totalcharge
  variable multiplicity
  variable checkfile
  variable nproc
  variable memory
  variable guess
  variable coordtype
  variable geometry
  variable atoms
  variable zmat
  variable calcesp
  variable calcresp
  variable calcnbo
  variable calcnboread
  variable solvent
  
  set fid [open $file w]
  puts $fid $route
  puts $fid ""
  ::ForceFieldToolKit::ORCA::setMemProc $fid $nproc $memory
  puts $fid "\%output"
  puts $fid "  PrintLevel Mini"
  puts $fid "  Print \[ P_Mulliken \] 1"
  puts $fid "  Print \[ P_AtCharges_M \] 1"
  puts $fid "end"
  puts $fid ""
  puts $fid "\%coords"
  puts $fid "  CTyp xyz"
  puts $fid "  Charge $totalcharge"
  puts $fid "  Mult $multiplicity"
  puts $fid "  Units Angs"
  puts $fid "  coords"
  puts $fid ""
  
  for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
      set temp [atomselect top "index $i"]
      lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
      $temp delete
  }
  
  # write the coordinates
  foreach atom_entry $atom_info {
     puts $fid "[lindex $atom_entry 0] [format %16.8f [lindex $atom_entry 1]] [format %16.8f [lindex $atom_entry 2]] [format %16.8f [lindex $atom_entry 3]]"
  }
  # end lines to terminate
  puts $fid "  end"
  puts $fid "end"

  # clean up
  close $fid
    
}