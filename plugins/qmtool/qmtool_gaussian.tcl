# Build route for Gaussian
proc ::QMtool::build_route_Gaussian { args } {
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
  variable route "#"
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
    append route "[lindex $method 1]"
    append route "/$basisset"
  } elseif {[regexp {SE} $method]} {
    # The basisset is not needed in semiempirical methods, its always STO-3G.
    append route "HF"
    append route " [lindex $method 1]"    
  } elseif {[regexp {MP2} $method]} {
    append route "HF"
    append route " $method"
    append route "/$basisset"
  } elseif {[regexp {HF} $method]} {
    append route "HF"
    append route "/$basisset"
  }
    
  if {$geometry=="Checkpoint file"} {
     if {$guess=="Read geometry and wavefunction from checkfile"} {
        if {$coordtype=="ModRedundant" || $coordtype=="Internal (explicit)"} {
           append route " Geom=(AllCheck,ModRedundant)"
        } elseif {$coordtype=="Internal (auto)"} {
           append route " Geom=(AllCheck,NewRedundant)"
        } else {
           append route " Geom=AllCheck"
        }
     } else {
           append route " Geom=Checkpoint"
     }
  } 

  if {$guess=="Take guess from checkpoint file" && $simtype!="Coordinate transformation"} {
     append route " Guess=Read"
  }
  if {$simtype=="Geometry optimization" || $simtype=="Relaxed potential scan" || $simtype=="Rigid potential scan"} {
     variable optmaxcycles
     if {$coordtype=="Internal (auto)"} {
  append route " Opt=(Redundant"
  if {[llength $optmaxcycles]} { append route ",MaxCycle=$optmaxcycles" }
  append route ")"
     } else {
  append route " Opt=(ModRedundant"
  if {[llength $optmaxcycles]} { append route ",MaxCycle=$optmaxcycles" }
  append route ")"
     }
     if {$simtype=="Relaxed potential scan"} {
  lappend popkeys "None"
  append autotitle "<qmtool> simtype=\"Relaxed potential scan\" </qmtool>"
     } elseif {$simtype=="Rigid potential scan"} {
  lappend popkeys "None"
  append autotitle "<qmtool> simtype=\"Rigid potential scan\" </qmtool>"
     } else {
  append autotitle "<qmtool> simtype=\"Geometry optimization\" </qmtool>"
     }
  } elseif {$simtype=="Frequency"} {
     variable hinderedrotor
     if {$coordtype=="Internal (auto)"} {
  if {$hinderedrotor} { 
     append route " Freq=(HinderedRotor)" 
  } else {
     append route " Freq"
  }
     } else {
  if {$hinderedrotor} { 
     append route " Freq=(ModRedundant,HinderedRotor)" 
  } else {
     append route " Freq=(ModRedundant)"
  }
     }
     lappend iops "7/33=1"
     #set otherkey [string trim [regsub -nocase {iop\(7/33=1\)} $otherkey ""]]
     append autotitle "<qmtool> simtype=\"Frequency analysis\" </qmtool>"
  } elseif {$simtype=="Coordinate transformation"} {
     append route " Freq=(Modredundant,ReadFC)"
     lappend iops "7/33=2"
     #set otherkey [string trim [regsub -nocase -all {iop\(7/33=[12]\)} $otherkey ""]]   
     append autotitle "<qmtool> simtype=\"Transformation of force constants from cartesian to internal coordinates\" </qmtool>"
  #} elseif {$simtype=="Rigid potential scan"} {
  #   append route " Scan NoSymm"
  #   append autotitle "<qmtool> simtype=\"Rigid potential scan\" </qmtool>"
  } else {
     append autotitle "<qmtool> simtype=\"Single point calculation\" </qmtool>"
  }

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
     lappend popkeys "NBO"
  }
  if {$calcnboread} {
     lappend popkeys "NBORead"
     if {![llength $havelewis]} {
  tk_messageBox -icon error -type ok -title Message -parent .qm_setup \
     -message "NBO analysis: No Lewis structure defined!\nUse Molefacture to define it." 
  return 0
     }
  }

  variable PCMmethod
  variable solvent
  if {$solvent!="None"} {
     append route " SCRF=($PCMmethod,Solvent=$solvent"
     if {![regexp "SCIPCM|IPCM" $PCMmethod]} {
  append route ",Read"
     }
     append route ")"
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
  
}


########################################
# write input for Gaussian software    #
########################################
proc ::QMtool::write_input_Gaussian { file } {

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
   variable route
   variable title 
   variable coordtype
   variable geometry
   variable atoms
   variable zmat
   variable calcesp
   variable calcresp
   variable calcnbo
   variable calcnboread
   variable solvent
   puts "checkfile2=$checkfile"

   set fid [open $file w]
   puts $fid "%chk=$checkfile"
   puts $fid "%nproc=$nproc"
   puts $fid "%mem=$memory"
   puts $fid $route
   puts $fid ""

   # title and coordinates are not put out if Geom=AllCheck:
   if {!($guess=="Read geometry and wavefunction from checkfile")} {
      set alltitle {}
      if {[llength $autotitle]}  {append alltitle "${autotitle}"}
      if {[llength $extratitle]} {append alltitle "\n${extratitle}"}
      if {[llength $title]}      {append alltitle "\n${title}"}
      
      # title line may not be longer than 80 chars, break it if necessary
      if {[llength $alltitle]} { 
	 set newtitle [break_lines $alltitle 5 79]
	 puts $fid $newtitle
      }
      puts $fid ""
      puts $fid "$totalcharge $multiplicity"
      
      if {$geometry!="Checkpoint file"} {
 	 variable molid
 	 set sel [atomselect $molid all]
 	 foreach coord [$sel get {x y z}] atom [$sel get index] {
 	    # We have to get the name from atomproplist instead of $sel 
 	    # since Paratool might change the names in QMtool.
 	    set name [get_atomprop Name $atom]
	    puts $fid "$name $coord"
 	 }
 	 $sel delete
	 puts $fid ""
      }  
   }

   if {$coordtype=="Internal (explicit)" && [llength $zmat]} {
      # First delete all existing coordinates
      puts $fid "B * * R"
      puts $fid "A * * * R"
      puts $fid "L * * * R"
      puts $fid "D * * * * R"
      puts $fid "O * * * * R"
   }

   if {$coordtype=="ModRedundant" || $coordtype=="Internal (explicit)"} {
      # Fixed cartesian coordinates
      variable atomproplist
      set num 1
      foreach atom $atomproplist {
	 if {[string match {*F*} [get_atomprop Flags [expr {$num-1}]]]} {
	    puts $fid "X $num B"
	    puts $fid "X $num F"
	 }
	 incr num
      }
      # Add internal coordinates explicitly:
      set num 0
      foreach entry $zmat {
	 if {$num==0} { incr num; continue }
	 set indexes {}
	 foreach ind [lindex $entry 2] {
	    lappend indexes [expr {$ind+1}]
	 }
	 set type [string toupper [string index [lindex $entry 1] 0]]
	 # We must model impropers as dihedrals because Gaussian ignores
	 # out-of-plane bends.
	 if {$type=="I"} { set type "D" }
	 if {$type=="O"} { set type "D" }
	 # Something is weird with the linear bend format in Gaussian:
	 if {$type=="L"} { 
	    set val {}
	 } else {
	    set val [lindex $entry 3]
	 }
	 set scan {}
	 if {[string match {*S*} [lindex $entry 5]]} {
	    variable scansteps
	    variable scanstepsize
	    set val "[expr {[lindex $entry 3]-0.5*$scansteps*$scanstepsize}]"
	    set scan "$scansteps $scanstepsize"
	 } else { 
	 }
	 puts $fid "$type $indexes $val [regsub {[QCRM]} [lindex $entry 5] {}]  $scan"
	 incr num
   
      }
      puts $fid ""
   }

   variable PCMmethod
   variable calcdGsolv
   variable solvent
   if {![regexp "SCIPCM|IPCM" $PCMmethod] && $solvent!="None"} {
      puts $fid "RADII=UAHF"
      if {$calcdGsolv} {
	 puts $fid "SCFVAC"
      }
      puts $fid ""
   }

   if {$calcesp} {
      variable molid
      set sel [atomselect $molid all]
      foreach atomicnum [$sel get atomicnumber] {
	 # Merz-Kollman radii in Gaussian are only defined from H - Cl
	 # thus we have to specify the others explicitely
	 if {$atomicnum>17} { 
	    set element [atomnum2element $atomicnum]
	    set radius  [get_esp_radius $element]
	    puts $fid "$element $radius"
	 }
      }
      $sel delete
      puts $fid ""
   }

   if {$calcnboread} {
      write_nbo_input $fid
      puts $fid ""
   }

   close $fid
}
