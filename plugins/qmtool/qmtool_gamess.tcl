
# Build route for GAMESS
proc ::QMtool::build_route_GAMESS { args } {
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
  variable route " \$CONTRL SCFTYP="
  variable totalcharge
  variable multiplicity
  variable ncoords
  variable natoms
  variable calcesp
  variable calcnpa
  variable calcnbo
  variable calcnboread
  variable havelewis
  variable solvent
  variable PCMmethod
  set iops {}
  set popkeys {}
  
  if {[regexp {Closed} $wavefunction]} {
    append route "RHF"
  } elseif {[regexp {open} $wavefunction]} {
    append route "UHF"
  } elseif {[regexp {Restricted} $wavefunction]} {
    append route "ROHF"
  }
  
  if {[regexp {DFT} $method]} {
    append route " DFTTYP="
    if {[regexp {LDA} $method]} {
      append route "SVWN"
    } else {
      append route "[lindex $method 1]"
    }
  } elseif {[regexp {MP2} $method]} {
    append route " MPLEVL=2"    
  }
  
  if {$simtype=="Geometry optimization"} {
    append route " RUNTYP=OPTIMIZE"
  } elseif {$simtype=="Frequency"} {
    append route " RUNTYP=HESSIAN"
  }
  
  append route " ICHARG=$totalcharge MULT=$multiplicity"
  append route " \$END\n"

  # Basis group
  append route " \$BASIS"
  if {[regexp {STO} $basisset]} {
    append route " NGAUSS=3 GBASIS=STO"  
  } elseif {[regexp {3-21G} $basisset]} {
    append route " NGAUSS=3 GBASIS=N21"  
  } elseif {[regexp {6-31} $basisset]} {
    append route " NGAUSS=6 GBASIS=N31"  
    if {[regexp {\*} $basisset]} {
      append route " NDFUNC=1"
    }
    if {[regexp {\+} $basisset]} {
      append route " NPFUNC=1"
    }
  } elseif {[regexp {SVP} $basisset] || [regexp {TZVP} $basisset] } {
    append route " GBASIS=$basisset"
   # Special case for SE methods, defined in BASIS    
  } elseif {[regexp {SE} $method]} {
    append route " GBASIS=[lindex $method 1]"
  }
  append route " \$END\n"
  #-------------------------------------#
  
  if {$solvent!="None"} {
    # solvent names are changed except for water, benzene, toluene, aniline and Acetone
    append route " \$PCM "
    if {$solvent=="Methanol"} {
      append route "SOLVNT=CH3OH"
    } elseif {$solvent=="Ethanol"} {
      append route "SOLVNT=C2H5OH"
    } elseif {$solvent=="DiMethylSulfoxide"} {
      append route "SOLVNT=DMSO"
    } elseif {$solvent=="Chloroform"} {
      append route "SOLVNT=CLFORM"
    } elseif {$solvent=="DiChloroMethane"} {
      append route "SOLVNT=METHYCL"
    } elseif {$solvent=="DiChloroEthane"} {
      append route "SOLVNT=12DCLET"
    } elseif {$solvent=="CarbonTetrachloride"} {
      append route "SOLVNT=CTCL"
    } elseif {$solvent=="ChloroBenzene"} {
      append route "SOLVNT=CLBENZ"
    } elseif {$solvent=="NitroMethane"} {
      append route "SOLVNT=NITMET"
    } elseif {$solvent=="Heptane"} {
      append route "SOLVNT=NEPTANE"
    } elseif {$solvent=="CycloHexane"} {
      append route "SOLVNT=CYCHEX"
    } elseif {$solvent=="TetraHydroFuran"} {
      append route "SOLVNT=THF"
    } else {
      append route "SOLVNT=$solvent"
    }
    if {$PCMmethod=="SMD"} {
       append route " SMD=.TRUE."
    }
    append route " \$END\n"
  }
  set words [expr {125*$memory}] ;# In GAMESS, 1 word = 64 bits, 125Mwords = 1Gbyte
  append route " \$SYSTEM MWORDS=$words \$END"
  if {[llength $otherkey]} {
    puts "otherkey"
    append route "\n$otherkey"
  }
  
}

proc ::QMtool::write_input_GAMESS { file } {
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

   for {set i 0} {$i < [molinfo top get numatoms]} {incr i} {
       set temp [atomselect top "index $i"]
       lappend atom_info [list [$temp get element] [$temp get x] [$temp get y] [$temp get z]]
       $temp delete
   }

   set fid [open $file w]
   puts $fid $route
   puts $fid " \$DATA"
   puts $fid "$title"
   puts $fid "C1"
   #puts $fid ""
   # write the coordinates
   foreach atom_entry $atom_info {
      puts $fid "[lindex $atom_entry 0] [format "%3.1f %16.8f %16.8f %16.8f" [lsearch $Molefacture::periodic [lindex $atom_entry 0]] [lindex $atom_entry 1] [lindex $atom_entry 2] [lindex $atom_entry 3]]"
   }   
   puts $fid " \$END"
   
   close $fid

}