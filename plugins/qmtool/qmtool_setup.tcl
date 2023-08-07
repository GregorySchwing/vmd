#
# Graphical interfaces for setting up a Gaussian job
#
# $Id: qmtool_setup.tcl,v 1.29 2020/08/13 21:27:02 mariano Exp $
#

#####################################################
# Open window to edit gaussian setup file.          #
#####################################################

proc ::QMtool::setupQM {args} {
   set writemode "ask"
   # Scan for single options
   set argnum 0
   set arglist $args
   foreach i $args {
      if {$i=="-force"}  then {
         set writemode "force"
         set arglist [lreplace $arglist $argnum $argnum]
         continue
      }
      incr argnum
   }
   variable basename
   variable checkfile
   variable fromcheckfile
   variable memory
   variable nproc
   variable wavefunction
   variable method
   variable basisset
   variable simtype
   variable qmsofttype
   variable otherkey
   variable geometry
   variable guess
   variable title
   variable coordtype
   variable totalcharge
   variable multiplicity
   variable availwavefunctions
   variable availmethods
   variable ncoords

   set oldvalues [list $basename $checkfile $fromcheckfile $memory $nproc $method $basisset \
		     $simtype $geometry $guess $title $coordtype $totalcharge $multiplicity]
   if {![llength $basename]}  { 
      if {[molinfo top]>=0} {
	 set basename [file rootname [molinfo top get name]]
      }
   }
   if {![llength $checkfile] && [llength $basename]} { set checkfile "$basename.chk" }
   if {![llength $memory]}    { set memory 1 }
   if {![llength $nproc]}     { set nproc 1 }
   if {![llength $wavefunction]} { set wavefunction "Closed-shell (R)"}
   if {![llength $method] || $method=="HF"} { set method HF }
   if {![llength $basisset]}  { set basisset "6-31G*" }
   if {![llength $simtype]}   { set simtype "Geometry optimization" }
   if {![llength $geometry]}  { set geometry "Z-matrix" }
   if {![llength $guess]}     { set guess "Guess (Harris)" }
   if {![llength $coordtype]} { 
      set coordtype "Internal (explicit)" 
      if {$ncoords==0} {set coordtype "Internal (auto)"}
   }
   if {![llength $title]}     { set title $basename }

   # New version is embedded into frame "args"
   set v $args
   
   trace add variable ::QMtool::wavefunction write ::QMtool::update_route
   trace add variable ::QMtool::method write ::QMtool::update_route
   trace add variable ::QMtool::basisset write ::QMtool::update_route
   trace add variable ::QMtool::simtype write ::QMtool::update_route
   trace add variable ::QMtool::otherkey write ::QMtool::update_route
   trace add variable ::QMtool::qmsofttype write ::QMtool::update_route
   trace add variable ::QMtool::totalcharge write ::QMtool::update_route
   trace add variable ::QMtool::multiplicity write ::QMtool::update_route
   trace add variable ::QMtool::solvent write ::QMtool::update_route
   trace add variable ::QMtool::PCMmethod write ::QMtool::update_route
   
   grid [labelframe $v.edit -text "QM input settings" -bd 2 -padx 8 -pady 8 -width 450] -column 0 -row 0 -padx 7
   
   ttk::style configure TCombobox -background white
   ttk::style map TCombobox -fieldbackground [list readonly #ffffff]

   # NEW QM SELECTION
   grid [frame $v.edit.qmsoft] -column 0 -row 0 -sticky w
   grid [label $v.edit.qmsoft.label -text "QM software:"] -column 0 -row 0 -sticky w
   set availsoft {"Gaussian" "ORCA" "GAMESS"}
   grid [ttk::combobox $v.edit.qmsoft.combQMSoft -values $availsoft -width 12 -state readonly -justify left -textvariable ::QMtool::qmsofttype] -column 1 -row 0 -pady 0 -sticky ns -padx 8

   grid [frame $v.edit.top] -column 0 -row 1 -columnspan 4 -sticky w 
   grid [label $v.edit.top.memorylabel -text "Requested memory (GB):"] -column 0 -row 0 -sticky w
   grid [entry $v.edit.top.memoryentry -textvariable ::QMtool::memory -width 3] -column 1 -row 0 -sticky w
   
   grid [label $v.edit.top.nproclabel -text "Number of processors:"] -column 2 -row 0 -sticky w
   grid [entry $v.edit.top.nprocentry -textvariable ::QMtool::nproc -width 3] -column 3 -row 0 -sticky w
   
   grid [label $v.edit.top.chargelabel -text "Total charge:"] -column 0 -row 1 -sticky w
   grid [entry $v.edit.top.chargeentry -textvariable ::QMtool::totalcharge  -width 3] -column 1 -row 1 -sticky w
   
   grid [label $v.edit.top.multiplabel -text "Multiplicity:"] -column 2 -row 1 -sticky w
   grid [entry $v.edit.top.multipentry -textvariable ::QMtool::multiplicity  -width 3] -column 3 -row 1 -sticky w
   
   grid [frame $v.edit.details] -column 0 -row 2 -sticky w
   grid [label $v.edit.details.wflabel -text "Wavefunction:"] -column 0 -row 0 -sticky w
   grid [ttk::combobox $v.edit.details.combWF -values $availwavefunctions -width 20 -state readonly -justify left -textvariable ::QMtool::wavefunction] -column 1 -row 0 -pady 0 -sticky ns -padx 3
   
   grid [label $v.edit.details.methodlabel -text "Method:"] -column 0 -row 1 -sticky w
   grid [ttk::combobox $v.edit.details.combMethod -values $availmethods -state readonly -width 20 -justify left -textvariable ::QMtool::method] -column 1 -row 1 -pady 0 -sticky ns -padx 3
   
   grid [label $v.edit.details.bslabel -text "Basis set:"] -column 0 -row 2 -sticky w
   set availbasis {"STO-3G" "3-21G" "6-31G" "6-31G*" "6-31+G*" "SVP" "TZVP"}
   grid [ttk::combobox $v.edit.details.combBasis -values $availbasis -state readonly -width 20 -justify left -textvariable ::QMtool::basisset] -column 1 -row 2 -pady 0 -sticky ns -padx 3
   
   grid [label $v.edit.details.typelabel -text "Simulation type:"] -column 0 -row 3 -sticky w
   
   set availsimtype {"Single point" "Geometry optimization" "Frequency" \
	      "Coordinate transformation"}
   grid [ttk::combobox $v.edit.details.combSimType -values $availsimtype -width 20 -state readonly -justify left -textvariable ::QMtool::simtype] -column 1 -row 3 -pady 0 -sticky ns -padx 3
  
   grid [label $v.edit.details.solvlabel -text "Solvent:"] -column 0 -row 4 -sticky w
   set availsolv {"None" "Water" "Methanol" "Ethanol" "DiMethylSulfoxide" \
      "Chloroform" "DiChloroMethane" "DiChloroEthane" "CarbonTetrachloride" "Benzene" "Toluene" \
      "ChloroBenzene" "NitroMethane" "Heptane" "CycloHexane" "Aniline" "Acetone" "TetraHydroFuran"}
   grid [ttk::combobox $v.edit.details.combSolv -values $availsolv -width 20 -state readonly -justify left -textvariable ::QMtool::solvent] -column 1 -row 4 -pady 0 -sticky ns -padx 3

   grid [label $v.edit.details.pcmlabel -text "PCM Method:"] -column 0 -row 5 -sticky w 
   set availpcm {"None" "CPCM" "SMD"}
   grid [ttk::combobox $v.edit.details.combPCM -values $availpcm -width 20 -state readonly -justify left -textvariable ::QMtool::PCMmethod] -column 1 -row 5 -pady 0 -sticky ns -padx 3

   grid [frame $v.edit.nbo] -column 0 -row 13  -sticky w
   grid [label $v.edit.nbo.label -text "Natural Bond Orbitals:"] -column 0 -row 0 -sticky w
   grid [checkbutton $v.edit.nbo.pop -text "NBO population analysis" -variable ::QMtool::calcnbo] -column 1 -row 0 -sticky w
   
   grid [frame $v.edit.text] -column 0 -row 14  -sticky w
   grid [label $v.edit.text.otherlabel -text "Other keywords:"] -column 0 -row 0 -sticky w
   grid [entry $v.edit.text.otherentry -textvariable ::QMtool::otherkey -width 30] -column 1 -row 0 -sticky w;#-width 72
   
   variable autotitle
   variable extratitle
   grid [label $v.edit.text.titlelabel -text "Title string:"] -column 0 -row 1 -sticky nw
   grid [text $v.edit.text.titleentry -width 30 -height 2] -column 1 -row 1 -sticky w
   $v.edit.text.titleentry insert 0.0 $title

   grid [label $v.edit.text.routelabel -text "QM Command\nRoute:"] -column 0 -row 2 -sticky w
   grid [text $v.edit.text.routetext -wrap none -bg white -height 4 -width 30 -font TkFixedFont -relief flat -foreground black \
        -yscrollcommand [list $v.edit.text.scr1 set] -xscrollcommand [list $v.edit.text.scr2 set]] -row 2 -column 1 -sticky wens
   $v.edit.text.routetext insert 0.0 $::QMtool::route

   ##Scroll_BAr V
   scrollbar $v.edit.text.scr1  -orient vertical -command [list $v.edit.text.routetext yview]
   grid $v.edit.text.scr1  -row 2 -column 2  -sticky ens
   ## Scroll_Bar H
   scrollbar $v.edit.text.scr2  -orient horizontal -command [list $v.edit.text.routetext xview]
   grid $v.edit.text.scr2 -row 2 -column 1 -sticky swe


   # Invoke some GUI updates
   if {$simtype=="Coordinate transformation"} { 
    #  .qmtool.fp.input.edit.guessbutton configure -state disabled
    #  .qmtool.fp.input.edit.type.hindrot configure -state disabled
   }
   if {$simtype=="Single point"} { 
    #  .qmtool.fp.input.edit.coordbutton configure -state disabled
  #    .qmtool.fp.input.edit.type.hindrot configure -state disabled
   }
   if {$simtype=="Geometry optimization"} { 
	 #.qmtool.fp.input.edit.type.hindrot configure -state disabled
   }

   #$v.edit.method.entry validate

}

proc ::QMtool::edit_input_cancel { oldvalues } {
   variable basename
   variable checkfile
   variable fromcheckfile
   variable memory
   variable nproc
   variable method
   variable basisset
   variable simtype
   variable otherkey
   variable geometry
   variable guess
   variable title
   variable coordtype
   variable totalcharge
   variable multiplicity

   set basename  [lindex $oldvalues 0]
   set fromcheckfile [lindex $oldvalues 1]
   set checkfile [lindex $oldvalues 2]
   set memory    [lindex $oldvalues 3]
   set nproc     [lindex $oldvalues 4]
   set method    [lindex $oldvalues 5]
   set basisset  [lindex $oldvalues 6]
   set simtype   [lindex $oldvalues 7]
   set otherkey  [lindex $oldvalues 8]
   set geometry  [lindex $oldvalues 9]
   set guess     [lindex $oldvalues 10]
   set title     [lindex $oldvalues 11]
   set coordtype [lindex $oldvalues 12]
   set totalcharge  [lindex $oldvalues 13]
   set multiplicity [lindex $oldvalues 14]
   destroy .qm_setup
}


#NEW update route 
proc ::QMtool::update_route { args } {
   # testing builiding route at each change from scratch
   variable qmsofttype
   variable route
   if {$qmsofttype == "Gaussian" || $qmsofttype == "GAMESS" || $qmsofttype == "ORCA"} {
      build_route_${qmsofttype}   
   } elseif {$qmsofttype == "QCEngine"} {
      set route ""
   } 
   
   .qmtool.fp.input.edit.text.routetext delete 0.0 end
   .qmtool.fp.input.edit.text.routetext insert 0.0 $::QMtool::route
}

proc ::QMtool::edit_input_ok { {action "close"}} {

   variable qmsofttype
   variable w
   variable basename
   variable checkfile
   variable fromcheckfile
   variable memory
   variable nproc
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

   if {$action!="close"} {
      if {$guess!="Read geometry and wavefunction from checkfile" && \
	     (![llength $totalcharge] || ![llength $multiplicity])} {
	 tk_messageBox -icon error -type ok -title Message -parent .qmtool.fp.input \
	    -message "Total charge and multiplicity must be defined!" 
	 focus .qmtool.fp.input.edit.chargeentry
	 return 0
      }
      
   #   if {($simtype=="Coordinate transformation" || $geometry=="Checkpoint file" || \
	      $guess!="Guess (Harris)") && ![llength $fromcheckfile]} {tk_messageBox -icon error -type ok -title Message -parent .qmtool.fp.input \
	    -message "Checkpoint file must be defined!" ;focus .qmtool.fp.input.edit.fcheckfileentry;return 0}

   #   if {($simtype=="Coordinate transformation" || $geometry=="Checkpoint file" || \
	      $guess!="Guess (Harris)") && [llength $fromcheckfile] && ![file exists $fromcheckfile]} {	 tk_messageBox -icon error -type ok -title Message -parent .qmtool.fp.input \
	    -message "Couln't find checkpoint file $fromcheckfile! \nCopy file into working directory manually or choose \"Geometry from: Z-matrix\" and \"Initial wavefunction: Guess\"." ; return 0  }
      
      if {$coordtype=="Internal (explicit)" && $ncoords<3*$natoms-6} {
	 tk_messageBox -icon error -type ok -title Message -parent .qmtool.fp.input \
	    -message "Must have at least 3*natoms-6=[expr {3*$natoms-6}] coordinates.\nCurrently only $ncoords coordinates are defined."
	 return 0
      }
   }

   # saving route that could have been modified by the user
   set customroute [.qmtool.fp.input.edit.text.routetext get 0.0 end]

   update_route "args"

   # check if the new route matches the customroute
   if {$::QMtool::route == $customroute} {
      #do nothing
   } else {
      #replace with customroute???
      set ::QMtool::route $customroute
   }

   #if {$action=="close"} { after idle {destroy .qm_setup} }
   
   if {[llength $fromcheckfile]} {
      append extratitle "<qmtool> parent=$fromcheckfile </qmtool>"
   }

   if {[llength $fromcheckfile]} {
      if {[file exists $basename.chk]} {
	 file rename -force $basename.chk $basename.chk.BAK
      }
      file copy -force $fromcheckfile $basename.chk
      file mtime $basename.chk [clock seconds]
   }
   set checkfile $basename.chk

   set ::QMtool::title [string trimright [.qmtool.fp.input.edit.text.titleentry get 0.0 end]]

   if {$action=="ask"} { 
      opendialog writecom $basename.com
   } elseif {$action=="force"} {
      write_input_$qmsofttype $basename.input
      #destroy .qm_setup
   }
}

################################################################################################
# Merz-Kollman radii (actually Pauling's radii) are only defined from H - Cl and for Br.       #
# Mopac has defined default radii to complete the table trough Bi and uses them whenever       #
# there are no MK-radii. We'll follow the same strategy here and shamelessly steal th radii    #
# from Mopacs manual.                                                                          #
# This is the complete list of radii for all elements from H-Bi.                               #
# The function returns the radius for the given element.                                       #
################################################################################################

proc ::QMtool::get_esp_radius { element } {
   array set radii {H 1.20 He 1.20 Li 1.37 Be 1.45 B 1.45 C 1.50 N 1.50 O 1.40 F 1.35 Ne 1.30 Na 1.57 Mg 1.36 Al 1.24 Si 1.17 P 1.80 S 1.75 Cl 1.70 Ar 1.88 K 2.75 Ca 2.17 Sc 2.26 Ti 2.26 V 2.15 Cr 2.05 Mn 2.10 Fe 2.06 Co 2.05 Ni 1.63 Cu 1.40 Zn 1.39 Ga 1.87 Ge 2.10 As 1.85 Se 1.90 Br 1.80 Kr 2.02 Rb 3.23 Sr 2.94 Y 2.90 Zr 2.85 Nb 2.80 Mo 2.20 Tc 2.20 Ru 2.20 Rh 2.20 Pd 1.63 Ag 1.72 Cd 1.58 In 1.93 Sn 2.17 Sb 2.16 Te 2.06 I 1.98 Xe 2.16 Cs 3.42 Ba 2.97 La 2.40 Hf 2.20 Ta 2.20 W 2.20 Re 2.20 Os 2.20 Ir 2.20 Pt 1.75 Au 1.66 Hg 1.55 Tl 1.96 Pb 2.02 Bi 2.26}

   return [lindex [array get radii $element] 1]
}

proc ::QMtool::write_nbo_input { fid } {
   variable method
   variable zmat
   variable lewislonepairs {}
   variable molid
   array set octet {{} 0 H 2 HE 2 \
			LI 8 BE 8 B 8 C 8 N 8 O 8 F 8 NE 8 NA 8 MG 8 AL 8 SI 8 P 12 S 8 CL 8 AR 8 \
			K 18 CA 18 SC 18 TI 18 V 18 CR 18 MN 18 FE 18 CO 18 NI 18 CU 18 ZN 18 \
			GA 18 GE 18 AS 18 SE 18 BR 18 KR 18 RB 18 SR 18 Y 18 ZR 18 \
			NB 18 MO 18 TC 18 RU 18 RH 18 PD 18 AG 18 CD 18 IN 18 SN 18 SB 18 \
			TE 18  I 18 XE 18 CS 32 BA 0 LA 0 CE 0 PR 0 ND 0 PM 0 SM 0 EU 0 GD 0 TB 0 \
			DY 0 HO 0 ER 0 TM 0 YB 0 LU 0 HF 0 TA 0 W 0 RE 0 OS 0 \
			IR 0 PT 0 AU 0 HG 0 TL 0 PB 0 BI 0 PO 0 AT 0 RN 0 }

   set sel [atomselect $molid all]
   foreach i [$sel list] bonds [$sel getbonds] bos [$sel getbondorders] {
      set lewischarge [get_atomprop Lewis $i]
      set element     [string toupper [get_atomprop Elem  $i]]
      set nbonds 0 
      foreach bo $bos {
	 set bo [expr {int($bo)}]
	 if {$bo < 0} { set bo 1 } 
	 incr nbonds $bo
      }
      set valence [expr {$nbonds-$lewischarge}]
      set numlonepairs [expr {[lindex [array get octet $element] 1]/2-$valence}]
      lappend lewislonepairs $numlonepairs
   }
   $sel delete

   variable molid
   variable molidlist  
   variable molnamelist  
   set filename [lindex $molnamelist [lsearch $molidlist $molid] 1]
   set nbofilename [regsub {(_opt)|(_sp)$} [file rootname $filename] {}]_nbo
   

   puts $fid "\$NBO RESONANCE NLMO PLOT"
   if {[llength $nbofilename]} { puts $fid "     FILE=$nbofilename" }
   puts $fid "\$END"
   puts $fid "\$CHOOSE"
   
   if {[string index $method 0]=="U"} { puts $fid "  ALPHA" }
   
   # Lone pairs
   puts -nonewline $fid "  LONE "
   set index 0
   foreach lp $lewislonepairs {
      puts -nonewline $fid "$index $lp "
      incr index
      if {!($index%8)} { puts -nonewline $fid "\n       " }
   }
   puts $fid " END"
   
   # Bonds
   puts -nonewline $fid "  BOND "
   set i 0
   foreach entry $zmat {
      set type [lindex $entry 1]
      if {![string match "*bond" $type]} { continue }
      set atom0 [lindex $entry 2 0]
      set atom1 [lindex $entry 2 1]
      #set bo [string index $type 0]
      set sel0 [atomselect $molid "index $atom0"]
      set sel1 [atomselect $molid "index $atom1"]
      set pos1in0 [lsearch [join [$sel0 getbonds]] $atom1]
      if {$pos1in0<0} { error "::QMtool::write_nbo_input: Didn't find $atom1 in [$sel0 getbonds]!" }
      set bo0 [lindex [join [$sel0 getbondorders]] $pos1in0]
puts "find $atom1 in [$sel0 getbonds] - $pos1in0 - [$sel0 getbondorders]"
      set pos0in1 [lsearch [join [$sel1 getbonds]] $atom0]
      if {$pos0in1<0} { error "::QMtool::write_nbo_input: Didn't find $atom0 in [$sel1 getbonds]!" }
      set bo1 [lindex [join [$sel1 getbondorders]] $pos0in1]
puts "find $atom0 in [$sel1 getbonds] - $pos0in1 - [$sel1 getbondorders]"
      if {$bo0!=$bo1} { 
	 error "::QMtool::write_nbo_input: Bad bondorder $bo0:$bo1!"
      }
      set bo {}
      if {$bo0==1} { 
	 set bo S
      } elseif {$bo0==2} {
	 set bo D
      } elseif {$bo0==3} {
	 set bo T
      } else { incr i; continue }
      puts -nonewline $fid "$bo $atom0 $atom1 "
      incr i
      if {!($i%4)} { puts -nonewline $fid "\n       " }
   }
   puts $fid " END"
   
   # End of alpha
   if {[string index $method 0]=="U"} { puts $fid "  END" }

#    if {[string index $method 0]=="U"} { 
#       variable lewislonepairsbeta
#       puts $fid "  BETA" 

#       # Lone pairs
#       puts -nonewline $fid "  LONE "
#       set index 0
#       foreach lp $lewislonepairsbeta {
# 	 puts -nonewline $fid "$index $lp "
# 	 incr index
#       }
#       puts $fid " END"
      
#       # Bonds (same as for alpha)
#       puts -nonewline $fid "  BOND "
#       set i 0
#       foreach entry $zmat {
# 	 set type [lindex $entry 1]
# 	 if {![string match "*bond" $type]} { continue }
# 	 set t "S"
# 	 switch $type {
# 	    dbond { set t D }
# 	    dbond { set t T }
# 	 }
# 	 set atom0 [lindex [lindex $entry 2] 0]
# 	 set atom1 [lindex [lindex $entry 2] 1]
# 	 puts -nonewline $fid "$t $atom0 $atom1 "
# 	 incr i
# 	 if {!($i%4)} { puts $fid "" }
#       }
#       puts $fid " END"      
#    }

    puts $fid "\$END"
}
