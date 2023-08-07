#
# Analysis routines
#
# $Id: qmtool_analysis.tcl,v 1.13 2021/11/19 16:58:40 mariano Exp $
#

#######################################
# NEW: Create Analysis window and add #
# QMEnergiesFrame from QwikMD         # 
#######################################
proc ::QMtool::createAnalysisTab { v } {
  
  #provisory call to parse_molinfo1
  #::QMtool::parse_molinfo1 $::QMtool::molid
  set type_selection ""
 
  #grid [frame $v.fp ] -row 0 -column 0 -sticky nsew -pady 2 -padx 2
  ttk::style configure TNotebook -tabposition nw
  ttk::style configure TNotebook.Tab -width 15
  ttk::style configure TNotebook.Tab -anchor left
  grid [ttk::notebook $v.fp -style TNotebook -width 300] -row 0 -column 0 -sticky nsew -pady 2 -padx 2
  grid columnconfigure $v.fp 0 -weight 1
  #grid rowconfigure $v.fp 0 -weight 0
  grid rowconfigure $v.fp 1 -weight 1
  #grid rowconfigure $v.fp 2 -weight 2
  
  frame $v.fp.energies
  frame $v.fp.nma 
  frame $v.fp.orbitals 

  #image create photo imag1 -format gif -file [file join $::env(QMTOOLDIR) energies.gif]
  #image create photo imag2 -format gif -file [file join $::env(QMTOOLDIR) frequencies.gif]
  #image create photo imag3 -format gif -file [file join $::env(QMTOOLDIR) orbitals.gif]
  #$v.fp add $v.fp.energies -image imag1 -padding 5 -sticky news
  #$v.fp add $v.fp.nma -image imag2 -padding 5 -sticky news
  #$v.fp add $v.fp.orbitals -image imag3 -padding 5 -sticky news

  $v.fp add $v.fp.energies -text "Plot Energies" -padding 5 -sticky news
  $v.fp add $v.fp.nma -text "Normal Modes" -padding 5 -sticky news
  $v.fp add $v.fp.orbitals -text "Molecular Orbitals" -padding 5 -sticky news

  #normalmode_gui $v.fp.nma ; #testing button option to create table
   
  grid [frame $v.fp.energies.plot] -row 0 -column 0 -sticky news
  grid [button $v.fp.energies.plot.button -text "Plot Energies" -pady 2 -padx 2 -width 15 -command {
    ::QMtool::plot_scf_energies .qmtool.fp.analysis.fp.energies.plot
    }] -row 0 -column 0 -sticky e -pady 2 -padx 2
  #grid [frame $v.fp.energies.plot.plot] -row 1 -column 0 -sticky news
  
  grid [frame $v.fp.orbitals.header] -row 0 -column 0 -sticky news -pady 2 -padx 2 
  grid [button $v.fp.orbitals.header.button -text "MO Table" -pady 2 -padx 2 -width 15 -command {
    ::QMtool::call_qm_orbitals .qmtool.fp.analysis.fp.orbitals.header
    }] -row 0 -column 0 -sticky e -pady 2 -padx 2
  grid columnconfigure $v.fp.orbitals.header 0 -weight 1
  grid [button $v.fp.orbitals.header.reset -text "Reset Table" -pady 2 -padx 2 -width 15 -command {
    if {[winfo exist .qmtool.fp.analysis.fp.orbitals.header.tableresults]} {
       destroy .qmtool.fp.analysis.fp.orbitals.header.tableresults
    }
    #::QMtool::call_qm_orbitals $v.fp.orbitals.header
    }] -row 0 -column 1 -sticky e -pady 2 -padx 2
  grid columnconfigure $v.fp.orbitals.header 0 -weight 1
  
 
  grid [frame $v.fp.nma.header] -row 0 -column 0 -sticky news -pady 2 -padx 2 
  grid [button $v.fp.nma.header.button -text "Frequencies Table" -pady 2 -padx 2 -width 15 -command {
    ::QMtool::normalmode_gui .qmtool.fp.analysis.fp.nma.header 
    }] -row 0 -column 0 -sticky e -pady 2 -padx 2
  grid [button $v.fp.nma.header.reset -text "Reset Table" -pady 2 -padx 2 -width 15 -command {
   if {[winfo exist .qmtool.fp.analysis.fp.nma.header.mode]} {
      destroy .qmtool.fp.analysis.fp.nma.header.mode
   }
  }] -row 0 -column 1 -sticky e -pady 2 -padx 2
  grid columnconfigure $v.fp.nma.header 0 -weight 1

  #.qmtool.analysis.fp select $tabid
  
}

################################################
# Plot Molecular orbitals                      #
################################################
proc ::QMtool::call_qm_orbitals { f } {
   
   variable plothandle2
   variable qmorbtable
   variable qmspan
   variable hasMOs
   global vmd_frame
   
   if {$hasMOs == "No"} {
   ## Create an error message to return saying no orbital energies
   tk_messageBox -message "Error loading orbital information. Check your molecule." -title "Orbitals not found" -icon error \
    -type ok -parent .qmtool.fp.analysis.fp.orbitals   
     return
   }
   
   #if {[winfo exist .qmtool.fp.analysis.fp.orbitals] == 0 } { return }
   

   # check the table exist already, if so redo it
   if {[winfo exist .qmtool.fp.analysis.fp.orbitals.header.tableresults]} {
      ::QMtool::updateOrbitalsTable
   } else {
   ## Creates table with orbitals, labels and energies
   grid [ttk::frame $f.tableresults] -row 1 -column 0 -sticky nswe -padx 2 -pady 4 -columnspan 8
   grid columnconfigure $f.tableresults 0 -weight 1
   grid rowconfigure $f.tableresults 0 -weight 1

   grid [ttk::frame $f.tableresults.ptrfrm] -row 0 -column 0 -sticky nswe -pady 2 
   grid columnconfigure $f.tableresults.ptrfrm 0 -weight 1
   grid rowconfigure $f.tableresults.ptrfrm 0 -weight 1
   
   set str "Orbitals Table"
   #grid [ttk::label $f.tableresults.ptrfrm.prt -text "$str" -image $QWIKMD::arrowDown -compound left] -row 0 -column 0 -stick w -pady 2   
  
   grid [ttk::frame $f.tableresults.ptrfrm.plot] -row 0 -column 1 -sticky ns -padx 0 
   grid columnconfigure $f.tableresults.ptrfrm.plot 0 -weight 1
   grid rowconfigure $f.tableresults.ptrfrm.plot 0 -weight 1

   grid [ttk::button $f.tableresults.ptrfrm.plot.btt -text "Plot" -padding "2 2 2 2" -width 10 -command {::MultiPlot::plotOrbitals QMtool}] -row 0 -column 0 -sticky e -pady 2 -padx 2
   #puts "Inside call_qm_orbitals: qmorbtable = $qmorbtable"
   grid [ttk::button $f.tableresults.ptrfrm.plot.bttclear -text "Clear Sel." -padding "2 2 2 2" -width 10 -command {
           if {[llength $::QMtool::qmorbtable] == 0} {
               return
           }
           $::QMtool::qmorbtable selection clear 0 end
           #set ::QMtool::qmorbrep ""
           ::QMtool::showOrbitals
       }] -row 0 -column 1 -sticky e -pady 2 -padx 2
  
   grid [ttk::frame $f.tableresults.ptrfrm.span] -row 0 -column 0 -sticky ns -padx 2 
   grid columnconfigure $f.tableresults.ptrfrm.span 0 -weight 1
   grid rowconfigure $f.tableresults.ptrfrm.span 0 -weight 1

   grid [ttk::label $f.tableresults.ptrfrm.span.lbl -text "Orbitals around HOMO:"] -row 0 -column 0 -stick ns -padx 2
   grid [ttk::entry $f.tableresults.ptrfrm.span.val -width 2 -textvariable ::QMtool::qmspan] -row 0 -column 1 -stick nsw -padx 4
  
   set qmspan 5
   #set qmspan $f.tableresults.ptrfrm.span.val
   bind .qmtool.fp.analysis.fp.orbitals.header.tableresults.ptrfrm.span.val <Return> {
       if {[$::QMtool::qmorbtable size] != 0} {
           set val [string trim $::QMtool::qmspan]
           if {[string is integer $val] != 1 || $val < 0} {
               tk_messageBox -message "Please select a valid number of orbitals >=0" -title "Number of Orbitals" -icon error -parent .qmtool.fp.analysis.fp.orbitals.header
            set ::QMtool::qmspan 5
            return
           }
           ::QMtool::updateOrbitalsTable
           ::MultiPlot::plotOrbitals QMtool
       }
       return
   }  

   grid [ttk::frame $f.tableresults.fcolapse] -row 2 -column 0 -sticky nswe -padx 4

   grid columnconfigure $f.tableresults.fcolapse 0 -weight 1
   grid rowconfigure $f.tableresults.fcolapse 0 -weight 1

   grid [ttk::frame $f.tableresults.fcolapse.t1] -row 0 -column 0 -sticky nswe -padx 4

   grid columnconfigure $f.tableresults.fcolapse.t1 0 -weight 1
   grid rowconfigure $f.tableresults.fcolapse.t1 0 -weight 1

   option add *Tablelist.activeStyle       frame
   
   set fro2 $f.tableresults.fcolapse.t1

   option add *Tablelist.movableColumns    no
   option add *Tablelist.labelCommand      tablelist::sortByColumn

       tablelist::tablelist $fro2.tb -columns {\
           0 "Number" center
           0 "Description" center
           0 "Energy (a.u.)" center
       }\
       -yscrollcommand [list $fro2.scr1 set] -xscrollcommand [list $fro2.scr2 set] \
               -showseparators 0 -labelrelief groove  -labelbd 1 -selectforeground black\
               -foreground black -background white -state normal -selectmode extended -height 10 -stretch all -stripebackgroun white -exportselection true\
               

   $fro2.tb columnconfigure 0 -selectbackground cyan -name OrbNum -maxwidth 5 -sortmode integer
   $fro2.tb columnconfigure 1 -selectbackground cyan -name Descr -maxwidth 10 -sortmode dictionary
   $fro2.tb columnconfigure 2 -selectbackground cyan -name Energy -maxwidth 10 -sortmode real

   grid $fro2.tb -row 0 -column 0 -sticky news 
   
   ##Scroll_BAr V
   scrollbar $fro2.scr1 -orient vertical -command [list $fro2.tb  yview]
    grid $fro2.scr1 -row 0 -column 1  -sticky ens

   ## Scroll_Bar H
   scrollbar $fro2.scr2 -orient horizontal -command [list $fro2.tb xview]
   grid $fro2.scr2 -row 1 -column 0 -sticky swe

   }

   set qmorbtable .qmtool.fp.analysis.fp.orbitals.header.tableresults.fcolapse.t1.tb
   
   bind .qmtool.fp.analysis.fp.orbitals.header.tableresults.fcolapse.t1.tb <<TablelistSelect>>  {
       ::QMtool::showOrbitals
   }  
  
   # If already initialized, just turn on
   if { [winfo exists .qmtool.plotorb] } {
      wm deiconify .qmtool.plotorb
      raise .qmtool.plotorb
      return
   }

   set mo [toplevel ".qmtool.plotorb"]
   wm title $mo "QMtool 2.0 - Molecular Orbitals Plot"
   wm resizable $mo 1 1
   set width 520 ;# in pixels
   set height 440 ;# in pixels
   wm geometry $mo ${width}x${height}

   grid [frame $mo.plot] -row 0 -column 0 -sticky news
  
   set qmID $::QMtool::molid
   set molname $::QMtool::molselected
   set title "Molecule: $molname ID: $qmID Orbitals Energy"
   set ylab "Orbital\nEnergy\n(a.u.)"
   set xlab ""
   set plothandle2 [multiplot embed $mo.plot -xsize 500 -ysize 400 -lines -title $title \
     -xlabel $xlab -ylabel $ylab -marker circle]
   set plotwindow [$plothandle2 getpath]

   ## Add more menus to clear and close plots not included by
   ## default in the multiplot windows.
   menubutton $plotwindow.menubar.clear -text "Clear" \
     -underline 0 -menu $plotwindow.menubar.clear.menu

   $plotwindow.menubar.clear config -width 5

   menu $plotwindow.menubar.clear.menu -tearoff 0

   $plotwindow.menubar.clear.menu add command -label "Clear Plot"
   $plotwindow.menubar.clear.menu entryconfigure 0 -command {
       #destroy .qmtool.analysis.fp.orbitals.plot
       #$::QMtool::qmorbtable delete 0 end
       #set ::QMtool::plothandle {}
   }

   menubutton $plotwindow.menubar.close -text "Close" \
     -underline 0 -menu $plotwindow.menubar.close.menu
   $plotwindow.menubar.close config -width 5
   menu $plotwindow.menubar.close.menu -tearoff 0
   $plotwindow.menubar.close.menu add command -label "Close Plot"
   $plotwindow.menubar.close.menu entryconfigure 0 -command {
      destroy .qmtool.plotorb.plot
      $::QMtool::qmorbtable delete 0 end
      set ::QMtool::plothandle2 {}
   }

   #pack $plotwindow.menubar.clear -side left    
   grid $plotwindow.menubar.clear -row 0 -column 2 -sticky nwes
   #pack $plotwindow.menubar.close -side left
   grid $plotwindow.menubar.close -row 0 -column 3 -sticky nwes

   grid $plotwindow -row 0 -column 0 -sticky nwes
   
   trace add variable vmd_frame($::QMtool::molid) write ::QMtool::updateOrbitalsTable
   mol top $::QMtool::molid
   ::QMtool::updateOrbitalsTable
  
}

proc ::QMtool::updateOrbitalsTable {args} {
    #puts "In QMtool::updateOrbitalsTable"
    variable qmorbtable
    variable molid
    variable qmspan
    ### save table selection to be kept after updating the values
    set tableselection [$qmorbtable curselection]

    $qmorbtable delete 0 end
    set descr ""
    set norbs [molinfo $molid get numorbitals]
    set homo [lindex [molinfo $molid get homo] 0 0]

    set lowerlimit 0
    set lowerlimit [expr $homo - $qmspan]
    if {$lowerlimit < 0} {
        set lowerlimit 0
    }

    set uperlimit [expr $homo + $qmspan]
    if {$uperlimit > $norbs} {
        set uperlimit [expr $norbs - 1]
    }

    set energies [lrange [join [lindex [molinfo $molid get orbenergies] 0]] $lowerlimit $uperlimit]
    set j 0
    for {set i $lowerlimit} {$i <= $uperlimit} {incr i} {
        set diff [expr abs($homo - $i)]
        if {$i <= [expr $homo -1]} {
            set descr "HOMO-$diff"
        } elseif {$i == $homo} {
            set descr "HOMO"
        } elseif {$i == [expr $homo + 1]} {
            set descr "LUMO"
        } elseif {$i > [expr $homo + 1]} {
            incr diff -1
            set descr "LUMO+$diff"
        }
        $qmorbtable insert end [list $i $descr [lindex $energies $j]]
        incr j
    }
    $qmorbtable selection set $tableselection
    return
}

################################################################################
### Proc to display the orbitals representation
################################################################################
proc ::QMtool::showOrbitals {} {
    variable qmorbtable
    variable molid
    variable qmorbrep
    variable orbprevrep
    variable orbprevreplist
    
    #puts "...................."
    #puts "In ::QMtool::showOrbitals, intial orbprevreplist = $orbprevreplist"
    
    set tbindex [$qmorbtable curselection]
    if {$tbindex == -1 || [llength $tbindex] == 0} {
        for {set i 0} {$i < [llength $orbprevreplist]} {incr i} {
            set index [lsearch -index 0 $qmorbrep [lindex $orbprevreplist $i] ]
            if {[lindex $qmorbrep $index 0] == -1} {
                continue
            } 
            foreach repname [lindex $qmorbrep $index 1] {
                #puts "$repname being deleted (1) for molid $molid"
                mol delrep [QWIKMD::getrepnum $repname $molid] $molid
            }
        }
        set orbprevrep [list]
        return
    }
    for {set i 0} {$i < [llength $orbprevreplist]} {incr i} {
        if {[lsearch $tbindex [lindex $orbprevreplist $i]] == -1} {
            set index [lsearch -index 0 $qmorbrep [lindex $orbprevreplist $i] ] 

            foreach repname [lindex $qmorbrep $index 1] {
                #puts "$repname being deleted (2) for molid $molid"
                mol delrep [QWIKMD::getrepnum $repname $molid] $molid
            }
            set qmorbrep [lreplace $qmorbrep $index $index]
        }
    }

    set frm [molinfo $molid get frame]
    set orblist [list]
    foreach ind $tbindex {
        lappend orblist [$qmorbtable cellcget $ind,0 -text]
        set orb [lindex $orblist end]
        if {[llength $qmorbrep] > 0 && [lsearch -index 0 $qmorbrep $orb] != -1} {
            continue
        }
        set listrep [list]

        mol addrep $molid
        set rep [expr [molinfo $molid get numreps] -1]
        lappend listrep [mol repname $molid $rep]

        mol modcolor $rep $molid ColorID 0
        mol modstyle $rep $molid "Orbital 0.050000 $orb 0 0 0.125 1 0 0 0 1"
        mol modselect $rep $molid all
        mol modmaterial $rep $molid Glossy
        
        mol addrep $molid
        set rep [expr [molinfo $molid get numreps] -1]
        lappend listrep [mol repname $molid $rep]

        mol modcolor $rep $molid ColorID 3
        mol modstyle $rep $molid "Orbital -0.050000 $orb 0 0 0.125 1 0 0 0 1"
        mol modselect $rep $molid all
        mol modmaterial $rep $molid Glossy

        lappend qmorbrep [concat $orb [list $listrep] ]
    }
    
    #puts "In ::QMtool::showOrbitals, final orblist = $orblist" 
    set orbprevreplist $orblist
    #puts "In ::QMtool::showOrbitals, final orbprevreplist = $orbprevreplist" 
}

################################################
# Plot SCF energies.                           #
################################################

proc ::QMtool::plot_scf_energies { notused } {
   global vmd_frame
   variable ::QMtool::plothandle
   variable molid
   variable molselected
   #puts "In plot_scf_energies, plothandle: $plothandle"
   #if {[winfo exists $f.plot] == 0} { return }
   
   variable scfenergies
   if {[llength $scfenergies]<=1} {
      tk_messageBox -message "Error loading SCF energies. Check your molecule." -title "Energies not found" -icon error \
    -type ok -parent .qmtool.fp.analysis.fp.energies
      #if { [winfo exists .qmtool.menu.analysis] } { .qmtool.menu.analysis entryconfigure 0 -state disabled }
      return 
   }
   
   # If already initialized, just turn on
   if { [winfo exists .qmtool.plot] } {
      wm deiconify .qmtool.plot
      raise .qmtool.plot
      return
   }

   set f [toplevel ".qmtool.plotmain"]
   wm title $f "QMtool 2.0 - Molecule $molselected ID $molid Energies Plot"
   wm resizable $f 1 1
   set width 630 ;# in pixels
   set height 450 ;# in pixels
   wm geometry $f ${width}x${height}

   grid [frame $f.plot] -row 1 -column 0 -sticky news
   
   set x {}
   set i 0
   foreach e $scfenergies {
      lappend x $i
      lappend y [lindex $e 1]
      incr i
   }
   set title "SCF energies (relative to the first value)"
   set plothandle [multiplot embed $f.plot -xsize 600 -ysize 400 -x $x -y $y -lines -title $title \
      -xlabel "Geometry Optimization Step" -ylabel "E(kcal/mol)" -marker circle -callback ::QMtool::marker_clicked]
   #$plothandle replot
   set plotwindow [$plothandle getpath]

   ## Add more menus to clear and close plots not included by
   ## default in the multiplot windows.
   #menubutton $plotwindow.menubar.clear -text "Clear" \
   -underline 0 -menu $plotwindow.menubar.clear.menu

   #$plotwindow.menubar.clear config -width 5
   #menu $plotwindow.menubar.clear.menu -tearoff 0
   #$plotwindow.menubar.clear.menu add command -label "Clear Plot"
   
   
   menubutton $plotwindow.menubar.close -text "Close" \
   -underline 0 -menu $plotwindow.menubar.close.menu
   $plotwindow.menubar.close config -width 5
   menu $plotwindow.menubar.close.menu -tearoff 0
   $plotwindow.menubar.close.menu add command -label "Close Plot"
   $plotwindow.menubar.close.menu entryconfigure 0 -command {
     destroy .qmtool.plotmain.plot
     set ::QMtool::plothandle {}
   }

   #pack $plotwindow.menubar.clear -side left
   #grid $plotwindow.menubar.clear -row 0 -column 2 -sticky nwes
   #pack $plotwindow.menubar.close -side left
   grid $plotwindow.menubar.close -row 0 -column 3 -sticky nwes

   grid $plotwindow -row 0 -column 0 -sticky nwes
   
   # Update frame to display frame marker in new plot
   trace add variable vmd_frame($::QMtool::molid) write ::QMtool::update_frame
   $plothandle replot
   #update_frame internal [molinfo top] w

}

proc ::QMtool::update_frame { name molid op } {
  # name == vmd_frame
  # molid == molecule id of the newly changed frame
  # op == w

  if { $molid != [molinfo top] } {
    return
  }
  set f [molinfo $molid get frame]
  # refresh the frame marker in the plot
  display_marker $f
}


# Display frame marker in plot at given frame
proc ::QMtool::display_marker { f } {
  variable ::QMtool::plothandle
  #puts "In display_marker, plothandle: $plothandle"
  if [info exists ::QMtool::plothandle] {
    # detect if plot was closed
    if [catch {$plothandle getpath}] {
      unset plothandle
    } else {
      # we tinker a little with Multiplot's internals to get access to its Tk canvas
      # necessary because Multiplot does not expose an interface to draw & delete
      # objects without redrawing the whole plot - which takes too long for this
      set ns [namespace qualifiers $::QMtool::plothandle]
      puts $ns
      set xmin [set ${ns}::xmin]
      set xmax [set ${ns}::xmax]
      # Move plot boundaries if necessary
      if { $f < $xmin } {
        set xmax [expr { $xmax + $f - $xmin }]
        set xmin $f
        $plothandle configure -xmin $xmin -xmax $xmax -plot
      }
      if { $f > $xmax } {
        set xmin [expr { $xmin + $f - $xmax }]
        set xmax $f
        $plothandle configure -xmin $xmin -xmax $xmax -plot
      }

    set y1 [set ${ns}::yplotmin]
    set y2 [set ${ns}::yplotmax]
    set xplotmin [set ${ns}::xplotmin]
    set scalex [set ${ns}::scalex]
    #puts "[set ${ns}::scalex]"
    #puts "In display_marker: scalex $scalex"
    set x [expr $xplotmin+($scalex*($f-$xmin))]

    set canv "[set ${ns}::w].f.cf"
    $canv delete frame_marker
    $canv create line  $x $y1 $x $y2 -fill blue -tags frame_marker
    }
  }
}

# Callback for click on marker in 2 cv plot
proc ::QMtool::marker_clicked { index x y color marker } {

  animate goto [expr {$index - 1 }]
}

########################################################
# Plots the spectrum by replacing each spectral line   #
# with a lorentzian function with a half width at half #
# maximum of $hwhm.                                    #
########################################################

proc ::QMtool::plot_spectrum { {hwhm 2} {npoints 600}} {
   variable lineintensities
   variable linewavenumbers
   variable nimag
   set c 299792458; # lightspeed in m/s
   set pi 3.141592654
   set giga 1000000000.0
   set cutoff 5*$hwhm

   set minwvn [lindex $linewavenumbers 0]
   set maxwvn [lindex $linewavenumbers end]
   set rangewvn [expr $maxwvn-$minwvn]
   #set maxlam [expr 1.0e4/$minwvn]; # in microm
   #set minlam [expr 1.0e4/$maxwvn]; # in micromm
   #set rangelam [expr $maxlam-$minlam]
   puts "rangewvn: $maxwvn - $minwvn = $rangewvn"
   #puts "rangefreq [expr 100.0*$c*$maxwvn/$giga] - [expr 100.0*$c*$minwvn/$giga] = [expr 100.0*$c*$rangewvn/$giga]"
   puts "Computing spectrum..."

   set deltanu [expr $rangewvn/$npoints]
   
   set binnedwvn {}
   foreach wvn $linewavenumbers {
      set pos [expr ($wvn-$minwvn)]
      set n [expr int($pos/$deltanu)]
      set delta [expr $pos-$deltanu*$n]
      if {[expr $delta/$deltanu] >= 0.5} {
 	 lappend binnedwvn [expr $minwvn+($n+1)*$deltanu]; # NewFreq(int(position)+1)
      } else {
 	 lappend binnedwvn [expr $minwvn+($n)*$deltanu];   # NewFreq(int(position))
      }
   }

   set specint {}
   set specfreq {}
   for {set i 1} {$i<=$npoints} {incr i} {
      set intensity 0.0
      set wvn [expr $minwvn+$rangewvn*$i/$npoints]; # in 1/cm

      foreach lineint $lineintensities linefreq $binnedwvn {
	 set offset [expr $wvn-$linefreq]
	 if {$offset<$cutoff} {
	    set intensity [expr $intensity + $lineint * [lorentz [expr $offset/$hwhm]]]
	 }
      }
      lappend specint $intensity
      lappend specwvn $wvn
      #puts "$wvn $freq  $linefreq $intensity"
   }

   set plothandle [multiplot -x $specwvn -y $specint -lines -title "Harmonic spectrum" \
      -xlabel "wavenumber in 1/cm" -ylabel "intensity"]
   $plothandle add $linewavenumbers $lineintensities -nolines -marker circle -plot
}

proc ::QMtool::lorentz { offset } {
   return [expr 1.0 / (1.0 + $offset * $offset)]
}

proc ::QMtool::thermochemistry {} {
   variable selectcolor

   # If already initialized, just turn on
   if { [winfo exists .qmtool_thermo] } {
      wm deiconify .qmtool_thermo
      focus .qmtool_thermo
      return
   }

   set v [toplevel ".qmtool_thermo"]
   wm title $v "Thermochemistry"
   wm resizable $v 0 0

   frame $v.energy

   ############## frame for selected component energies #################
   labelframe $v.energy.comp -bd 2 -relief ridge -text "Energy of selected component" -padx 1m -pady 1m

   label $v.energy.comp.templabel -text "Temperature: "
   entry $v.energy.comp.tempentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::temperature
   grid $v.energy.comp.templabel -row 1 -column 0 -sticky w
   grid $v.energy.comp.tempentry -row 1 -column 1

   label $v.energy.comp.evaclabel -text "E(gas): "
   entry $v.energy.comp.evacentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::Evacuum
   grid $v.energy.comp.evaclabel -row 2 -column 0 -sticky w
   grid $v.energy.comp.evacentry -row 2 -column 1

   label $v.energy.comp.esolvlabel -text "E(solv): "
   entry $v.energy.comp.esolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::Esolv
   grid $v.energy.comp.esolvlabel -row 3 -column 0 -sticky w
   grid $v.energy.comp.esolventry -row 3 -column 1

   label $v.energy.comp.gvaclabel -text "G(gas): "
   entry $v.energy.comp.gvacentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::Gvacuum
   grid $v.energy.comp.gvaclabel -row 4 -column 0 -sticky w
   grid $v.energy.comp.gvacentry -row 4 -column 1

   label $v.energy.comp.gsolvlabel -text "G(solv): "
   entry $v.energy.comp.gsolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::Gsolv
   grid $v.energy.comp.gsolvlabel -row 5 -column 0 -sticky w
   grid $v.energy.comp.gsolventry -row 5 -column 1

   label $v.energy.comp.dgsolvlabel -text "Delta G_solvation: "
   entry $v.energy.comp.dgsolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGsolvation
   grid $v.energy.comp.dgsolvlabel -row 6 -column 0 -sticky w
   grid $v.energy.comp.dgsolventry -row 6 -column 1

   ############## frame for selected component energies #################
   labelframe $v.energy.react -bd 2 -relief ridge -text "Reaction Energies" -padx 1m -pady 1m

   label $v.energy.react.templabel -text "Delta G_solvation(Educts): "
   entry $v.energy.react.tempentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGsolvationE
   grid $v.energy.react.templabel -row 1 -column 0 -sticky w
   grid $v.energy.react.tempentry -row 1 -column 1

   label $v.energy.react.evaclabel -text "Delta G_solvation(Products): "
   entry $v.energy.react.evacentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGsolvationP
   grid $v.energy.react.evaclabel -row 2 -column 0 -sticky w
   grid $v.energy.react.evacentry -row 2 -column 1

   label $v.energy.react.dgsolvlabel -text "Total Delta G_solvation: "
   entry $v.energy.react.dgsolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGsolvationTotal
   grid $v.energy.react.dgsolvlabel -row 3 -column 0 -sticky w
   grid $v.energy.react.dgsolventry -row 3 -column 1

   label $v.energy.react.esolvlabel -text "Delta G_reaction(gas): "
   entry $v.energy.react.esolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGreactionGas
   grid $v.energy.react.esolvlabel -row 4 -column 0 -sticky w
   grid $v.energy.react.esolventry -row 4 -column 1

   label $v.energy.react.gvaclabel -text "Delta G_reaction(solution): "
   entry $v.energy.react.gvacentry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGreactionSol
   grid $v.energy.react.gvaclabel -row 5 -column 0 -sticky w
   grid $v.energy.react.gvacentry -row 5 -column 1

   label $v.energy.react.gsolvlabel -text "Delta G_reaction(cycle): "
   entry $v.energy.react.gsolventry -relief sunken -width 14 -justify right -state readonly \
      -font {tkFixed 9} -textvariable ::QMtool::dGreactionCycle
   grid $v.energy.react.gsolvlabel -row 6 -column 0 -sticky w
   grid $v.energy.react.gsolventry -row 6 -column 1

   pack $v.energy.comp $v.energy.react -side left -padx 1m


   ############## frame for educt file list #################
   labelframe $v.educt -bd 2 -relief ridge -text "Educts" -padx 1m -pady 1m
   frame $v.educt.list
   scrollbar $v.educt.list.scroll -command "$v.educt.list.list yview"
   listbox $v.educt.list.list -activestyle dotbox -yscroll "$v.educt.list.scroll set" \
      -width 60 -height 3 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::QMtool::eductlist
   frame  $v.educt.list.buttons
   button $v.educt.list.buttons.add -text "Add"    -command { ::QMtool::add_educt }

   button $v.educt.list.buttons.delete -text "Delete" -command {
      foreach i [.qmtool_thermo.educt.list.list curselection] {
	 ::QMtool::molecule_delete [lindex $::QMtool::eductlist $i 0]
	 .qmtool_thermo.educt.list.list delete $i
	 set ::QMtool::thermElist [lreplace $::QMtool::thermElist $i $i] 
      }
      ::QMtool::update_reaction 
   }
   pack $v.educt.list.buttons.add $v.educt.list.buttons.delete -expand 1 -fill x
   pack $v.educt.list.list -side left  -fill x -expand 1
   pack $v.educt.list.scroll $v.educt.list.buttons -side left -fill y -expand 1
   pack $v.educt.list -expand 1 -fill x

   ############## frame for product file list #################
   labelframe $v.product -bd 2 -relief ridge -text "Products" -padx 1m -pady 1m
   frame $v.product.list
   scrollbar $v.product.list.scroll -command "$v.product.list.list yview"
   listbox $v.product.list.list -activestyle dotbox -yscroll "$v.product.list.scroll set" \
      -width 60 -height 3 -setgrid 1 -selectmode browse -selectbackground $selectcolor \
      -listvariable ::QMtool::productlist
   frame  $v.product.list.buttons
   button $v.product.list.buttons.add -text "Add"    -command { ::QMtool::add_product  }
   button $v.product.list.buttons.delete -text "Delete" -command {
      foreach i [.qmtool_thermo.product.list.list curselection] {
	 ::QMtool::molecule_delete [lindex $::QMtool::productlist $i 0]
	 .qmtool_thermo.product.list.list delete $i
	 set ::QMtool::thermPlist [lreplace $::QMtool::thermPlist $i $i] 
      }
      ::QMtool::update_reaction 
   }
   pack $v.product.list.buttons.add $v.product.list.buttons.delete -expand 1 -fill x
   pack $v.product.list.list -side left  -fill x -expand 1
   pack $v.product.list.scroll $v.product.list.buttons -side left -fill y -expand 1
   pack $v.product.list -expand 1 -fill x

   pack $v.energy -pady 1m -fill x
   pack $v.educt $v.product -padx 1m -pady 1m -fill x

   bind $v.educt.list.list <<ListboxSelect>> {
      ::QMtool::update_educt
   }

   bind $v.product.list.list <<ListboxSelect>> {
      ::QMtool::update_product
   }
}

proc ::QMtool::add_educt {} {
   ::QMtool::opendialog log
   array set molecules [join $::QMtool::molnamelist]
   variable Evacuum
   variable Esolv
   variable Gvacuum
   variable Gsolv
   variable EGvacuum
   variable EGsolv
   variable dGsolvation
   variable thermElist
   variable eductlist
   lappend thermElist [list $Evacuum $Esolv $Gvacuum $Gsolv $EGvacuum $EGsolv $dGsolvation]
   lappend eductlist [array get molecules $::QMtool::molid]
   update_reaction
}

proc ::QMtool::update_educt {} {
   set i [.qmtool_thermo.educt.list.list curselection]
   variable thermElist
   set E [lindex $thermElist $i]
   variable Evacuum     [lindex $E 0]
   variable Esolv       [lindex $E 1]
   variable Gvacuum     [lindex $E 2]
   variable Gsolv       [lindex $E 3]
   variable EGvacuum    [lindex $E 4]
   variable EGsolv      [lindex $E 5]
   variable dGsolvation [lindex $E 6]
}

proc ::QMtool::add_product {} {
   ::QMtool::opendialog log
   array set molecules [join $::QMtool::molnamelist]
   variable Evacuum
   variable Esolv
   variable Gvacuum
   variable Gsolv
   variable EGvacuum
   variable EGsolv
   variable dGsolvation
   variable thermPlist
   variable productlist
   lappend thermPlist [list $Evacuum $Esolv $Gvacuum $Gsolv $EGvacuum $EGsolv $dGsolvation]
   lappend productlist [array get molecules $::QMtool::molid]
   update_reaction
}

proc ::QMtool::update_product {} {
   set i [.qmtool_thermo.product.list.list curselection]
   variable thermPlist
   set E [lindex $thermPlist $i]
   variable Evacuum     [lindex $E 0]
   variable Esolv       [lindex $E 1]
   variable Gvacuum     [lindex $E 2]
   variable Gsolv       [lindex $E 3]
   variable EGvacuum    [lindex $E 4]
   variable EGsolv      [lindex $E 5]
   variable dGsolvation [lindex $E 6]
}

proc ::QMtool::update_reaction {} {
   set EtotEGvacuum 0.0
   set EtotEGsolv   0.0
   variable dGsolvationE 0.0

   variable thermElist
   foreach E $thermElist {
      #if {[llength $E]} {
	 set EtotEGvacuum    [expr [lindex $E 4]+$EtotEGvacuum]
	 set EtotEGsolv      [expr [lindex $E 5]+$EtotEGsolv  ]
	 set dGsolvationE    [expr [lindex $E 6]+$dGsolvationE]
      #}
   }

   set PtotEGvacuum 0.0
   set PtotEGsolv   0.0
   variable dGsolvationP 0.0

   variable thermPlist
   foreach E $thermPlist {
      #if {[llength $E]} {
	 set PtotEGvacuum    [expr [lindex $E 4]+$PtotEGvacuum]
	 set PtotEGsolv      [expr [lindex $E 5]+$PtotEGsolv  ]
	 set dGsolvationP    [expr [lindex $E 6]+$dGsolvationP]
      #}
   }
   variable dGsolvationTotal [expr -$dGsolvationE+$dGsolvationP]
   variable dGreactionGas    [expr $PtotEGvacuum -$EtotEGvacuum]
   variable dGreactionSol    [expr $PtotEGsolv   -$EtotEGsolv]
   variable dGreactionCycle  [expr $dGsolvationTotal+$dGreactionGas]
   if {[llength $dGsolvationE]}     { set dGsolvationE [format "%14.2f" $dGsolvationE]}
   if {[llength $dGsolvationP]}     { set dGsolvationP [format "%14.2f" $dGsolvationP]}
   if {[llength $dGsolvationTotal]} { set dGsolvationTotal [format "%14.2f" $dGsolvationTotal]}
   if {[llength $dGreactionGas]}    { set dGreactionGas    [format "%14.2f" $dGreactionGas]}
   if {[llength $dGreactionSol]}    { set dGreactionSol    [format "%14.2f" $dGreactionSol]}
   if {[llength $dGreactionCycle]}  { set dGreactionCycle  [format "%14.2f" $dGreactionCycle]}
}

proc ::QMtool::normalmode_gui { v } {
   variable normalmodes
   variable qmnmatable
   variable lineintensities
   variable linewavenumbers
   variable formatnormalmodes {}

   if {[llength $lineintensities] == "0"} {
   tk_messageBox -message "Error loading frequencies information. Check your molecule." -title "Frequencies not found" -icon error \
    -type ok -parent .qmtool.fp.analysis.fp.nma.header  
     return
   }
 


   #array set modelist $normalmodes
   set i 0
   foreach intens $lineintensities wvn $linewavenumbers {
      lappend formatnormalmodes [format "%3i: %8.2f %8.2f" $i $wvn $intens]
      incr i
   }

   variable selectcolor

   ############## frame for molecule list #################
   grid [frame $v.mode] -column 0 -row 1 -sticky w -columnspan 2
      
   set fro2 $v.mode

   # NEW table for normalmode, from qwikmd
   option add *Tablelist.activeStyle       frame
   option add *Tablelist.background        gray98
   option add *Tablelist.stripeBackground  #e0e8f0
   option add *Tablelist.setGrid           no
   option add *Tablelist.movableColumns    no
  
   tablelist::tablelist $fro2.tb -width 40 -columns { 
                0 "Mode"  center 
                0 "Frequency (cm-1)"  center
                0 "IR Intensities" center
                } -yscrollcommand [list $fro2.scr1 set] -showseparators 0 -labelrelief groove  -labelbd 1 -selectbackground cyan \
                -selectforeground black -foreground black -background white -state normal -selectmode browse -stretch "0 1 2" -stripebackgroun white -height 10
   
   #-editstartcommand QWIKMD::cellStartEditPtcl -editendcommand QWIKMD::cellEndEditPtcl -forceeditendcommand true -editselectedonly true

   grid $fro2.tb -row 0 -column 0 -sticky news

   #$fro2.tb columnconfigure 0 -width 0
    ##Scroll_BAr V
   grid [scrollbar $fro2.scr1 -orient vertical -command [list $fro2.tb  yview]] -row 0 -column 1  -sticky ens

   set qmnmatable $fro2.tb

   bind $fro2.tb <<TablelistSelect>>  {
      set mode [.qmtool.fp.analysis.fp.nma.header.mode.tb curselection]
      ::QMtool::show_normalmode $mode
   }  

   # need to populate the table
   for {set i 0} {$i < [llength $formatnormalmodes]} {incr i} {
      $qmnmatable insert end [lindex $formatnormalmodes $i]
   }

   grid [frame $v.buttons] -column 0 -row 2 -sticky w
   grid [label $v.buttons.label -text "Scaling factor"] -column 0 -row 0 -sticky w
   grid [spinbox $v.buttons.spinb -from 0 -to 10 -increment 0.05 -width 5 \
      -textvariable ::QMtool::normalmodescaling -command {
	      set mode [$v.mode.list.list curselection]
	      if {[llength $mode]} {::QMtool::show_normalmode $mode 0}
       }] -column 1 -row 0 -sticky w
   grid [label $v.buttons.label2 -text "Number of animation steps"] -column 0 -row 1 -sticky w
   grid [entry $v.buttons.steps -textvariable ::QMtool::normalmodesteps -width 4] -column 1 -row 1 -sticky w
   grid [checkbutton $v.buttons.arrows -text "Show arrows" -variable ::QMtool::normalmodearrows] -column 0 -row 2 -sticky w

}

proc ::QMtool::show_normalmode { mode {ncycles 1} {movie "nomovie"} {arrows ""}} {
   variable normalmodes
   variable normalmodearrows
   if {![llength $arrows]} { 
      set arrows $normalmodearrows
   } elseif {$arrows!="arrows"} { set arrows 0 }

   draw delete all
   draw color yellow
   variable molid
   variable normalmodescaling
   set sel [atomselect $molid all]
   foreach atommode [lindex $normalmodes $mode] coord [$sel get {x y z}] {
      if {$arrows && [veclength $atommode]>0.2} {
	 ::QMtool::arrow $molid $coord [vecadd $coord [vecscale $atommode $normalmodescaling]] 0.1
      }
   }

   set initialcoords [$sel get {x y z}]
   set deg2rad [expr 3.14159265358979/180.0]
   variable normalmodesteps
   for {set cyc 0} {$cyc<$ncycles} {incr cyc} {
      for {set i 1} {$i<=$normalmodesteps} {incr i} {
	 set angle [expr 360.0*$i/double($normalmodesteps)]
	 set diff [expr $normalmodescaling*sin($deg2rad*$angle)]
	 puts "angle=$angle: $diff"
	 set poslist {}
	 foreach atommode [lindex $normalmodes $mode] coord $initialcoords {
	    lappend poslist [vecadd $coord [vecscale $atommode $diff]]
	 }
	 if {$movie=="movie"} { animate dup $molid }
	 $sel lmoveto $poslist
	 display update
      }
   }
   $sel lmoveto $initialcoords
}


proc ::QMtool::arrow {mol start end {rad 1} {res 6}} {
    # an arrow is made of a cylinder and a cone
    set middle [vecadd $start [vecscale 0.85 [vecsub $end $start]]]
    graphics $mol cone $middle $end radius [expr $rad*2.0] resolution $res
    graphics $mol cylinder $start $middle radius $rad resolution $res
    #puts "$middle $end"
}

