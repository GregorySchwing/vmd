
proc ::QMtool::qmtool_gui {} {
   variable w 
   variable selectcolor

   # If already initialized, just turn on
   if { [winfo exists .qmtool] } {
      wm deiconify .qmtool
      raise .qmtool
      return
   }

   set w [toplevel ".qmtool"]
   wm title $w "QMtool 2.0 - Setup and Analysis for QM Simulations"
   wm resizable $w 1 1
   set width 440 ;# in pixels
   set height 330 ;# in pixels
   wm geometry $w ${width}x${height}+68+125
   grid columnconfigure $w 0 -weight 1
   grid columnconfigure $w 1 -weight 0
   grid rowconfigure $w 0 -weight 0
   grid rowconfigure $w 1 -weight 1

   # trying menubar again
   grid [frame $w.menubar -relief raised -bd 2] -row 0 -column 0 -sticky nswe -pady 2 -padx 2
   grid columnconfigure $w.menubar 4 -weight 1
   grid rowconfigure $w.menubar 0 -weight 1

   grid [menubutton $w.menubar.file -text "File" -width 8 -menu $w.menubar.file.menu] -row 0 -column 0 -sticky ew
   grid [menubutton $w.menubar.edit -text "Edit" -width 8 -menu $w.menubar.edit.menu] -row 0 -column 3 -sticky ew
   grid [menubutton $w.menubar.help -text "Help" -width 8 -menu $w.menubar.help.menu] -row 0 -column 4 -sticky ew
      
   # File
   menu $w.menubar.file.menu -tearoff no
   $w.menubar.file.menu add command -label "Load Gaussian File" -command {::QMtool::opendialog log}
   $w.menubar.file.menu add command -label "Write Input File" -command {::QMtool::edit_input_ok ask}

   # Edit
   menu $w.menubar.edit.menu -tearoff no
   $w.menubar.edit.menu add command -label "Edit Structure" -command {::QMtool::molefactureStart}
   $w.menubar.edit.menu add command -label "Reset Data" -command {
      ::QMtool::clear_zmat;
      ::QMtool::init_variables [namespace current];
      ::QMtool::molecule_export}
   $w.menubar.edit.menu add command -label "Reset QMtool" -command {::QMtool::reset_all}
   
   # Help
   menu $w.menubar.help.menu -tearoff no
   $w.menubar.help.menu add command -label "Website, Tutorial and FAQs" \
       -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"

   global vmd_initialize_structure
   trace add variable vmd_initialize_structure write ::QMtool::updateMolMenu
   grid [label $w.menubar.mollabel -text "Molecule:" -anchor w -width 8] -row 0 -column 1 -sticky ew -padx 4
   grid [ttk::combobox $w.menubar.inputmol -values $::QMtool::mollistmenu -width 20 -state readonly -justify left -textvariable ::QMtool::molselected] -row 0 -column 2 -sticky ew
   trace add variable ::QMtool::molselected write ::QMtool::newSelect

   ttk::style configure new.TNotebook -tabposition wn
   ttk::style configure new.TNotebook.Tab -width 15
   ttk::style configure new.TNotebook.Tab -anchor center
   font create customfont -size 100 -weight bold
   ttk::style configure New.TNotebook.Tab -font customfont
   grid [ttk::notebook $w.fp -style new.TNotebook -width 430] -row 1 -column 0 -sticky nsew -pady 2 -padx 2
   grid columnconfigure $w.fp 0 -weight 1
   grid rowconfigure $w.fp 0 -weight 1
  
   frame $w.fp.molecules
   #frame $w.fp.info
   frame $w.fp.input
   frame $w.fp.analysis
   #frame $w.fp.help
   
   set fontarg "helvetica 20 bold"
   
   $w.fp add $w.fp.molecules -text "\nMolecule Info\n " -padding 5 -sticky news
   #$w.fp add $w.fp.info -text "\nInformation\n " -padding 5 -sticky news
   $w.fp add $w.fp.input -text "\nDefine QM Input\n" -padding 5 -sticky news
   $w.fp add $w.fp.analysis -text "\nAnalysis\n " -padding 5 -sticky news
   #$w.fp add $w.fp.help -text "\nHelp\n " -padding 5 -sticky news

   # Create manage molecules tab
   #createMoleculeTab $w.fp.molecules

   # Create info tab
   ::QMtool::createInfoTab $w.fp.molecules

   # Create input tab
   ::QMtool::setupQM $w.fp.input

   # Create analysis tab
   ::QMtool::createAnalysisTab $w.fp.analysis

   # Create help tab
   #createHelpTab $w.fp.help

   variable molid
   variable molidlist
   #puts "molid: $molid; molidlist $molidlist"
   if {$molid>=0} {
      if {[molinfo $::QMtool::molid get numframes]>0} {
	 ::QMtool::update_molidlist
	 ::QMtool::update_intcoorlist
      }
   }

}

proc ::QMtool::updateMolMenu {args} {
   variable mollistmenu
   # need to populate mollistmenu with molecules
        
   set mollistmenu {}
   foreach t [molinfo list] {
        lappend mollistmenu [molinfo $t get name]
   }
   # update values of combobox
   .qmtool.menubar.inputmol configure -values $::QMtool::mollistmenu

}

proc ::QMtool::createInfoTab {w} {

   ############## frame for file info #################
   #grid [labelframe $w.files -bd 2 -relief ridge -text "Info" -width 45] -column 0 -row 3
   grid [frame $w.files] -padx 5 -pady 5
   grid columnconfigure $w.files 0 -weight 1
   grid rowconfigure $w.files 0 -weight 1
      
   label $w.files.softwaretypelabel -text "Software: "
   label $w.files.softwaretypevar -textvariable ::QMtool::qmsofttype
   grid $w.files.softwaretypelabel -column 0 -row 0 -sticky w
   grid $w.files.softwaretypevar   -column 1 -row 0 -sticky w

   label $w.files.filetypelabel -text "Filetype: "
   label $w.files.filetypevar -textvariable ::QMtool::filetype
   grid $w.files.filetypelabel -column 0 -row 1 -sticky w
   grid $w.files.filetypevar   -column 1 -row 1 -sticky w

   #label $w.files.filenamelabel -text "Filename: "
   #label $w.files.filenamevar -textvariable ::QMtool::filename
   #grid $w.files.filenamelabel -column 0 -row 2 -sticky w
   #grid $w.files.filenamevar   -column 1 -row 2 -sticky w

   grid [label $w.files.filenamelabel -text "Orbital Information: "] -column 0 -row 2 -sticky w
   grid [label $w.files.filenamevar -textvariable ::QMtool::hasMOs] -column 1 -row 2 -sticky w

   #label $w.files.checklabel -text "Checkfile: "
   #label $w.files.checkvar -textvariable ::QMtool::checkfile
   #grid $w.files.checklabel -column 0 -row 3 -sticky w
   #grid $w.files.checkvar   -column 1 -row 3 -sticky w

   label $w.files.nproclabel -text "Nproc: "
   label $w.files.nprocvar -textvariable ::QMtool::nproc
   grid $w.files.nproclabel -column 0 -row 3 -sticky w
   grid $w.files.nprocvar   -column 1 -row 3 -sticky w

   label $w.files.memorylabel -text "Memory: "
   label $w.files.memoryvar -textvariable ::QMtool::memory
   grid $w.files.memorylabel -column 0 -row 4 -sticky w
   grid $w.files.memoryvar   -column 1 -row 4 -sticky w

   label $w.files.routelabel -text "Route: "
   label $w.files.routevar -wraplength 10c -justify left -textvariable ::QMtool::route
   grid $w.files.routelabel  -column 0 -row 5 -sticky wn
   grid $w.files.routevar    -column 1 -row 5 -sticky wn

   label $w.files.titlelabel -text "Title: "
   label $w.files.titlevar -wraplength 10c -justify left -textvariable ::QMtool::title
   grid $w.files.titlelabel -column 0 -row 6 -sticky wn
   grid $w.files.titlevar   -column 1 -row 6 -sticky wn

   label $w.files.chargelabel -text "Total charge: "
   label $w.files.chargevar -textvariable ::QMtool::totalcharge
   grid $w.files.chargelabel -column 0 -row 7 -sticky w
   grid $w.files.chargevar   -column 1 -row 7 -sticky w

   label $w.files.multiplabel -text "Multiplicity: "
   label $w.files.multipvar -textvariable ::QMtool::multiplicity
   grid $w.files.multiplabel -column 0 -row 8 -sticky w
   grid $w.files.multipvar   -column 1 -row 8 -sticky w

   label $w.files.nimaglabel -text "Imaginary freq: "
   label $w.files.nimagvar -textvariable ::QMtool::nimag
   grid $w.files.nimaglabel -column 0 -row 9 -sticky w
   grid $w.files.nimagvar   -column 1 -row 9 -sticky w
 
}

#################################################################
# NOT BEING USED AT THE MOMENT
#################################################################


##########################################
# Create manage molecule frame
##########################################
proc createMoleculeTab {w} {
   variable selectcolor
   variable mollisttable
   variable molselected

   #buttons
   #load top molecule
   grid [frame $w.main] -column 0 -row 0 -sticky w -padx 10 -pady 10

 
   #grid [button $w.main.useTopMol -text "Load molecule" -command {set newmolid [molinfo top]; if {$newmolid!=-1} {::QMtool::use_vmd_molecule $newmolid}}] -column 0 -row 0 -sticky news -pady 2 -padx 2
   #new from molefacture 
   #grid [button $w.main.useMolefature -text "Edit molecule" -command {::QMtool::molefactureStart}] -column 1 -row 0 -sticky news -pady 2 -padx 2
   #delete molecule
   #grid [button $w.main.delMolecule -text "Delete molecule" -command {::QMtool::molecule_delete}] -column 2 -row 0 -sticky news -pady 2 -padx 2
   #clear selection
   #grid [button $w.main.clearSelection -text "Clear selection" -command {}] -column 0 -row 1 -sticky news -pady 2 -padx 2
   #reset molecule data

   #molecule list
   ############## frame for molecule list #################
   #grid [labelframe $w.mol -bd 2 -relief ridge -text "Molecule list" -padx 1m -pady 1m] -column 0 -row 2 -sticky w
   #grid [frame $w.mol.list]
   #grid [scrollbar $w.mol.list.scroll -orient vertical -command "$w.mol.list.list yview"] -column 1 -row 0 -sticky ens
   #grid [listbox $w.mol.list.list -activestyle dotbox -yscroll "$w.mol.list.scroll set"  \
      -width 43 -height 5  -selectmode browse -listvariable ::QMtool::molnamelist] -column 0 -row 0 -sticky w
   # -setgrid 1  -font {tkFixed 9} -selectbackground $selectcolor
   # This will be executed when a new molecule is selected:   
   #bind $w.mol.list.list <<ListboxSelect>> {}
      #set a [lindex [.qmtool.fp.molecules.mol.list.list get [.qmtool.fp.molecules.mol.list.list curselection]] 0]
      #::QMtool::molecule_select $a
   #{}
   
      ############## frame for molecule list #################
   grid [frame $w.mol] -column 0 -row 2 -sticky w -padx 25 -pady 5
      
   set fro2 $w.mol

   # NEW table for normalmode, from qwikmd
   option add *Tablelist.activeStyle       frame
   option add *Tablelist.background        gray98
   option add *Tablelist.stripeBackground  #e0e8f0
   option add *Tablelist.setGrid           no
   option add *Tablelist.movableColumns    no
  
   tablelist::tablelist $fro2.tb -width 40 -columns { 
                0 "ID"  center 
                0 "File Name"  center
                0 "Software" center
                } -yscrollcommand [list $fro2.scr1 set] -showseparators 0 -labelrelief groove  -labelbd 1 -selectbackground cyan \
                -selectforeground black -foreground black -background white -state normal -selectmode browse -stretch "0 1 2" -stripebackgroun white -height 5
   
   #-editstartcommand QWIKMD::cellStartEditPtcl -editendcommand QWIKMD::cellEndEditPtcl -forceeditendcommand true -editselectedonly true

   grid $fro2.tb -row 0 -column 0 -sticky news

   #$fro2.tb columnconfigure 0 -width 0
    ##Scroll_BAr V
   grid [scrollbar $fro2.scr1 -orient vertical -command [list $fro2.tb  yview]] -row 0 -column 1  -sticky ens

   set mollisttable $fro2.tb

   bind $fro2.tb <<TablelistSelect>>  {
      set a [lindex [.qmtool.fp.molecules.mol.tb get [.qmtool.fp.molecules.mol.tb curselection]] 0]
      ::QMtool::molecule_select $a
   }  

   # need to populate the table
   for {set i 0} {$i < [llength $::QMtool::molnamelist]} {incr i} {$mollisttable insert end [lindex $::QMtool::molnamelist $i]}

   variable statuscolor
   #grid [frame $w.status] -column 0 -row 2 -sticky w
   #grid [label $w.status.text -textvariable ::QMtool::statustext -fg $statuscolor] -column 0 -row 0 -sticky w
   
   #createInfoTab $w

}

proc createHelpTab {f} {
   #link to website
   
   grid [labelframe $f.faq -bd 2 -relief ridge -text "Links to Frequenly Ask Questions" -width 400]
   grid [frame $f.faq.main]
   #list of frequenly tasks, with links to screen capture videos in website
   #  I)   How to load a QM file in VMD
   grid [button $f.faq.main.q1 -text "How to load a QM output molecule?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 0 -sticky w -pady 2 -padx 2
   #  II)  How to Edit a structure
   grid [button $f.faq.main.q2 -text "How to edit the molecular structure?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 1 -sticky w -pady 2 -padx 2
   #  III) How to create an input
   grid [button $f.faq.main.q3 -text "How to create a simulation input file?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 2 -sticky w -pady 2 -padx 2
   #  IV)  How to plot energies
   grid [button $f.faq.main.q4 -text "How to plot geometry optimization energies?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 3 -sticky w -pady 2 -padx 2
   #  V)   How to show normal modes
   grid [button $f.faq.main.q5 -text "How to show normal mode analysis?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 4 -sticky w -pady 2 -padx 2
   #  VI)  How to plot molecular orbitals
   grid [button $f.faq.main.q6 -text "How to show molecular orbitals?" \
      -command "vmd_open_url [string trimright [vmdinfo www] /]/plugins/qmtool"] -column 0 -row 5 -sticky w -pady 2 -padx 2
   

}