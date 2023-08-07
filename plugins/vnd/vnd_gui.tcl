##
## VND -- Visual Neuronal Dynamics graphical interface
##
## $Id: vnd_gui.tcl,v 1.6 2022/04/25 17:24:22 mariano Exp $
##
##
## Home Page
## ---------
##   http://www.ks.uiuc.edu/Research/vnd/
##
##

package require Tk
package require tablelist

package provide vnd 0.5

namespace eval ::NeuronVND:: {

    variable modellist 0
    variable modelselected ""
    variable repselected 0
    variable selRep ""
    variable styleRep ""
    variable materialRep Opaque
    variable colorRep Type
    variable colorID ""
    variable showRep true
    variable sphereScale 3
    variable sphereRes 5
    variable proxyantialias ""
    variable proxydepthcueing ""
    variable proxyfps ""
    variable proxylight0 on 
    variable proxylight1 off 
    variable proxylight2 off 
    variable proxylight3 off
    variable renderMethod snapshot
    variable renderFile "vmdscene" 
    variable objList ""
    variable colorObj white
    variable historyCalls ""
    variable mouseMode R
    variable objMouse 0
    # variables for object management
    variable movex 10
    variable movey 10
    variable movez 10
    variable aggoffset {0 0 0}
    variable rotarx 10
    variable rotary 10
    variable rotarz 10
    variable aggrot [transidentity]
    # testing toplevel
    variable topGui ".neuron"
    variable bindTop 0
    # auxiliary variable for example rep checkbutton
    variable exampleRep 0
    variable exampleRepID -1
    # variable to use in connectivity GUI
    variable listOfRepsForConnect ""

    proc initialize {} {
        global env
        #########################
        variable listmodels
        set listmodels(-1) ""
        set listmodels(0,name) ""
        variable indexmodel 0

        # source neuro from barry
        #source /Projects/barryi/vmd/scripts/neuro_viz/neuro_read.tcl

    }
    initialize
}

source [file join $env(VNDPLUGINDIR) vnd_read.tcl]

proc ::NeuronVND::resizeGUI {w} {
  # taken from fftk
  update idletasks
  regexp {([0-9]+)x[0-9]+[\+\-]+[0-9]+[\+\-]+[0-9]+} [wm geometry $w] all dimW
  set dimH [winfo reqheight $w]
  #set dimH [expr {$dimH + 10}]
  set dimW [winfo reqwidth $w]
  #set dimW [expr {$dimW + 5}]
  wm geometry $w [format "%ix%i" $dimW $dimH]
  update idletasks

}

proc ::NeuronVND::neuronRep { } {
    variable repselected
    variable selRep
    variable styleRep
    variable materialRep
    variable colorRep
    variable colorID
    variable sphereScale
    variable sphereRes

    set w .neuron.fp.systems.rep
    #wm title $w "Representations"
    #wm resizable $w 1 1
    #set width 288 ;# in pixels
    #set height 160 ;# in pixels
    #wm geometry $w ${width}x${height}+797+747   
    grid columnconfigure $w 0 -weight 1

    grid [labelframe $w.main -text "Representations" -labelanchor n] -row 0 -column 0 -sticky news
    grid columnconfigure $w.main 0 -weight 1
    #grid [ttk::combobox $w.main.modelsel.inp -width 37 -background white -values $::NeuronVND::modellist -state readonly -justify left -textvariable ::NeuronVND::modelselected] -row 0 -column 1 -sticky ew -padx 1
    #bind $w.main.modelsel.inp <<ComboboxSelected>> {set text [%W get]; %W selection clear}

    grid [frame $w.main.rep] -row 1 -column 0 -sticky news -padx 2 -pady 2
    grid [button $w.main.rep.add -text "Create Rep" -command {::NeuronVND::createRepArgs}] -row 0 -column 0 -sticky n;#ews
    grid [button $w.main.rep.show -text "Show / Hide" -command {::NeuronVND::showHideRep}] -row 0 -column 1 -sticky n;#ews
    grid [button $w.main.rep.del -text "Delete Rep" -command {::NeuronVND::delRep}] -row 0 -column 2 -sticky n;#e
    #grid columnconfigure $w.main.rep 2 -weight 1

    grid [frame $w.main.table] -row 2 -column 0 -sticky news -padx 4 -pady 2
    grid columnconfigure $w.main.table 0 -weight 1
    grid [tablelist::tablelist $w.main.table.tb -columns {
        0 "Style" 
        0 "Color"
        0 "Selection"
        } \
        -yscrollcommand [list $w.main.table.scr1 set] \
        -stretch all -background white -stretch 2 -height 6 -width 100 -exportselection false]
    
    ##Scroll_BAr V
    grid [scrollbar $w.main.table.scr1 -orient vertical -command [list $w.main.table.tb yview]] -row 0 -column 1  -sticky ens

    $w.main.table.tb columnconfigure 0 -width 10
    #$w.main.table.tb columnconfigure 2 -width 15
    #$w.main.table.tb columnconfigure 0 -width 0 -editable true -editwindow ttk::checkbutton

    bind $w.main.table.tb <<TablelistSelect>>  {
      set ::NeuronVND::repselected [%W curselection]  
      ::NeuronVND::updateRepMenu
    }

    grid [labelframe $w.main.sel -text "Selected Neurons" -labelanchor n -borderwidth 0] -row 3 -column 0 -sticky news -padx 0 -pady 2
    grid [entry $w.main.sel.entry -textvariable ::NeuronVND::selRep -width 100] -row 0 -column 0 -sticky news -padx 0
    bind $w.main.sel.entry <Return> {
      ::NeuronVND::editRep sel
      return
    }
    

    grid [frame $w.main.def] -row 4 -column 0 -sticky news -padx 5 -pady 1
    grid [label $w.main.def.colorlbl -text "Coloring Method" -anchor c] -row 0 -column 0
    grid [ttk::combobox $w.main.def.colorcb -width 12 -values {"Type" "Color"} -textvariable ::NeuronVND::colorRep -state readonly] -row 1 -column 0
    # button option for color not being used
    button $w.main.def.colorid -background white -width 1 -command {
      set auxcolor [tk_chooseColor -initialcolor $::NeuronVND::colorID -title "Choose color"]
      if {$auxcolor != ""} {
          set ::NeuronVND::colorID $auxcolor
          .neuron.fp.systems.rep.main.def.colorid configure -background $auxcolor}
    }

    grid [label $w.main.def.coloridlb -width 9] -row 1 -column 1 -sticky news
    ttk::combobox $w.main.def.coloridcb -width 7 -values [colorinfo colors] -textvariable ::NeuronVND::colorID -state readonly

    bind $w.main.def.colorcb <<ComboboxSelected>> {
        set text [%W get]
        switch $text {
            "Color" {
                grid .neuron.fp.systems.rep.main.def.coloridcb -row 1 -column 1 -sticky news
            }
            "default" {
                grid remove .neuron.fp.systems.rep.main.def.coloridcb
                set ::NeuronVND::colorID Type
                ::NeuronVND::editRep color
            }
        }
        #::NeuronVND::editRep color
        %W selection clear
    }

    bind $w.main.def.coloridcb <<ComboboxSelected>> {
        set text [%W get]
        ::NeuronVND::editRep color
        %W selection clear
    }    
    
    
    grid [label $w.main.def.matlbl -text "Material" -width 10 -anchor c] -row 0 -column 3
    set materiallist {"Opaque" "Transparent" "BrushedMetal" "Diffuse" "Ghost" "Glass1" "Glass2" "Glass3" "Glossy" "HardPlastic" "MetallicPastel" "Steel" \
        "Translucent" "Edgy" "EdgyShiny" "EdgyGlass" "Goodsell" "AOShiny" "AOChalky" "AOEdgy" "BlownGlass" "GlassBubble" "RTChrome"}
    grid [ttk::combobox $w.main.def.matcb -width 12 -values $materiallist -textvariable ::NeuronVND::materialRep -state readonly] -row 1 -column 3
    bind $w.main.def.matcb <<ComboboxSelected>> {
        set text [%W get]
        ::NeuronVND::editRep material
        %W selection clear
    }

    grid [label $w.main.def.stylbl -text "Style"] -row 2 -column 0
    grid [ttk::combobox $w.main.def.stycb -width 12 -values {"soma" "morphology"} -textvariable ::NeuronVND::styleRep -state readonly] -row 3 -column 0
    bind $w.main.def.stycb <<ComboboxSelected>> {
        set text [%W get]
        switch $text {
            "soma" {
                #.neuron.fp.systems.rep.main.arg.en1 configure -state normal
                #.neuron.fp.systems.rep.main.arg.en2 configure -state normal
            }
            "morphology" {
                #.neuron.fp.systems.rep.main.arg.en1 configure -state disabled
                #.neuron.fp.systems.rep.main.arg.en2 configure -state disabled
            }
        }
        ::NeuronVND::editRep style
        %W selection clear
    }

    grid [frame $w.main.arg] -row 5 -column 0 -sticky e -padx 2 -pady 2

    #grid [label $w.main.arg.lb1 -text "Sphere Scale" -anchor e] -row 0 -column 0 -sticky news
    #grid [spinbox $w.main.arg.en1 -width 3 -increment 1 -from 1 -to 10 -textvariable ::NeuronVND::sphereScale -background white -command {::NeuronVND::editRep sphere}]  -row 0 -column 1 -padx 2 -sticky w
    #grid [entry $w.main.arg.en1 -text "5" -width 5] -row 0 -column 1 -sticky news
    #grid [label $w.main.arg.lb2 -text "Sphere Resolution"] -row 1 -column 0 -sticky news
    #grid [spinbox $w.main.arg.en2 -width 3 -increment 5 -from 5 -to 30 -textvariable ::NeuronVND::sphereRes -background white -command {::NeuronVND::editRep sphere}]  -row 1 -column 1 -padx 2 -sticky w

    #NeuronVND::resizeGUI $w

}

proc ::NeuronVND::createPages { orientation } {
    variable listmodels
    set w .neuron
    if {[winfo exists .neuron.fp]} { destroy .neuron.fp }
    ttk::style configure new.TNotebook -tabposition $orientation
    ttk::style configure new.TNotebook.Tab -width 12
    ttk::style configure new.TNotebook.Tab -anchor center
    #font create customfont2 -size 100 -weight bold
    ttk::style configure New.TNotebook.Tab -font customfont2
    grid [ttk::notebook $w.fp -style new.TNotebook -width 420] -row 1 -column 0 -sticky nsew -pady 2 -padx 2
    grid columnconfigure $w.fp 0 -weight 1
    grid rowconfigure $w.fp 0 -weight 1

    frame $w.fp.systems
    #frame $w.fp.systems.rep
    frame $w.fp.navigation

    set fontarg "helvetica 20 bold"

    if {$orientation == "wn"} {
        set text1 "\nMain\n"
        set text2 "\nGraphics\n"
        set text3 "\nNavigation\n"
        set width 350
        set height 300
    } else {
        set text1 "Main"
        set text2 "Graphics"
        set text3 "Navigation"
        set width 532
        set height 352
    }

    $w.fp add $w.fp.systems -text $text1 -padding 2 -sticky news
    #$w.fp add $w.fp.systems.rep -text $text2 -padding 2 -sticky news
    $w.fp add $w.fp.navigation -text $text3 -padding 2 -sticky news

    grid [labelframe $w.fp.systems.main -text "Systems" -labelanchor n] -row 0 -column 0 -sticky news
    grid [tablelist::tablelist $w.fp.systems.main.tb -columns {
        0 "ID" 
        0 "T"
        0 "D"
        0 "Name"
        0 "Neurons"
        } \
        -yscrollcommand [list $w.fp.systems.main.scr1 set] \
        -stretch all -background white -stretch all -height 4 -width 40 -exportselection false]  

    ##Scroll_BAr V
    grid [scrollbar $w.fp.systems.main.scr1 -orient vertical -command [list $w.fp.systems.main.tb yview]] -row 0 -column 1  -sticky ens

    #$w.main.tb insert end [list "0" "T" "D" "V1" [::neuro::cmd_query num_neurons] ""]

    #$w.main.tb insert end [list "0" "T" "D" "event30K" "30.000" "100"]
    #$w.main.tb insert end [list "1" "" "D" "test500K" "500.000" "0"]

    # Testing add dropside menu
    # First add thin vertical button to increase window with
    grid [label $w.fp.systems.sidebutton1 -text "<" -width 1] -row 0 -rowspan 2 -column 1 -sticky ens
    label $w.fp.systems.sidebutton2 -text ">" -width 1

    # set mouse click bindings to expand/contract window
    bind $w.fp.systems.sidebutton1 <Button-1> {
        grid remove .neuron.fp.systems.sidebutton1
        grid .neuron.fp.systems.sidebutton2 -row 0 -rowspan 2 -column 1 -sticky ens
        # resize to hide representations
        wm geometry .neuron 281x352
    }

    bind $w.fp.systems.sidebutton2 <Button-1> {
        grid remove .neuron.fp.systems.sidebutton2
        grid .neuron.fp.systems.sidebutton1
        # resize to show representations
        wm geometry .neuron 532x352
    }

    # frame for representations
    grid [frame $w.fp.systems.rep] -row 0 -rowspan 2 -column 2 -sticky news
    grid columnconfigure $w.fp.systems 2 -weight 1

    # Add info frame
    grid [labelframe $w.fp.systems.info -text "Information" -labelanchor n] -row 1 -column 0 -sticky news
    set rowid 0
    grid [label $w.fp.systems.info.header1 -text "Tree view"] -row $rowid -column 1 -sticky news -padx 1
    #grid [label $w.fp.systems.info.header2 -text "Query"] -row $rowid -column 2 -sticky news -padx 1
    incr rowid
if {0} {
    grid [tablelist::tablelist $w.fp.systems.info.tb2 -columns {
        0 "Results" 
        } \
        -yscrollcommand [list $w.fp.systems.info.scr2 set] \
        -stretch all -background white -stretch all -height 9 -width 12 -showlabels 0] -row $rowid -column 2 -padx 1

    ##Scroll_BAr V
    grid [scrollbar $w.fp.systems.info.scr2 -orient vertical -command [list $w.fp.systems.info.tb2 yview]] -row $rowid -column 3  -sticky ens -padx 1
}  

    # Testing treeview
    grid [ttk::treeview $w.fp.systems.info.tv -show tree -height 6 -yscrollcommand [list $w.fp.systems.info.scr1 set]] -row $rowid -column 1 -padx 1
    ##Scroll_BAr V
    grid [scrollbar $w.fp.systems.info.scr1 -orient vertical -command [list $w.fp.systems.info.tv yview]] -row $rowid -column 0  -sticky wens -padx 1

    set tv $w.fp.systems.info.tv
    #$tv heading #0 -text "Model"
    $tv column #0 -width 200

    incr rowid
    grid [button $w.fp.systems.info.button1 -text "Create Representation" -command {::NeuronVND::createExampleRep [.neuron.fp.systems.info.tv selection]}] -row $rowid -column 1 -sticky w

    # auto rep creation replaced by button
    if {0} {  
    bind $w.fp.systems.info.tv <<TreeviewSelect>> {
        set sel [.neuron.fp.systems.info.tv selection]
        if {$::NeuronVND::exampleRep} {
            ::NeuronVND::createExampleRep $sel
        }
    }
    }

    #$tv configure -columns "results"
    #$tv heading results -text "Results"
    #$tv column results -width 125

    if {$listmodels(0,name) != ""} {
    .neuron.fp.systems.main.tb insert end [list "0" "T" "D" $listmodels(0,name) $listmodels(0,neurons)]
    }

    grid [frame $w.fp.navigation.main] -row 0 -column 0 -sticky news
       
    grid [labelframe $w.fp.navigation.main.mmode -text "Mouse Mode" -labelanchor n -width 20] -row 0 -column 0 -sticky news
    grid [radiobutton $w.fp.navigation.main.mmode.rot -text "Rotate (R)" -variable ::NeuronVND::mouseMode -value R] -row 0 -column 0 -sticky nws
    grid [radiobutton $w.fp.navigation.main.mmode.trans -text "Translate (T)" -variable ::NeuronVND::mouseMode -value T] -row 1 -column 0 -sticky nws
    grid [radiobutton $w.fp.navigation.main.mmode.scale -text "Scale (S)" -variable ::NeuronVND::mouseMode -value S] -row 2 -column 0 -sticky nws

    grid [labelframe $w.fp.navigation.main.obj -text "Object Management" -labelanchor n] -row 0 -column 1 -sticky nes
    grid [label $w.fp.navigation.main.obj.sel -text "Select:"] -row 0 -column 0 -sticky news
    grid [ttk::combobox $w.fp.navigation.main.obj.cb -width 30 -values "" -state readonly] -row 0 -column 1 -columnspan 4
    grid [radiobutton $w.fp.navigation.main.obj.move -text "Use Mouse Mode to arrange object" -variable ::NeuronVND::objMouse -value 1 -command {
        # need to fix all except for objid
        foreach m [molinfo list] {mol fix $m}
        mol free [lindex $::NeuronVND::objList 0 0]
    }] -row 1 -column 0 -sticky nws -columnspan 5
    grid [radiobutton $w.fp.navigation.main.obj.move2 -text "Use buttons to arrange object" -variable ::NeuronVND::objMouse -value 0 -command {
        foreach m [molinfo list] {mol free $m}    
    }] -row 2 -column 0 -columnspan 5 -sticky nws
    grid [label $w.fp.navigation.main.obj.movex -text "Translate in X:"] -row 3 -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.movex1 -text "-" -command {::NeuronVND::moveGraphs x neg}] -row 3 -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.movex2 -text "-" -command {::NeuronVND::moveGraphs x -100}] -row 3 -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.movex2 -width 3 -increment 10 -from 10 -to 1000 -textvariable ::NeuronVND::movex -background white]  -row 3 -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.movex3 -text "+" -command {::NeuronVND::moveGraphs x pos}] -row 3 -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.movex4 -text "+++" -command {::NeuronVND::moveGraphs x 500}] -row 3 -column 4 -sticky news
    grid [label $w.fp.navigation.main.obj.movey -text "Translate in Y:"] -row 4 -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.movey1 -text "-" -command {::NeuronVND::moveGraphs y neg}] -row 4 -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.movey2 -text "-" -command {::NeuronVND::moveGraphs y -100}] -row 4 -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.movey2 -width 3 -increment 10 -from 10 -to 1000 -textvariable ::NeuronVND::movey -background white]  -row 4 -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.movey3 -text "+" -command {::NeuronVND::moveGraphs y pos}] -row 4 -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.movey4 -text "+++" -command {::NeuronVND::moveGraphs y 500}] -row 4 -column 4 -sticky news
    grid [label $w.fp.navigation.main.obj.movez -text "Translate in Z:"] -row 5 -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.movez1 -text "-" -command {::NeuronVND::moveGraphs z neg}] -row 5 -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.movez2 -text "-" -command {::NeuronVND::moveGraphs z -100}] -row 5 -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.movez2 -width 3 -increment 10 -from 10 -to 1000 -textvariable ::NeuronVND::movez -background white]  -row 5 -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.movez3 -text "+" -command {::NeuronVND::moveGraphs z pos}] -row 5 -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.movez4 -text "+++" -command {::NeuronVND::moveGraphs z 500}] -row 5 -column 4 -sticky news
    
    grid [ttk::separator $w.fp.navigation.main.obj.sep1] -row 6 -column 0 -sticky news -columnspan 5 -pady 2 -padx 2
    set row 7
    grid [label $w.fp.navigation.main.obj.rotx -text "Rotate around X:"] -row $row -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.rotx1 -text "-" -command {::NeuronVND::rotGraphs x neg}] -row $row -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.rotx2 -text "-" -command {::NeuronVND::rotGraphs x -15}] -row $row -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.rotx2 -width 3 -increment 5 -from 5 -to 180 -textvariable ::NeuronVND::rotarx -background white]  -row $row -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.rotx3 -text "+" -command {::NeuronVND::rotGraphs x pos}] -row $row -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.rotx4 -text "+++" -command {::NeuronVND::rotGraphs x 45}] -row $row -column 4 -sticky news
    incr row
    grid [label $w.fp.navigation.main.obj.roty -text "Rotate around Y:"] -row $row -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.roty1 -text "-" -command {::NeuronVND::rotGraphs y neg}] -row $row -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.roty2 -text "-" -command {::NeuronVND::rotGraphs y -15}] -row $row -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.roty2 -width 3 -increment 5 -from 5 -to 180 -textvariable ::NeuronVND::rotary -background white]  -row $row -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.roty3 -text "+" -command {::NeuronVND::rotGraphs y pos}] -row $row -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.roty4 -text "+++" -command {::NeuronVND::rotGraphs y 45}] -row $row -column 4 -sticky news
    incr row
    grid [label $w.fp.navigation.main.obj.rotz -text "Rotate around Z:"] -row $row -column 0 -sticky news
    grid [button $w.fp.navigation.main.obj.rotz1 -text "-" -command {::NeuronVND::rotGraphs z neg}] -row $row -column 1 -sticky news
    #grid [button $w.fp.navigation.main.obj.rotz2 -text "-" -command {::NeuronVND::rotGraphs z -15}] -row $row -column 2 -sticky news
    grid [spinbox $w.fp.navigation.main.obj.rotz2 -width 3 -increment 5 -from 5 -to 180 -textvariable ::NeuronVND::rotarz -background white]  -row $row -column 2 -sticky news
    grid [button $w.fp.navigation.main.obj.rotz3 -text "+" -command {::NeuronVND::rotGraphs z pos}] -row $row -column 3 -sticky news
    #grid [button $w.fp.navigation.main.obj.rotz4 -text "+++" -command {::NeuronVND::rotGraphs z 45}] -row $row -column 4 -sticky news
    incr row
    grid [ttk::separator $w.fp.navigation.main.obj.sep2] -row $row -column 0 -sticky news -columnspan 5 -pady 2 -padx 2
    incr row
    grid [label $w.fp.navigation.main.obj.color -text "Color:"] -row $row -column 0 -sticky news 
    grid [ttk::combobox $w.fp.navigation.main.obj.colorcb -width 7 -values [colorinfo colors] -textvariable ::NeuronVND::colorObj -state readonly] -row $row -column 1 -sticky news -columnspan 2

    bind $w.fp.navigation.main.obj.colorcb <<ComboboxSelected>> {
        set ::NeuronVND::colorObj [%W get]
        set objid [lindex $::NeuronVND::objList 0 0]
        if {$objid != ""} {
            graphics $objid replace 0
            graphics $objid color $::NeuronVND::colorObj
        }
        %W selection clear
    }
    #grid [button $w.scaleframe.left -image left -command { }] -row 0 -column 0
    #grid [entry $w.fp.navigation.main.entry -textvariable timeentry -width 5] -row 0 -column 0
    #grid [scale $w.fp.navigation.main.scale -from 0 -to 100 -orien horizontal -length 220 -variable timeentry -sliderlength 10 -command {}] -row 0 -column 1 -sticky news -columnspan 3
    #grid [label $w.fp.navigation.main.steplbl -text "Step" -width 5] -row 1 -column 0 -sticky news
    #grid [entry $w.fp.navigation.main.stepentry -textvariable stepval -width 5] -row 1 -column 1 -sticky news
    #grid [label $w.fp.navigation.main.speedlbl -text "Speed" -width 5] -row 1 -column 2 -sticky news
    #grid [scale $w.fp.navigation.main.speed -from 1 -to 5 -orien horizontal -length 20 -label "" -variable speedvar -sliderlength 10 -command {}] -row 1 -column 3 -sticky news

    ::NeuronVND::neuronRep
    ::NeuronVND::renderPage
    ::NeuronVND::connectGUI

    wm geometry $w ${width}x${height}

}

proc ::NeuronVND::renderPage { } {
    set w .neuron
    frame $w.fp.render
    $w.fp add $w.fp.render -text "Render" -padding 2 -sticky news

    grid [frame $w.fp.render.main] -row 0 -column 0 -sticky news
    set gr 0
    grid [label $w.fp.render.main.reslbl -text "Resolution:"] -row $gr -column 0 -sticky news -padx 1 -pady 2
    set reslist {"SD (480p)" "HD (720p)" "FullHD (1080p)" "QuadHD (1440p)" "2K (1080p)" "4K (2160p)" "8K (4320p)"}
    grid [ttk::combobox $w.fp.render.main.rescb -width 25 -values $reslist -state readonly] -row $gr -column 1 -sticky news -padx 1 -pady 2
    bind $w.fp.render.main.rescb <<ComboboxSelected>> {
        set text [%W get]
        switch $text {
            "SD (480p)" { display resize 640 480 }
            "HD (720p)" { display resize 1280 720 }
            "FullHD (1080p)" { display resize 1920 1080 }
            "QuadHD (1440p)" { display resize 2560 1440 }
            "2K (1080p)" { display resize 2048 1080 }
            "4K (2160p)" { display resize 3840 2160 }
            "8K (4320p)" { display resize 7680 4320 }
        }
        %W selection clear
    }
    incr gr
    
    grid [label $w.fp.render.main.renderlbl -text "Render using:"] -row $gr -column 0 -sticky news -padx 1 -pady 2
    set renderlist [render list]
    grid [ttk::combobox $w.fp.render.main.rendercb -width 25 -values $renderlist -state readonly] -row $gr -column 1 -sticky news -padx 1 -pady 2
    bind $w.fp.render.main.rendercb <<ComboboxSelected>> {
        set text [%W get]
        set ::NeuronVND::renderMethod $text
        %W selection clear
    }
    incr gr
    grid [label $w.fp.render.main.filelbl -text "Filename:"] -row $gr -column 0 -sticky news -padx 1 -pady 2
    grid [entry $w.fp.render.main.fileentry -textvariable ::NeuronVND::renderFile] -row $gr -column 1 -sticky news -padx 1 -pady 2
    incr gr
    grid [button $w.fp.render.main.renderbut -text "Start Rendering" -command {
        render $::NeuronVND::renderMethod $::NeuronVND::renderFile [render default $::NeuronVND::renderMethod]
    }] -row $gr -column 0 -columnspan 2 -sticky news -padx 1 -pady 2

}

proc ::NeuronVND::connectGUI { } {
    variable listOfRepsForConnect

    set w .neuron
    frame $w.fp.connect
    $w.fp add $w.fp.connect -text "Connectivity" -padding 2 -sticky news

    grid [frame $w.fp.connect.main] -row 0 -column 0 -sticky news
    set gr 0
    set aframe $w.fp.connect.main
    grid [label $aframe.title -text "Display connections (edges) between selection of neurons"] -row $gr
    incr gr
    grid [labelframe $aframe.lbl1 -text "Source" -labelanchor n ] -row $gr -column 0 -sticky news 
    grid [label $aframe.lbl1.sel -text "Select from existing rep:"] -row 0 -column 0 -sticky news
    grid [ttk::combobox $aframe.lbl1.cb -width 30 -values $listOfRepsForConnect -state readonly] -row 0 -column 1 -columnspan 4
    incr gr
    grid [labelframe $aframe.lbl2 -text "Target" -labelanchor n] -row $gr -column 0 -sticky news 
    grid [label $aframe.lbl2.sel -text "Select from existing rep:"] -row 0 -column 0 -sticky news
    grid [ttk::combobox $aframe.lbl2.cb -width 30 -values $listOfRepsForConnect -state readonly] -row 0 -column 1 -columnspan 4
    incr gr
    grid [frame $aframe.but] -row $gr -column 0
    grid [button $aframe.but.create -text "Create connection rep" -command {}] -row 0 -column 0
    grid [button $aframe.but.show -text "Show/Hide connection rep" -command {}] -row 0 -column 0
    grid [button $aframe.but.del -text "Delete connection rep" -command {}] -row 0 -column 2

    incr gr
    grid [frame $aframe.edgeslist] -row $gr -column 0
    grid [tablelist::tablelist $aframe.edgeslist.tb -columns {
        0 "Style" 
        0 "Source"
        0 "Target"
        } \
        -yscrollcommand [list $aframe.edgeslist.scr1 set] \
        -stretch all -background white -stretch 2 -height 6 -width 100 -exportselection false]
    
    ##Scroll_BAr V
    grid [scrollbar $aframe.edgeslist.scr1 -orient vertical -command [list $aframe.edgeslist.tb yview]] -row 0 -column 1  -sticky ens

    $aframe.edgeslist.tb columnconfigure 0 -width 10
}

proc ::NeuronVND::neuronGui { } {

   variable timeentry
   variable proxyantialias
   variable proxydepthcueing
   variable proxyfps
   variable proxylight0 
   variable proxylight1 
   variable proxylight2 
   variable proxylight3 

   set w [toplevel $::NeuronVND::topGui]
   wm title $w "Visual Neuronal Dynamics"
   wm resizable $w 1 1
   set width 532 ;# in pixels
   set height 300 ;# in pixels 290x200+782+454
   wm geometry $w ${width}x${height}+782+454
   grid columnconfigure $w 0 -weight 1
   grid columnconfigure $w 1 -weight 0
   grid rowconfigure $w 0 -weight 0
   grid rowconfigure $w 1 -weight 1

   wm protocol $::NeuronVND::topGui WM_DELETE_WINDOW ::NeuronVND::exit

   grid [frame $w.menubar -relief raised -bd 2] -row 0 -column 0 -sticky nswe -pady 2 -padx 2
   grid columnconfigure $w.menubar 4 -weight 1
   grid rowconfigure $w.menubar 0 -weight 1

   grid [menubutton $w.menubar.file -text "File" -width 5 -menu $w.menubar.file.menu] -row 0 -column 0 -sticky ew
   #grid [menubutton $w.menubar.system -text "System" -width 8 -menu $w.menubar.system.menu] -row 0 -column 1 -sticky ew
   grid [menubutton $w.menubar.display -text "Display" -width 8 -menu $w.menubar.display.menu] -row 0 -column 2 -sticky ew
   grid [menubutton $w.menubar.analysis -text "Analysis" -width 8 -menu $w.menubar.analysis.menu] -row 0 -column 3 -sticky ew
   grid [menubutton $w.menubar.help -text "Help" -width 5 -menu $w.menubar.help.menu] -row 0 -column 4 -sticky e
      
   # File
   menu $w.menubar.file.menu -tearoff no
   $w.menubar.file.menu add command -label "Open File" -command { 
        set cfgfile [tk_getOpenFile -initialdir "." -title "Choose config file"]
        if {$cfgfile != ""} {::NeuronVND::loadFiles $cfgfile}
   }
   $w.menubar.file.menu add command -label "Open File with Edges" -command { 
        set cfgfile [tk_getOpenFile -initialdir "." -title "Choose config file"]
        if {$cfgfile != ""} {::NeuronVND::loadFiles $cfgfile true true}
   }
   $w.menubar.file.menu add command -label "Add Object" -command { 
        set file [tk_getOpenFile -initialdir "." -title "Choose object file"]
        if {$file != ""} {::NeuronVND::loadObject $file}
   }
   $w.menubar.file.menu add separator
   $w.menubar.file.menu add command -label "Load Visualization State" -command {::NeuronVND::visState load}
   $w.menubar.file.menu add command -label "Save Visualization State" -command {::NeuronVND::visState save}

   $w.menubar.file.menu add separator
   $w.menubar.file.menu add command -label "Quit" -command {::NeuronVND::exit}

   #$w.menubar.file.menu add command -label "Write Input File" -command { }

   # System
   #menu $w.menubar.system.menu -tearoff no
   #$w.menubar.system.menu add command -label "System Information" -command { ::NeuronVND::neuronInfo } -state disabled
   #$w.menubar.system.menu add command -label "Representations" -command { ::NeuronVND::neuronRep }

    # Display
   menu $w.menubar.display.menu -tearoff no
   #menu $w.menubar.display.menu.orient -tearoff no -title "Menu"
   #$w.menubar.display.menu add cascade -label "Menu" -menu $w.menubar.display.menu.orient
   #$w.menubar.display.menu.orient add radiobutton -label "Horizontal" -variable orient -value "nw" -command { ::NeuronVND::createPages nw }
   #$w.menubar.display.menu.orient add radiobutton -label "Vertical" -variable orient -value "wn" -command { ::NeuronVND::createPages wn }
   
   $w.menubar.display.menu add command -label "Reset View" -command { display resetview; scale by 0.021; translate by 0.0 1.0 0.0}
   $w.menubar.display.menu add separator
   $w.menubar.display.menu add radiobutton -label "Perspective" -variable perps -value on -command {display projection Perspective}
   $w.menubar.display.menu add radiobutton -label "Orthographic" -variable perps -value off -command {display projection Orthographic}
   set perps on
   $w.menubar.display.menu add separator   
   $w.menubar.display.menu add checkbutton -label "Antialiasing" -variable ::NeuronVND::proxyantialias -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxyantialias {
           "on"  { display antialias on }
           "off" { display antialias off }
       }
   }
   $w.menubar.display.menu add checkbutton -label "Depth Cueing" -variable ::NeuronVND::proxydepthcueing -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxydepthcueing {
           "on"  { display depthcue on }
           "off" { display depthcue off }
       }
   }
   $w.menubar.display.menu add checkbutton -label "FPS Indicator" -variable ::NeuronVND::proxyfps -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxyfps {
           "on"  { display fps on }
           "off" { display fps off }
       }
   }
   $w.menubar.display.menu add separator   
   $w.menubar.display.menu add checkbutton -label "Light 0" -variable ::NeuronVND::proxylight0 -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxylight0 {
           "on"  { light 0 on }
           "off" { light 0 off }
       }
   }
   $w.menubar.display.menu add checkbutton -label "Light 1" -variable ::NeuronVND::proxylight1 -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxylight1 {
           "on"  { light 1 on }
           "off" { light 1 off }
       }
   }
   $w.menubar.display.menu add checkbutton -label "Light 2" -variable ::NeuronVND::proxylight2 -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxylight2 {
           "on"  { light 2 on }
           "off" { light 2 off }
       }
   }
   $w.menubar.display.menu add checkbutton -label "Light 3" -variable ::NeuronVND::proxylight3 -onvalue on -offvalue off -command { 
       switch $::NeuronVND::proxylight3 {
           "on"  { light 3 on }
           "off" { light 3 off }
       }
   }      
   $w.menubar.display.menu add separator   
   menu $w.menubar.display.menu.axes -tearoff no -title "Axes"
   $w.menubar.display.menu add cascade -label "Axes" -menu $w.menubar.display.menu.axes
   $w.menubar.display.menu.axes add radiobutton -label "Off" -variable axes -value off -command { axes location Off }
   $w.menubar.display.menu.axes add radiobutton -label "Origin" -variable axes -value origin -command { axes location Origin }
   $w.menubar.display.menu.axes add radiobutton -label "Lower Left" -variable axes -value lowerleft -command { axes location LowerLeft }
   $w.menubar.display.menu.axes add radiobutton -label "Lower Right" -variable axes -value lowerright -command { axes location LowerRight }
   $w.menubar.display.menu.axes add radiobutton -label "Upper Left" -variable axes -value upperleft -command { axes location UpperLeft }
   $w.menubar.display.menu.axes add radiobutton -label "Upper Right" -variable axes -value upperright -command { axes location UpperRight }
   
   menu $w.menubar.display.menu.background -tearoff no -title "Background"
   $w.menubar.display.menu add cascade -label "Background" -menu $w.menubar.display.menu.background
   $w.menubar.display.menu.background add radiobutton -label "Solid Color" -variable bgsolid -value on -command { display backgroundgradient off }
   $w.menubar.display.menu.background add radiobutton -label "Gradient" -variable bgsolid -value off -command { display backgroundgradient on }
   $w.menubar.display.menu add separator
   menu $w.menubar.display.menu.rendermode -tearoff no -title "Render Mode"
   $w.menubar.display.menu add cascade -label "Render Mode" -menu $w.menubar.display.menu.rendermode
   $w.menubar.display.menu.rendermode add radiobutton -label "Normal" -variable render -value normal -command { display rendermode Normal }
   $w.menubar.display.menu.rendermode add radiobutton -label "GLSL" -variable render -value glsl -command { display rendermode GLSL }
   $w.menubar.display.menu.rendermode add radiobutton -label "Tachyon RTX RTRT" -variable render -value rtrt -command { display rendermode "Tachyon RTX RTRT" }
   $w.menubar.display.menu.rendermode add radiobutton -label "Acrobat3D" -variable render -value a3D -command { display rendermode Acrobat3D }
   $w.menubar.display.menu add separator
   $w.menubar.display.menu add command -label "Display Settings" -command { menu display off; menu display on }

    # Analysis
   menu $w.menubar.analysis.menu -tearoff no
   $w.menubar.analysis.menu add command -label "Timeline Analysis" -command { neuronTimeline } -state disabled
   
   # Help
   menu $w.menubar.help.menu -tearoff no
   $w.menubar.help.menu add command -label "Website, Tutorial and FAQs" \
       -command "vmd_open_url https://www.ks.uiuc.edu/Research/vnd/"
   $w.menubar.help.menu add checkbutton -label "Debug Mode" -variable ::neuro::debugMode -onvalue 1 -offvalue 0 -command { 
       puts "debugMode $::neuro::debugMode"
   }      

    ::NeuronVND::createPages nw

}

::NeuronVND::neuronGui




proc ::NeuronVND::neuronInfo { } {

   set w [toplevel ".neuron.info"]
   wm title $w "System Information"
   wm resizable $w 1 1
   set width 290 ;# in pixels
   set height 100 ;# in pixels
   wm geometry $w ${width}x${height}

   grid [frame $w.main] -row 0 -column 0 -sticky news

   grid [ttk::frame $w.main.t1] -row 0 -column 0 -sticky nswe -padx 4 -columnspan 8

   grid columnconfigure $w.main.t1 0 -weight 1
   grid rowconfigure $w.main.t1 0 -weight 1


   #grid columnconfigure $w.main 0 -weight 1
   #grid rowconfigure $w.main 0 -weight 1

   #option add *Tablelist.activeStyle       frame
   
   set fro2 $w.main.t1

   option add *Tablelist.movableColumns    no
   option add *Tablelist.labelCommand      tablelist::sortByColumn

       tablelist::tablelist $fro2.tb -columns {\
           0 "Type" center
           0 "Number" center
           0 "Events" center
           0 "Notes" center
       }\
       -yscrollcommand [list $fro2.scr1 set] \
               -showseparators 0 -labelrelief groove  -labelbd 1 -selectforeground black\
               -foreground black -background white -width 45 -height 6 -state normal -selectmode extended -stretch all -stripebackgroun white -exportselection true\
               
   grid $fro2.tb -row 0 -column 0 -sticky news 
   
   ##Scroll_BAr V
   scrollbar $fro2.scr1 -orient vertical -command [list $fro2.tb  yview]
    grid $fro2.scr1 -row 0 -column 1  -sticky ens

    $fro2.tb insert end [list "101" "20.000" "12" "mayority"]
    $fro2.tb insert end [list "102" "7.000" "5" ""]
    $fro2.tb insert end [list "103" "3.000" "1" ""]

}

proc ::NeuronVND::loadObject {f} {
    variable objList
    variable historyCalls
    # current shortcut
    set objid [mol new]
    graphics $objid color white
    mol addfile $f
    mol rename $objid $f
    lappend objList [list $objid $f]
    catch {.neuron.fp.navigation.main.obj.cb configure -values $::NeuronVND::objList}
    lappend ::NeuronVND::historyCalls "::NeuronVND::loadObject $f"
}

proc ::NeuronVND::loadFiles {cfgfile {createrep true} {loadedges false}} {
    variable listmodels
    variable indexmodel
    variable historyCalls
    
    set neuronrep .neuron.fp.systems.rep;#.neuron.rep.

    ####### For the future application 
    # read files
    set success 0
    # if succesful increase indexmodel
    if {$success} {}
    # populate main table required values
    set listmodels($indexmodel,id) 0
    set listmodels($indexmodel,name) ""
    set listmodels($indexmodel,neurons) ""
    ############################################
    
    ::neuro::cmd_load_model_config_file [pwd] $cfgfile $loadedges

    # preliminary naming for the models coming from file .h5
    set listmodels(0,name) [lindex [split [lindex [::neuro::cmd_query filesets] 0 0] /] end]
    set listmodels(0,neurons) [::neuro::cmd_query num_neurons]

    .neuron.fp.systems.main.tb insert end [list "0" "T" "D" $listmodels(0,name) $listmodels(0,neurons)]

    ::NeuronVND::populateTree

    lappend ::NeuronVND::historyCalls "::NeuronVND::loadFiles $cfgfile false $loadedges"

    # Checking if default rep works
    if {$createrep} {::NeuronVND::createRepArgs}

}

# unified createRep procs
proc ::NeuronVND::createRepArgs {args} {
    variable repselected
    variable styleRep
    variable colorRep
    variable colorID
    variable selRep
    variable materialRep
    variable showRep
    variable sphereScale
    variable sphereRes
    
    # if given args
    if {[llength $args]} {
        puts "args:$args"

        # new args language
        set auxpos [lsearch $args selection]
        if {$auxpos != -1} {
            set selRep [lindex $args [incr $auxpos]]
        }

        set auxpos [lsearch $args style]
        if {$auxpos != -1} {
            set styleRep [lindex $args [incr $auxpos]]
        }

        set auxpos [lsearch $args material]
        if {$auxpos != -1} {
            set materialRep [lindex $args [incr $auxpos]]
        }

        set auxpos [lsearch $args show]
        if {$auxpos != -1} {
            set showRep [lindex $args [incr $auxpos]]
        }
        
        set auxpos [lsearch $args color]
        if {$auxpos != -1} {
            set content [lindex $args [incr $auxpos]]
            if {$content == "Type"} {
                 set colorRep "Type"
                 set colorID "Type"
            } else {
                 set colorRep "Color"
                 set colorID $content
            }
        }

    } elseif {$repselected == "" || $::neuro::nrepList == ""} {    
        # no rep selected, create a default one
        set styleRep soma 
        set colorID Type
        set selRep "all"
        set sphereScale 3
        set sphereRes 5
        set showRep true
        set materialRep Opaque
        set repselected 0

        # limit crowding in the default preview
        if {[::neuro::cmd_query "num_neurons"] > 10000} {set selRep "stride 5"}
        if {[::neuro::cmd_query "num_neurons"] > 100000} {set selRep "stride 50"}
        if {[::neuro::cmd_query "num_neurons"] > 1000000} {set selRep "stride 500"}
    } else {
        #instead of defining the rep, use the selected to get a copy
        set auxrow [.neuron.fp.systems.rep.main.table.tb get [.neuron.fp.systems.rep.main.table.tb curselection]]
        puts "auxrow: $auxrow"
        if {[llength $auxrow]} {
            lassign $auxrow styleRep colorID selRep
            set repselected [.neuron.fp.systems.rep.main.table.tb index end]
        }
    }
    
    # main call
    set repselected [.neuron.fp.systems.rep.main.table.tb index end]
    puts "repselected $repselected, styleRep $styleRep, colorID $colorID, selRep $selRep"
    set repid [::neuro::cmd_create_rep_node_fullsel $styleRep $colorID $materialRep $selRep]
    
    # insert repid details in table
    set rowid [.neuron.fp.systems.rep.main.table.tb insert $repselected [list $styleRep $colorID $selRep]]
    # set table curselection and repselected
    .neuron.fp.systems.rep.main.table.tb selection clear 0 end
    .neuron.fp.systems.rep.main.table.tb selection set $rowid
    set repselected [.neuron.fp.systems.rep.main.table.tb curselection]
    # update GUI elements for rep
    ::NeuronVND::updateRepMenu

    # hide rep if status not true
    if {!$showRep} {
       ::neuro::cmd_hide_rep $repid
       .neuron.fp.systems.rep.main.table.tb rowconfigure $rowid -foreground red
       .neuron.fp.systems.rep.main.table.tb rowconfigure $rowid -selectforeground red
    }

}

proc ::NeuronVND::delRep {args} {
    variable repselected

    #set repselected [.neuron.fp.systems.rep.table.tb curselection]
    set repid [lindex $::neuro::nrepList $repselected 0]
    ::neuro::cmd_delete_rep $repid
    .neuron.fp.systems.rep.main.table.tb delete $repselected 
    puts "delRep: selRep=$::NeuronVND::selRep, styleRep=$::NeuronVND::styleRep, colorID=$::NeuronVND::colorID" 

    #after deletion, select repselected - 1 row in table
    .neuron.fp.systems.rep.main.table.tb selection clear 0 end
    .neuron.fp.systems.rep.main.table.tb selection set [expr $repselected - 1]
    incr repselected -1

}

proc ::NeuronVND::showHideRep {} {
    variable repselected

    #set repselected [.neuron.fp.systems.rep.table.tb curselection]
    set status [.neuron.fp.systems.rep.main.table.tb rowcget $repselected -foreground]
    set repid [lindex $::neuro::nrepList $repselected 0]
    if {$status == "red"} {
        # rep was hidden, show it
        ::neuro::cmd_show_rep $repid
        # update foreground color
        .neuron.fp.systems.rep.main.table.tb rowconfigure $repselected -foreground black
        .neuron.fp.systems.rep.main.table.tb rowconfigure $repselected -selectforeground black
    } else {
        # rep was showing, hide it
        ::neuro::cmd_hide_rep $repid
        # update foreground color
        .neuron.fp.systems.rep.main.table.tb rowconfigure $repselected -foreground red
        .neuron.fp.systems.rep.main.table.tb rowconfigure $repselected -selectforeground red
    }
}

proc ::NeuronVND::updateRepMenu {} {
    variable repselected
    variable styleRep
    variable materialRep
    variable colorRep
    variable colorID
    variable sphereScale
    variable sphereRes
    variable selRep

    # get rep details from neuro_read
    set repdetails [lindex $::neuro::nrepList $repselected]
    puts "updateRepMenu: repdetails = $repdetails"
    # make rep top molecule
    mol top [lindex $repdetails 2]
    # idea: when a rep is selected in the table, populate the selection, style and color entry/boxs
    set styleRep [lindex $repdetails 3]
    if {[lindex $repdetails 4] != "Type"} {
        set colorRep "Color"
        set colorID [lindex $repdetails 4]
        grid .neuron.fp.systems.rep.main.def.coloridcb -row 1 -column 1 -sticky news
    } else {
        set colorRep "Type"
        grid remove .neuron.fp.systems.rep.main.def.coloridcb
    }
    set materialRep [lindex $repdetails 5]
    set selRep [lindex $repdetails 6]
    if {$styleRep == "soma"} {
        set sphereScale [lindex $repdetails 10]
        set sphereRes [lindex $repdetails 11]
    }
    .neuron.fp.systems.rep.main.table.tb selection clear 0 end
    .neuron.fp.systems.rep.main.table.tb selection set $repselected
    puts "updateRepMenu: repselected = $repselected, curselection = [.neuron.fp.systems.rep.main.table.tb curselection]"
    #.neuron.fp.systems.rep.main.sel.entry delete 0 end
    #.neuron.fp.systems.rep.main.sel.entry insert 0 $selRep

    # update variable for connectivity purposes
    variable listOfRepsForConnect
    set listOfRepsForConnect ""
    foreach r [.neuron.fp.systems.rep.main.table.tb get 0 end] {
        lappend listOfRepsForConnect [lindex $r 2]
    }
    
}

proc ::NeuronVND::editRep {case} {
    variable repselected
    variable styleRep
    variable materialRep
    variable colorRep
    variable colorID
    variable selRep
    variable sphereScale
    variable sphereRes

    set t .neuron.fp.systems.rep.main.table.tb

    # check the selected rep has a different style, get rep details from neuro_read
    set repdetails [lindex $::neuro::nrepList $repselected]
    if {$repdetails == ""} {return}

    # define color to be added to table
    if {$colorRep == "Type"} {set color Type}
    if {$colorRep == "Color"} {set color $colorID}

    switch $case {
        "style" {
            puts "editRep: styleRep = $styleRep"
            if {$styleRep != [lindex $repdetails 3]} {
                # delete previous rep
                ::NeuronVND::delRep
                # create a copy rep with a different style
                ::NeuronVND::createRepArgs style $styleRep               
            }
        }
        "sel" {
            if {$selRep != [lindex $repdetails 6]} {
                # delete previous rep
                ::NeuronVND::delRep
                # create a copy rep with a different sel
                ::NeuronVND::createRepArgs selection $selRep
            }
        }
        "color" {
            if {$colorID != [lindex $repdetails 4]} {
                # delete previous rep
                ::NeuronVND::delRep
                # create a copy rep with a different color
                ::NeuronVND::createRepArgs color $color

            }
        }
        "material" {
            if {$materialRep != [lindex $repdetails 5]} {
                # update neuro::nrepList
                set ::neuro::nrepList [lreplace $::neuro::nrepList $repselected $repselected [lreplace [lindex $::neuro::nrepList $repselected] 5 5 $materialRep]]
                # changing material through draw command
                # needs to make mol top then call draw material xxxx
                mol top [lindex $repdetails 2]
                draw material $materialRep
            }
        }
        # SPHERE IS NOT WORKING AT THE MOMENT
        "sphere" {
            if {$sphereScale != [lindex $repdetails 10] || $sphereRes != [lindex $repdetails 11]} {
                # create a copy rep with a different style
                ::NeuronVND::createRepArgs
                # delete previous rep
                ::NeuronVND::delRep
            }
        }
    }
}

proc ::NeuronVND::moveGraphs {dim sign} {
  variable objList
  variable movex
  variable movey
  variable movez
  variable aggoffset
  set objid [lindex $objList 0 0]
  if {$objid == ""} {return}
  set numG [llength [graphics $objid list]]
  switch $dim {
    "x" { 
        set val $movex
        if {$sign == "neg"} {set val [expr -1*$val]}
        set offset [list $val 0.0 0.0] 
    }
    "y" { 
        set val $movey
        if {$sign == "neg"} {set val [expr -1*$val]}
        set offset [list 0.0 $val 0.0] 
    }
    "z" { 
        set val $movez
        if {$sign == "neg"} {set val [expr -1*$val]}
        set offset [list 0.0 0.0 $val] 
    }
    "default" {puts "error: dimension must be either x, y or z"}
  }
  # update aggregate offset to save state
  set aggoffset [vecadd $aggoffset $offset]

  display update off
  for {set i 1} {$i < $numG} {incr i} {
    lassign [graphics $objid info $i] t v1 v2 v3
    # offset v1 v2 v3
    set newv1 [vecadd $v1 $offset]
    set newv2 [vecadd $v2 $offset]
    set newv3 [vecadd $v3 $offset]
    # redraw graphics i
    graphics $objid replace $i
    graphics $objid triangle $newv1 $newv2 $newv3
  }
  display update on
}

proc ::NeuronVND::rotGraphs {axis sign} {
  variable objList
  variable rotarx
  variable rotary
  variable rotarz
  variable aggrot
  set objid [lindex $objList 0 0]
  if {$objid == ""} {return}
  set numG [llength [graphics $objid list]]
  if {$axis != "x" && $axis != "y" && $axis != "z"} {
    puts "error: axis must be either x, y or z"
    return
  }
  switch $axis {
    "x" { set val $rotarx }
    "y" { set val $rotary }
    "z" { set val $rotarz }
  }
  if {$sign == "neg"} {set val [expr -1*$val]}

  # update aggrot 
  set aggrot [transmult $aggrot [transaxis $axis $val]]

  display update off
  for {set i 1} {$i < $numG} {incr i} {
    lassign [graphics top info $i] t v1 v2 v3
    # offset v1 v2 v3
    set newv1 [vectrans [transaxis $axis $val] $v1]
    set newv2 [vectrans [transaxis $axis $val] $v2]
    set newv3 [vectrans [transaxis $axis $val] $v3]
    # redraw graphics i
    graphics top replace $i
    graphics top triangle $newv1 $newv2 $newv3
  }
  display update on
}


proc ::NeuronVND::revealVars {repdetails} {

    set show [lindex $repdetails 1] 
    set style [lindex $repdetails 3]
    set color [lindex $repdetails 4]
    set material [lindex $repdetails 5]
    set selection [lindex $repdetails 6]
    set scale [lindex $repdetails 10] 
    set resolution [lindex $repdetails 11]

    set result "show $show style $style color $color material $material selection $selection";#scale $scale resolution $resolution"

    return $result
}

proc ::NeuronVND::visState {mode} {
    variable historyCalls

    set types {
	 {{TCL files} {.tcl}   }
	 {{All Files}        *            }
    }

   # Check ns is a correct namespace

   #################
   
   switch $mode {
      "save" {
         set newpathfile [tk_getSaveFile \
			  -title "Choose file name" \
			  -initialdir [pwd] -filetypes $types]
         if {$newpathfile == ""} {return}
         set fid [open $newpathfile w]
         puts $fid "# Visual Neuronal Dynamics"
         puts $fid "# Visualization State"
         foreach call $historyCalls {
               puts $fid $call
         }
         puts $fid "# List of representations"
         foreach r $neuro::nrepList {
               set s [::NeuronVND::revealVars $r]
               puts $fid "::NeuronVND::createRepArgs $s" 
         }
         
         close $fid  
      }
      "load" {
         set newpathfile [tk_getOpenFile \
			  -title "Choose file name" \
			  -initialdir [pwd] -filetypes $types]
         if {$newpathfile == ""} {return}
         set fid [open $newpathfile r]
         #check this is a multiplot options file
         set line [gets $fid]
         if {[regexp {# Visual Neuronal Dynamics} $line] == 1} {
               source $newpathfile   
         } else {puts "VND) Not a valid state file"}
      }
   } 

}

proc ::NeuronVND::populateTree {} {

    set tv .neuron.fp.systems.info.tv
    # Populate model tree with population
    # file level
    foreach f [::neuro::cmd_query fileset_pop_groups] {
        # pop level
        foreach p [lindex $f 1] {
            set popname [lindex $p 0]
            $tv insert {} end -id $popname -text "population == $popname"
            # group level
            foreach g [lindex $p 1] {
                $tv insert $popname end -id ${popname}_$g -text "group == $g"
                # type level
                foreach t [::neuro::cmd_query node_types_in_group [lindex $f 0] $popname $g] {
                    $tv insert ${popname}_$g end -id ${popname}_${g}_$t -text "type == $t"
                }
            }
        }

    }
}

proc ::NeuronVND::createExampleRep {args} {
    # IDEA: create a temporary rep with selection defined by treeview
    # split str and check length
    set str [split [lindex $args end] _]
    switch [llength $str] {
        "1" {
            set sel "population == [lindex $str 0]"
            puts "# selected: $sel"
        }
        "2" {
            set sel "population == [lindex $str 0] && group == [lindex $str 1]"
            puts "# selected: $sel"
        }
        "3" {
            set sel "population == [lindex $str 0] && group == [lindex $str 1] && type == [lindex $str 2]"
            puts "# selected: $sel"
        }
    }

    ::NeuronVND::createRepArgs style soma selection $sel color yellow

}

proc ::NeuronVND::exit {} {
    # prompt exit confirmation
    set answer [tk_messageBox -message "Do you want to quit VND?" -type yesno -title "Closing VND" -icon info -parent $::NeuronVND::topGui]
    if {$answer == "no"} {
        return
    }
    destroy .neuron

    # Trying to quit the whole VND/VMD program from here
    # but when using "exit" or "quit" it prompts 'bad window path name' error on .neuron
    #exit

}
