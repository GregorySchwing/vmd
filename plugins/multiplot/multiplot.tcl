#
# 2-D plotting tool
#
# $Id: multiplot.tcl,v 1.44 2021/11/19 17:08:33 mariano Exp $
#
# Author:
# Jan Saam
# Beckman Institute
# University of Illinois
# saam@ks.uiuc.edu

# You can have several independent instances of "multiplot" runnning at the 
# same time because the data are kept in different namespaces.
# It returns a plothandle which you can use to control an existing plot,
# add new datasets or quit it.

# Usage:
# set plothandle [multiplot ?reset|list|embed <path>? ?options?]
# reset --- Closes all windows and deletes all namespaces and plothandles
# list  --- Lists all existing plothandles
# embed --- Provide path to a parent widget, so that multiplot can be 
#           embedded instead of being a toplevel window.

# You can use the returned plothandle to control the plot:
# $plothandle add|replot|namespace|configure|data|export|quit ?options?

# $plothandle add X Y ?options?   --- Adds a dataset to the plot
# $plothandle replot              --- Replots the current data
# $plothandle namespace           --- Returns the current namespace
# $plothandle configure ?options? --- Modifies the existing plot according to the options
#                                     These modifications are silent until you call 'replot'
#                                     unless you specify the flag -plot
# $plothandle nsets               --- Returns the number of datasets in the plot
# $plothandle data                --- Returns all x and y datasets of the plot
# $plothandle xdata               --- Returns all x datasets of the plot
# $plothandle ydata               --- Returns all y datasets of the plot
# $plothandle getpath             --- Returns the widget path to the plot window
# $plothandle export program filename  --- Export plot to external program
# $plothandle draw <item> ?options?    --- Draws arbitrary item into canvas. 
#                                          The item can be one of the items defined in Tk's
#                                          'canvas' command (arc, bitmap, image, line, oval,
#                                          polygon, rectangle, text, window). The options are
#                                          the same as for the corresponding 'canvas create <item>'
#                                          command (they are directly passed to canvas create).
#                                          Don't dare to ask me about these options, read the Tk manual!
#                                          For your convenience I have added wrappers for the oval
#                                          and rectangle items where you can specify the center 
#                                          coordinates instead of the two opposite corners of the
#                                          bounding box. These two additional items are called
#                                          'circle' and 'square'. The size is controlled by the
#                                          additional flag -radius (which in case of square denotes
#                                          the half side length).
# $plothandle undraw <tag>        --- Deletes an item that is associated with a tag.
#                                     The item can be one of the items defined in Tk's 'canvas'
#                                     command (see drawing-command above).
# $plothandle clear               --- Removes all datasets but does not delete the plothandle
# $plothandle quit                --- Destroy the window of the plot and delete all data

# Options for the plothandle:
# ===========================
# Switches:
# -lines       --- Connect datapoint with lines
# -nolines     --- Don't connect datapoint with lines
# -stats       --- Print some statistics of the last added dataset
# -nostats     --- Get rid of the statistics
# -plot        --- Actually plot the data otherwise only the canvas and the axes are drawn
#                  This is equivalent to '$plothandle replot'.
#                  If you have multiple datasets it is a good idea to add all data first and 
#                  then plot them all at once since this will be a lot faster.
# -autoscale   --- Automatically scale plot to fit all data points
# -xanglescale --- Use 90 degree as major tic unit for the x-axis
# -yanglescale --- Use 90 degree as major tic unit for the y-axis

# Options with argument:
# -set i           --- Apply all dataset specific modifications to dataset i
# -x X             --- Supply the x-values for a dataset in a list
# -y Y             --- Supply the y-values for a dataset in a list
# -title text      --- Title of the plot
# -xlabel text     --- Text for the x-axis label
# -ylabel text     --- Text for the y-axis label
# -xmajortics i    --- Distance between two x-axis ticlabels
# -ymajortics i    --- Distance between two y-axis ticlabels
# -xminortics i    --- Distance between two x-axis minor tic marks
# -yminortics i    --- Distance between two y-axis minor tic marks
# -xsize s         --- Width of the canvas
# -ysize s         --- Height of the canvas
# -xmin s          --- Minimum x value; use "auto" to take the minimum x value of all datasets
# -xmax s          --- Maximum x value; use "auto" to take the maximum x value of all datasets
# -ymin s          --- Minimum y value; use "auto" to take the minimum y value of all datasets
# -ymax s          --- Maximum y value; use "auto" to take the maximum y value of all datasets
# -marker type     --- Draw markers at datapoints (none|point|circle|square)
# -radius r        --- Data point marker (radius of circle and point, size of square)
# -fillcolor color --- Fill color of datapoint markers (option can be abbreviated with -fill)
# -linewidth w     --- Width of the lines connecting datapoints
# -linecolor color --- Color of the lines connecting datapoints
# -bkgcolor color  --- Color of the canvas background
# -dash pattern    --- Draw dashed lines. The dash pattern is specified by one of the
#                      following characters "-,._" (uses the same format as -dash for Tk canvas)
#                      Note that each line segment is dashed. Hence you'll get a solid line when
#                      the datapoints are so dense that the line segments are shorter than the dashes!
# -legend text     --- Add an entry for this dataset to the legend
#                      Note that the legend is drawn in the upper left corner of the plot
#                      but you can drag the legend anywhere you want using the mouse.
# -hline {y args}  --- Draw a horizontal line at position y, args are arguments for the Tk canvas 
#                      'create line' command. Through args like '-width 2 -fill red -dash "-"' you can 
#                      determine the line style.
# -vline {x args}  --- Draw a vertical line at position x, args are arguments for the Tk canvas 
#                      'create line' command. Through args like '-width 2 -fill red -dash "-"' you can 
#                      determine the line style.
# XXX initial support for NMWiz callbacks
# -callback proc   --- Bind a callback procedure to mouse click actions on markers. 
#                      Five arguments, {index x y color marker}, will be passed to the procedure.
#                      Callback procedures will be called only for datasets associated with the procedure,
#                      i.e. when "-callback proc" is specified during addition of a dataset.
#                      A usage example is provided below.

# Examples:
# ---------
#   package require multiplot
#   set x {-2 -1 0 1 2 3 4 5 6 7 8 9 10}
#   set y {-2  0 2 3 4 5 5 4 3 2 1 0 1}
#   proc click {args} {puts "[lindex $args 4] at x=[lindex $args 1] and y=[lindex $args 2] is colored [lindex $args 3]"}

### This plot will be immediately created because we specified -plot 
#   set plothandle [multiplot -x $x -y $y -title "Example plot" -lines -linewidth 3 -marker point -callback click -plot]

### Now we change the appearence of the existing plot.
### BUT WE WON'T SEE THIS change until the next replot is requested!
#   $plothandle configure -fillcolor yellow -radius 6 -plot

### Let's add a vertical dotted line at x=3
#   $plothandle configure -vline {3 -width 2 -fill red -dash "."} -plot

### And now redraw the plot so that the changes become visible:
#   $plothandle replot;
  
### It's time to add a second dataset to the same plot
#   set y2 {6 7 8 7 6 6 5 4 4 3 2 3 4}
#   $plothandle add $x $y2 -marker circle  -fillcolor green -radius 4 -callback click -plot

### Of course we can change the appearence of the the two sets independently:
#   $plothandle configure -set 1 -lines -linewidth 4 -dash "," -plot

### Export to xmgrace, load with 'xmgrace -nxy /tmp/foo.plot'
#   $plothandle export xmgrace /tmp/foo.plot

### Close the plot
#   $plothandle quit

package require exectool  1.2
package provide multiplot 1.7

namespace eval ::MultiPlot:: {
   proc initialize {} {
      variable plotlist {}
      variable plotcount -1
      variable parent
      variable verbose 0
      variable ptable ""
      variable plugin ""
      variable qmspan 5
      ############ new variables for editplot
      variable editplotysize ""
      variable editplotxsize ""
      variable editlist ""
      variable editselected "" ;# the plot selected for edition
      variable editlastselected ""; # additional workaround for spinbox overwriting
      variable editdatasetlist ""
      variable editdataset ""  ;# the dataset selected for edition
      variable editauxvalues "" ;# workaround for spinbox overwriting given values
      variable edittitle ""
      variable edittitlefont "Helvetica"
      variable edittitlefontsize "10"
      variable editylabel ""
      variable editylabelfont "Helvetica"
      variable editylabelfontsize "10"
      variable editxlabel ""
      variable editxlabelfont "Helvetica"
      variable editxlabelfontsize "10"
      variable editymin ""
      variable editymax ""
      variable editxmin ""
      variable editxmax ""
      variable editlinedash ""
      variable editlinewidth ""
      variable editlinecolor ""
      variable editmarkersize ""
      variable editmarkercolor ""
      variable editmarkertype ""
   }
   initialize
}

proc ::MultiPlot::init_plot {args} {
   variable parent
   variable verbose

   set  parent [lindex $args 0]
   incr ::MultiPlot::plotcount
   set ns "::MultiPlot::Plot${::MultiPlot::plotcount}"

   if {$verbose} {
     if {[namespace exists $ns]} {
       puts "Reinitializing namespace $ns."
     } else {
       puts "Creating namespace $ns"
     }
   }

   namespace eval $ns {
      # Default values
      variable nsets 0
      variable datasets
      array unset datasets
      lappend args -set 0
      variable curset 0
      variable title {}
      variable titlefontsize 10
      variable ticfontsize   8
      variable ylabelfontsize 10
      variable xlabelfontsize 10
      variable titlefont  "Helvetica"
      variable ticfont    "Helvetica $ticfontsize"
      variable ylabelfont  "Helvetica"
      variable xlabelfont  "Helvetica"
      variable infoFont   {Courier 9}
      variable postscript "multiplot.ps"
      variable printstats 0
      variable replot 0
      variable canh   700;   # canvas height
      variable canw   1000;  # canvas width
      variable resize 0
      variable ticlen 10; # length of the major tic marks
      variable rim    8;  # extra space around the plot
      variable xlabeloffset
      variable xlabeltext ""
      variable ylabeloffset
      variable ylabeltext ""
      variable ylabelorient ""
      variable lines     1;        # connect data points with lines?
      variable marker    "none";   # display data points [circle|point|square|none]
      variable radius    2;        # radius of circles and points , size of squares
      variable linewidth 1;        # width of the line connecting data points
      variable fillcolor Skyblue2; # fill color of data point markers   
      variable linecolor black;    # color of lines connecting data points
      variable bkgcolor white;    # color of the canvas background
      variable dashed    {{}};     # Draw dashed lines (uses the same format as -dash for Tk canvas)
      variable legend    {{}};     # legend string for current dataset
      variable colorlist {black red green blue magenta orange OliveDrab2 cyan maroon gold2 yellow gray60 SkyBlue2 orchid3 ForestGreen PeachPuff LightSlateBlue}
                        # XXX initial support for NMWiz callbacks
      variable callback "none";        # Click callback procedure 

      variable predefRange 0
      variable givenXmin auto
      variable givenYmin auto
      variable givenXmax auto
      variable givenYmax auto

      variable xmin   0
      variable xmin_y 0
      variable ymin   0
      variable ymin_x 0
      variable xmax   0
      variable xmax_y 0
      variable ymax   0
      variable ymax_x 0
      variable spanx  0
      variable spany  0
      variable anglescalex 0
      variable anglescaley 0
      variable xmajortics {}
      variable ymajortics {}
      variable xminortics {}
      variable yminortics {}

      variable hline    {}
      variable vline    {}
      variable xplotmin {}
      variable yplotmin {}
      variable xplotmax {}
      variable yplotmax {}
      variable scalex {}
      variable scalex {}
      variable dimx 0
      variable dimy 0
      variable minorticx 5
      variable minorticy 5

      variable objectlist {}; # other drawn objects like circles, text, lines,...
      variable redraw 0;      # redraw all objects

      variable w ${::MultiPlot::parent}.plotwindow${::MultiPlot::plotcount}
      variable istoplevel 1;  # set to 0 for embedded widgets
      variable c
      variable namespace ::Plothandle${::MultiPlot::plotcount}

      variable gridy "none" ;# to define grid in the vertical y axis: none, major, all
      variable gridx "none" ;# to define grid in the horizontal x axis: none, major, all

      if {${::MultiPlot::parent} != ""} {
        set istoplevel 0
        set canh 350
        set canw 500
      } else {
        set istoplevel 1
      } 

      catch {destroy $w}
      if {$istoplevel} {
        toplevel $w
        wm title $w "MultiPlot"
        wm iconname $w "MultiPlot"
        wm protocol $w WM_DELETE_WINDOW "[namespace current]::plothandle quit"
        wm withdraw $w
      } else {
        frame $w -bd 0
        grid $w 
      }

      frame $w.menubar -relief raised -bd 2
      menubutton $w.menubar.file -text "File" -underline 0 \
         -menu $w.menubar.file.menu
      $w.menubar.file config -width 3 

      menu $w.menubar.file.menu -tearoff 0
      
      $w.menubar.file.menu add command -label "Export to PostScript" -command "[namespace current]::savedialog "
      $w.menubar.file.menu add command -label "Export to Xmgrace" -command "[namespace current]::xmgracedialog "
      $w.menubar.file.menu add command -label "Export to ASCII matrix..." -command "[namespace current]::savematrix "
      $w.menubar.file.menu add command -label "Export to ASCII vectors..." -command "[namespace current]::savevectors "
      $w.menubar.file.menu add command -label "Print plothandle in Console" -command "vmdcon [namespace current]::plothandle"
      if {$istoplevel} {
        $w.menubar.file.menu add command -label "Quit" -command "[namespace current]::plothandle quit"
      }

      ## TESTING NEW EDIT GUI
      menubutton $w.menubar.edit -text "Edit" -underline 0 \
         -menu $w.menubar.edit.menu
      $w.menubar.edit config -width 4

      menu $w.menubar.edit.menu -tearoff 0
      
      $w.menubar.edit.menu add command -label "Edit Plot" -command "[namespace current]::editplot"
      $w.menubar.edit.menu add command -label "Load Plot Style" -command "[namespace current]::managestyle load"
      $w.menubar.edit.menu add command -label "Save Plot Style" -command "[namespace current]::managestyle save"
      grid $w.menubar.file $w.menubar.edit
      ########################################

      grid $w.menubar -sticky news


      if {![winfo exists $w.f.cf]} {
         variable canw
         variable canh
         variable c $w.f.cf
         frame $w.f 
         grid $w.f -sticky news

         canvas $c -relief flat -borderwidth 0 -width $canw -height $canh -bg $bkgcolor 
         scrollbar $w.f.y -orient vertical   -command [namespace code {$c yview}]
         scrollbar $w.f.x -orient horizontal -command [namespace code {$c xview}]
         $c configure  -yscrollcommand [namespace code {$w.f.y set}] -xscrollcommand [namespace code {$w.f.x set}]
         $c configure  -scrollregion   "0 0 $canw $canh"
         grid $c $w.f.y 
         grid $w.f.x    
         grid rowconfigure    $w.f 0 -weight 1
         grid columnconfigure $w.f 0 -weight 1
         grid configure $w.f.y  -sticky ns
         grid configure $w.f.x  -sticky we
      }

      # Create a plothandle procedure that provides some commands to control the plot.
      # It's full name will be returned when you invoke multiplot.
      proc plothandle { command args } {
         variable w
         switch $command {
            namespace { return [namespace current] }
            replot    { variable replot 1; plot_update; return }
            add       { 
               set newX [lindex $args 0]
               set newY [lindex $args 1]

               set lenX [llength $newX]
               set lenY [llength $newY]
               if {!$lenX} { error "X data vector is empty!" }
               if {!$lenY} { error "Y data vector is empty!" }

               # Check, if we have several coordinate sets:
               if {[llength [join $newX]]>$lenX || [llength [join $newY]]>$lenY} {
                  if {$lenX != $lenY} {
                     error "Different number of datasets for x and y ($lenX!=$lenY)"
                  }
                  foreach x $newX y $newY {
                     eval add_data [list $x] [list $y] [lrange $args 2 end]
                  }
               } else {
                  eval add_data [list $newX] [list $newY] [lrange $args 2 end]
               }
               plot_update 
            }
	         clear {
               variable datasets
               variable nsets 0
	            array unset datasets
	         }
            draw {
               # Register the new object
               variable objectlist
               lappend objectlist $args

               # Make sure that the plot geometry was calculated already and draw
               variable xplotmin
               if {![llength $xplotmin]} {
                  variable redraw 1;
                  variable resize 1;
                  plot_update; # implicitely draws all objects
               } else {
                  draw_object $args
               }
            }
            undraw {
                 undraw_object $args
            }
            configure { 
               variable datasets
               variable nsets
               variable resize

               variable curset {} 
               set pos [lsearch $args "-set"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  variable curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
               }
               if {![llength $curset]} { set curset 0 }

               set havedata 0
               set pos [lsearch $args "-x"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  if {$nsets==0} {
                     lappend datasets(X) {}
                     lappend datasets(Y) {}
                     lappend datasets(xmin) {}
                     lappend datasets(xmax) {}
                     lappend datasets(xmin_y) {}
                     lappend datasets(xmax_y) {}
                     lappend datasets(ymin) {}
                     lappend datasets(ymax) {}
                     lappend datasets(ymin_x) {}
                     lappend datasets(ymax_x) {}
                     incr nsets
                  }
                  lset datasets(X) $curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
                  variable resize 1
                  set havedata 1
               }

               set pos [lsearch $args "-y"]
               if {$pos>=0 && $pos+1<[llength $args]} { 
                  if {$nsets==0} {
                     lappend datasets(X) {}
                     lappend datasets(Y) {}
                     lappend datasets(xmin) {}
                     lappend datasets(xmax) {}
                     lappend datasets(xmin_y) {}
                     lappend datasets(xmax_y) {}
                     lappend datasets(ymin) {}
                     lappend datasets(ymax) {}
                     lappend datasets(ymin_x) {}
                     lappend datasets(ymax_x) {}
                     incr nsets
                  }
                  lset datasets(Y) $curset [lindex $args [expr $pos+1]]
                  set args [lreplace $args $pos [expr $pos+1]]
                  variable resize 1
                  set havedata 1
               }

               plot_scan_options $args; 

               if {$resize && $havedata} {
                  init_dataset
               }

               plot_update 
            }
            nsets      {
               variable datasets;
               if { [info exists datasets(Y)] } {
                  return [llength $datasets(Y)]
               } else {
                  return 0
               }
            }
            all      {
               variable datasets;
               return [array get datasets]
            }
            xdata      {
               variable datasets;
               return $datasets(X)
            }
            ydata      {
               variable datasets;
               return $datasets(Y)
            }
            data      { 
               variable datasets;
               return [list $datasets(X) $datasets(Y)]
            }
            getpath   {
               variable w;
               return $w
            }
            export    {
              variable datasets
              variable title
              variable legend
              variable nsets
              if { [llength $args] < 2} {
                vmdcon -err "Incorrect export syntax"
                return
              }
              set progname [lindex $args 0]
              set filename [lindex $args 1] 

              switch $progname {
                grace - 
                xmgr  -
                xmgrace {
                  vmdcon -info "Exporting plot in xmgrace format as filename $filename"
                  set fd [open $filename "w"]
                  puts $fd "@type xy"
                  puts $fd "@title \"$title\""
                  set ylen [llength $datasets(Y)]
                  for {set s 0} {$s < $nsets} {incr s} {
                    if {[lindex $legend $s] != ""} {
                      puts $fd "@s$s legend \"[lindex $legend $s]\""
                    }
                  }
                  for {set s 0} {$s < $nsets} {incr s} {
                    set xlen [llength [lindex $datasets(X) $s]]
                    for {set i 0} {$i < $xlen} {incr i} {
                      puts $fd "[lindex $datasets(X) $s $i] [lindex $datasets(Y) $s $i]"
                    }
                    puts $fd "&"
                  }                    
                  close $fd
                }
 	            vectors {	# TONI
                  vmdcon -info "Exporting plot as ASCII vectors in file $filename"
                  set fd [open $filename "w"]
                  for {set s 0} {$s < $nsets} {incr s} {
                     set xlen [llength [lindex $datasets(X) $s]]
                     for {set i 0} {$i < $xlen} {incr i} {
                        puts $fd "[lindex $datasets(X) $s $i] [lindex $datasets(Y) $s $i]"
                     }
                     puts $fd ""
                  }                    
                  close $fd
               }
               matrix {	# TONI
                  if {![cansavematrix]} {
                     vmdcon -err "Not in matrix form"
                     return
                  }
                  vmdcon -info "Exporting plot as ASCII matrix in file $filename"
                  set fd [open $filename "w"]
                  set xlen [llength [lindex $datasets(X) 0]]
                  for {set i 0} {$i < $xlen} {incr i} {
                     puts -nonewline $fd "[lindex $datasets(X) 0 $i] "; # X value
                     for {set s 0} {$s < $nsets} {incr s} {
                     puts -nonewline $fd "[lindex $datasets(Y) $s $i] "
                     }
                     puts $fd ""
                  }                    
                  close $fd
               }
              }

               return ;# export option
            }
            quit   { 
               destroy $w;
               namespace delete [namespace current]
               return
            }
            # Testing adding more commands to retrieve info about plothandle
            title { return [set [namespace current]::title] }
            titlefont { return [set [namespace current]::titlefont] }
            titlefontsize { return [set [namespace current]::titlefontsize] }
   
         }
      }

      proc init_dataset {} {
         variable datasets
         variable curset
         set minx [lindex $datasets(X) $curset 0]
         set maxx [lindex $datasets(X) $curset end]
         set miny [lindex $datasets(Y) $curset 0]
         set maxy [lindex $datasets(Y) $curset end]
         set minx_y [lindex $datasets(Y) $curset 0]
         set maxx_y [lindex $datasets(Y) $curset end]
         set miny_x [lindex $datasets(X) $curset 0]
         set maxy_x [lindex $datasets(X) $curset end]
         foreach x [lindex $datasets(X) $curset] y [lindex $datasets(Y) $curset] {
            if {$x<$minx} {
               set minx   $x
               set minx_y $y
            }
            if {$x>$maxx} {
               set maxx   $x
               set maxx_y $y
            }
            if {$y<$miny} {
               set miny   $y
               set miny_x $x
            }
            if {$y>$maxy} {
               set maxy   $y
               set maxy_x $x
            }
         }
         lset datasets(xmin)   $curset $minx
         lset datasets(xmin_y) $curset $minx_y
         lset datasets(xmax)   $curset $maxx
         lset datasets(xmax_y) $curset $maxx_y
         lset datasets(ymin)   $curset $miny
         lset datasets(ymin_x) $curset $miny_x
         lset datasets(ymax)   $curset $maxy
         lset datasets(ymax_x) $curset $maxy_x
      }

      proc plot_scan_options { arg } {
         set drawlines {}
         set points 0
         variable printstats
         variable anglescalex
         variable anglescaley

         # Scan for single options
         set argnum 0
         set arglist $arg
         foreach i $arg {
            if {$i=="-lines"}  then {
               set drawlines 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-nolines"}  then {
               set drawlines 0
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-stats"}  then {
               set printstats 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-plot"}  then {
               variable replot 1
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-nostats"}  then {
               set printstats 0
               set arglist [lreplace $arglist $argnum $argnum]
               continue
            }
            if {$i=="-xanglescale"}  then {
               set anglescalex 1
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-yanglescale"}  then {
               set anglescaley 1
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1 
               continue
            }
            if {$i=="-autoscale"} then {
               variable predefRange 0
               variable givenXmin auto
               variable givenXmax auto
               variable givenYmin auto
               variable givenYmax auto
               set arglist [lreplace $arglist $argnum $argnum]
               variable resize 1
               continue
            }
            incr argnum
         }

         # must search for the dataset option first
         variable nsets
         variable curset 
         foreach {i j} $arglist {
            if {$i=="-set"}       then { 
               if {$j>=$nsets} {
                  error "Dataset $j doesn't exist"
               }
               variable curset $j;
            }
         }

         #variable curset
         if {[llength $drawlines]} {
            variable lines
            if {![llength $curset]} {
               for {set s 0} {$s<$nsets} {incr s} {
                  lset lines $s $drawlines
               }
            } else {
               lset lines $curset $drawlines
            }
         }

         # Scan for options with one argument
         variable hline
         variable vline
         variable datasets
         foreach {i j} $arglist {
#           if {$i=="-x"}          then { 
#              if {![llength [array get datasets X]]} {
#                 lappend datasets(X) $curset $j;
#              } else {
#              lset datasets(X) $curset $j
#              }
#              variable resize 1
#           }
#           if {$i=="-y"}          then { 
#              if {![llength [array get datasets Y]]} {
#                 lappend datasets(Y) $curset $j;
#              } else {
#                 lset datasets(Y) $curset $j
#              }
#              variable resize 1
#           }
            if {$i=="-title"}      then { variable title $j; variable resize 1 }
            if {$i=="-titlefont"}      then { variable titlefont $j; variable resize 1 }
            if {$i=="-titlefontsize"}      then { variable titlefontsize $j; variable resize 1 }
            if {$i=="-labelfont"}      then { variable ylabelfont $j; variable xlabelfont $j; variable resize 1 }
            if {$i=="-labelfontsize"}      then { variable ylabelfontsize $j; variable xlabelfontsize $j; variable resize 1 }
            if {$i=="-ylabelfont"}      then { variable ylabelfont $j; variable resize 1 }
            if {$i=="-xlabelfont"}      then { variable xlabelfont $j; variable resize 1 }
            if {$i=="-ylabelfontsize"}      then { variable ylabelfontsize $j; variable resize 1 }
            if {$i=="-xlabelfontsize"}      then { variable xlabelfontsize $j; variable resize 1 }
            if {$i=="-xlabel"}     then { variable xlabeltext $j; variable resize 1 }
            if {$i=="-ylabel"}     then { variable ylabeltext $j; variable resize 1 }
            if {$i=="-xmajortics"} then { variable xmajortics $j; variable resize 1 }
            if {$i=="-ymajortics"} then { variable ymajortics $j; variable resize 1 }
            if {$i=="-xminortics"} then { variable xminortics $j; variable resize 1 }
            if {$i=="-yminortics"} then { variable yminortics $j; variable resize 1 }
            if {$i=="-xsize"}      then { variable canw $j; variable resize 1 }
            if {$i=="-ysize"}      then { variable canh $j; variable resize 1 }

            if {$i=="-xmin"}       then { variable givenXmin $j; variable predefRange 1; variable resize 1 }
            if {$i=="-xmax"}       then { variable givenXmax $j; variable predefRange 1; variable resize 1 }
            if {$i=="-ymin"}       then { variable givenYmin $j; variable predefRange 1; variable resize 1 }
            if {$i=="-ymax"}       then { variable givenYmax $j; variable predefRange 1; variable resize 1 }

            if {$i=="-hline"}      then { lappend hline $j; variable resize 1 }
            if {$i=="-vline"}      then { lappend vline $j; variable resize 1 }
            if {$i=="-bkgcolor"}   then { variable bkgcolor $j}
            if {$i=="-radius"}     then { 
               variable radius 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset radius $s $j
                  }
               } else {
                  lset radius $curset $j
               }
            }
            if {$i=="-dash"}     then { 
               variable dashed
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset dashed $s $j
                  }
               } else {
                  lset dashed $curset $j
               }
            }
            if {[string match "-fill*" $i]} then { 
               variable fillcolor 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset fillcolor $s $j;
                  }
               } else {
                  lset fillcolor $curset $j
               }
            }
            if {$i=="-linewidth"} then { 
               variable linewidth 
               variable datasets
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset linewidth $s $j
                  }
               } else {
                  lset linewidth $curset $j
               }
            }
            if {$i=="-linecolor"} then {
               variable linecolor 
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset linecolor $s $j;
                  }
               } else {
                  lset linecolor $curset $j
               }
            }
            if {[string match "-mark*" $i]} then { 
               variable marker
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset marker $s $j
                  }
               } else {
                  lset marker $curset $j
               }
            }
            if {$i=="-legend"}      then { 
               variable legend
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset legend $s $j
                  }
               } else {
                  lset legend $curset $j
               }
            }

            # XXX initial support for NMWiz callbacks
            if {$i=="-callback"}     then {
              variable callback
               if {![llength $curset]} {
                  for {set s 0} {$s<$nsets} {incr s} {
                     lset callback $s $j
                  }
               } else {
                  lset callback $curset $j
               }
            } 
         }
      }

      proc undraw_object {args} {
          variable c
          
          $c delete $args
      }
      
      proc add_data {x y args} {
         if {[llength $x] != [llength $y]} {
            error "Different number of x and y coordinates ([llength $x]!=[llength $y])"
         }
         variable datasets
         variable nsets
         variable curset $nsets
         variable lines 
         variable linewidth 
         variable linecolor
         variable marker
         variable fillcolor 
         variable dashed
         variable radius
         variable legend
         variable colorlist
         # XXX initial support for NMWiz callbacks
         variable callback

         lappend datasets(X) $x
         lappend datasets(Y) $y
         lappend datasets(xmin)   {}
         lappend datasets(xmax)   {}
         lappend datasets(xmin_y) {}
         lappend datasets(xmax_y) {}
         lappend datasets(ymin)   {}
         lappend datasets(ymax)   {}
         lappend datasets(ymin_x) {}
         lappend datasets(ymax_x) {}
         lappend lines     1
         lappend linewidth 1
         lappend linecolor [lindex $colorlist [expr {(1+$nsets)%[llength $colorlist]}]]
         lappend marker    "none"
         lappend fillcolor [lindex $colorlist [expr {(1+$nsets)%[llength $colorlist]}]]
         lappend radius    2
         lappend dashed    {}
         lappend legend    {}
         # XXX initial support for NMWiz callbacks
         lappend callback  "none"

         # Evaluate the command line options
         lappend args -set $nsets
         incr nsets
         plot_scan_options $args

         #variable replot 1
         init_dataset
      }

      proc plot_update {} {
         variable datasets
         set lenx [llength [lindex [array get datasets X] 1 0]]
         set leny [llength [lindex [array get datasets Y] 1 0]]
         if {!$leny} {
            vmdcon -warn "multiplot: Data vector empty, ignoring plot!"
            variable replot 0; return
         }
         if {$lenx && $lenx!=$leny} {
            vmdcon -warn "multiplot: Different size of X and Y data, ignoring plot!"
            variable replot 0; return
         }

         variable replot 
         variable redraw

         if {!$replot && !$redraw} { return }

         # Use index number if no X-coordinate was specified
         set j 0
         foreach X $datasets(X) Y $datasets(Y) {
            if {![llength $X]} {
               set x {}
               for {set i 0} {$i<[llength $Y]} {incr i} {
                  lappend x $i
               }
               lset datasets(X) $j $x
               init_dataset
            }
            incr j
         }

         variable w
         variable c
         variable resize
         variable istoplevel

         # Display some statistics in an info frame
         variable printstats
         variable bkgcolor
         if {![winfo exists $w.info] && $printstats} { draw_infobox }
         if {[winfo exists $w.info] && !$printstats} { destroy $w.info; grid $w.cf }

         if {[winfo exists $c] && $resize} {
            variable canw
            variable canh
            $c configure -width $canw -height $canh -bg $bkgcolor
         }

         calculate_range
         calculate_ticspacing

         if {$resize} {
            # Clear the canvas
            $c addtag all all
            $c delete all
            calculate_labelsize
            calculate_plot_geometry
            draw_periphery
            redraw_objects 
            variable redraw 0
            variable resize 0
         }

         if {$replot} {
            if {$istoplevel} {
              wm deiconify $w
            }
            plot_data
            variable replot 0
         }
      }

      proc plot_data {} {
         # Plot the values
         variable c
         variable datasets
         variable marker
         # XXX initial support for NMWiz callbacks
         variable callback
         variable radius
         variable dashed
         variable fillcolor
         variable lines
         variable linewidth
         variable linecolor
         variable xmin
         variable ymin
         variable xmax
	 variable ymax
	 variable xplotmin
         variable xplotmax
         variable yplotmin
         variable yplotmax
         variable scalex
         variable scaley
         variable legend
         variable legendheight
         $c delete legend
         $c delete point
         $c delete lines

         set ds 0
         foreach X $datasets(X) Y $datasets(Y) dsxmin $datasets(xmin) dsxmin_y $datasets(xmin_y) {
            set stride 1
            #set len [llength $dset(X)]
            #if {$len>[expr 2*($xplotmax-$xplotmin)]} { 
            #   set stride [expr int([llength $dset(X)]/[expr $xplotmax-$xplotmin])]
            #   puts "Using stride $stride for data set $ds ([llength $dset(X)]/[expr $xplotmax-$xplotmin])"
            #}
            set fc   [lindex $fillcolor $ds]
            set lc   [lindex $linecolor $ds]
            set rad  [lindex $radius $ds]
            set dash [lindex $dashed $ds]
            set leg  [lindex $legend $ds]

            if {[lindex $lines $ds]} {
               set i 0
               foreach cx $X cy $Y {
                  set cxf [format "%10g" $cx]
                  set cyf [format "%10g" $cy]

                  incr i
                  if {[expr $i%$stride]} { continue }
                  set x [expr {$xplotmin + ($scalex*($cx-$xmin))}]
                  set y [expr {$yplotmin + ($scaley*($cy-$ymin))}]

                  set outofBounds 1
                  if { $cxf < $xmin } { set outofBounds [expr $outofBounds * 2] }
                  if { $cxf > $xmax } { set outofBounds [expr $outofBounds * 3] }
                  if { $cyf < $ymin } { set outofBounds [expr $outofBounds * 5] }
                  if { $cyf > $ymax } { set outofBounds [expr $outofBounds * 7] }

                  if { $i == 1 } { 
                     set oldcx $cx
                     set oldcy $cy
                     set oldx $x
                     set oldy $y
                     set oldoutofBounds $outofBounds
                     continue
                  }

                  if { $outofBounds == 1 } {
                     if { $oldoutofBounds == 1 } {
                        set item [$c create line $oldx $oldy $x $y -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     } else {
                        set xyinters [calcIntersect $oldcx $oldcy $cx $cy $oldoutofBounds]
                        set xinterplot [expr {$xplotmin + ($scalex*([lindex $xyinters 0] -$xmin))}]
                        set yinterplot [expr {$yplotmin + ($scaley*([lindex $xyinters 1] -$ymin))}]
                        set item [$c create line $xinterplot $yinterplot $x $y -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     }
                  } else {
                     if { $oldoutofBounds == 1 } {
                        set xyinters [calcIntersect $oldcx $oldcy $cx $cy $outofBounds]
                        set xinterplot [expr {$xplotmin + ($scalex*([lindex $xyinters 0] -$xmin))}]
                        set yinterplot [expr {$yplotmin + ($scaley*([lindex $xyinters 1] -$ymin))}]
                        set item [$c create line $oldx $oldy $xinterplot $yinterplot -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                        $c addtag lines withtag $item
                     } else {
                        if { ($outofBounds % 2 != 0 || $oldoutofBounds % 2 != 0) && ($outofBounds % 3 != 0 || $oldoutofBounds % 3 != 0) && \
                              ($outofBounds % 5 != 0 || $oldoutofBounds % 5 != 0) && ($outofBounds % 7 != 0 || $oldoutofBounds % 7 != 0) } {
                           set xyinters1 [calcIntersect $oldcx $oldcy $cx $cy $oldoutofBounds]
                           if {$xmin <= [lindex $xyinters1 0] && [lindex $xyinters1 0] <= $xmax && $ymin <= [lindex $xyinters1 1] && [lindex $xyinters1 1] <= $ymax} {
                              set xinterplot1 [expr {$xplotmin + ($scalex*([lindex $xyinters1 0] -$xmin))}]
                              set yinterplot1 [expr {$yplotmin + ($scaley*([lindex $xyinters1 1] -$ymin))}]
                              set xyinters2 [calcIntersect $oldcx $oldcy $cx $cy $outofBounds]
                              if {$xmin <= [lindex $xyinters2 0] && [lindex $xyinters2 0] <= $xmax && $ymin <= [lindex $xyinters2 1] && [lindex $xyinters2 1] <= $ymax} {               
                                 set xinterplot2 [expr {$xplotmin + ($scalex*([lindex $xyinters2 0] -$xmin))}]
                                 set yinterplot2 [expr {$yplotmin + ($scaley*([lindex $xyinters2 1] -$ymin))}]
                                 set item [$c create line $xinterplot1 $yinterplot1 $xinterplot2 $yinterplot2 -width [lindex $linewidth $ds] -fill $lc -dash $dash]
                                 $c addtag lines withtag $item
                              }
                           }
                        }
                     }
                  }
                  set oldcx $cx
                  set oldcy $cy
                  set oldx $x
                  set oldy $y
                  set oldoutofBounds $outofBounds
               }
            }

            if {[lindex $marker $ds]!="none"} {
               set i 0
               foreach cx $X cy $Y {
                  set cxf [format "%10g" $cx]
                  set cyf [format "%10g" $cy]
                  if { $cxf >= $xmin && $cxf <= $xmax && $cyf >= $ymin && $cyf <= $ymax } {
                     incr i
                     if {[expr $i%$stride]} { continue }
                     set x [expr {$xplotmin + ($scalex*($cx-$xmin))}]
                     set y [expr {$yplotmin + ($scaley*($cy-$ymin))}]
                     if {[string match "point*" [lindex $marker $ds]]} {
                        set item [$c create oval [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 0 -fill $fc]
                                  
                        $c addtag point withtag $item
		        $c bind $item <Any-Enter> "puts \"$cx $cy\""; # TONI - store real coordinates
                     } elseif {[lindex $marker $ds]=="circle"} {
                        set item [$c create oval [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 1 -outline $lc \
                                     -fill $fc]
                        $c addtag point withtag $item
		        $c bind $item <Any-Enter> "puts \"$cx $cy\""
                     } elseif {[lindex $marker $ds]=="square"} {
                        set item [$c create rectangle [expr {$x-$rad}] [expr {$y-$rad}] \
                                     [expr {$x+$rad}] [expr {$y+$rad}] -width 1 -outline $lc \
                                     -fill $fc]
                        $c addtag point withtag $item
		        $c bind $item <Any-Enter> "puts \"$cx $cy\""
                     }
                     # XXX initial support for NMWiz callbacks
                     if {[lindex $callback $ds]!="none"} {
                        $c bind $item <Any-ButtonPress> "[lindex $callback $ds] $i $cx $cy $fc [lindex $marker $ds]"
                     }
                  }
               }
            }

            # Draw the legend
            if {[llength $leg]} {
               variable ticfont
               variable ticfontsize
               set ylegpos [expr $yplotmax+2*$ticfontsize+$ds*2.2*$ticfontsize]
               set xlegpos [expr $xplotmin+30]
               set item [$c create line $xlegpos $ylegpos [expr $xlegpos+30] $ylegpos \
                  -width [lindex $linewidth $ds] -fill $lc -dash $dash]
               $c addtag legend withtag $item
               set item [$c create text [expr $xlegpos+30+$ticfontsize] $ylegpos -text $leg \
                            -font $ticfont -anchor w]
               $c addtag legend withtag $item
               if {[lindex $marker $ds]=="points"} {
                  set item [$c create oval [expr {$xlegpos-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+$rad}] [expr {$ylegpos+$rad}] -width 1 -fill $fc]
                  $c addtag legend withtag $item
                  set item [$c create oval [expr {$xlegpos+30-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+30+$rad}] [expr {$ylegpos+$rad}] -width 1 -fill $fc]
                  $c addtag legend withtag $item
               } elseif {[lindex $marker $ds]=="circle"} {
                  set item [$c create oval [expr {$xlegpos-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+$rad}] [expr {$ylegpos+$rad}] -width 1 -outline $lc \
                               -fill $fc]
                  $c addtag legend withtag $item
                  set item [$c create oval [expr {$xlegpos+30-$rad}] [expr {$ylegpos-$rad}] \
                               [expr {$xlegpos+30+$rad}] [expr {$ylegpos+$rad}] -width 1 -outline $lc \
                               -fill $fc]
                  $c addtag legend withtag $item
               }
            }

            incr ds
         }
#         TONI - disabled printing of mouse-derived coordinates in favour of dataset ones
#         $c bind point <Any-Enter> [namespace code {
#            #$c itemconfig current -fill red; 
#            print_datapoint %x %y
#         }]
         #$c bind point <Any-Leave> "$c itemconfig current -fill $fc"
         $c bind legend <1> "[namespace current]::grab_legend $c %x %y"
         $c bind legend <B1-Motion> "[namespace current]::move_legend $c %x %y"
      }


      proc calcIntersect { oldcx oldcy cx cy outofBounds } {
 
         variable xmin
         variable ymin
         variable xmax
         variable ymax
  
         set slope [expr {($cy - $oldcy)*1.0/($cx - $oldcx)}]
         if { $outofBounds % 2 == 0 } {
            set xinter $xmin
            set yinter [expr {$slope*($xmin - $oldcx) + $oldcy}]
         } elseif { $outofBounds % 3 == 0 } {
            set xinter $xmax
            set yinter [expr {$slope*($xmax - $oldcx) + $oldcy}]
         }
         if { $outofBounds % 5 == 0 } {
            if { $outofBounds % 2 == 0 || $outofBounds % 3 == 0 } {
               if { $yinter < $ymin } {
                  set yinter $ymin
                  set xinter [expr {($ymin - $oldcy)/$slope + $oldcx}]
               }
            } else {
               set yinter $ymin
               set xinter [expr {($ymin - $oldcy)/$slope + $oldcx}]
            }
         } elseif { $outofBounds % 7 == 0 } {
            if { $outofBounds % 2 == 0 || $outofBounds % 3 == 0 } {
               if { $yinter > $ymax } {
                  set yinter $ymax
                  set xinter [expr {($ymax - $oldcy)/$slope + $oldcx}]
               }
            } else {
               set yinter $ymax
               set xinter [expr {($ymax - $oldcy)/$slope + $oldcx}]
            }
         }
         return [list $xinter $yinter]
      }

      # Transforms coordinates from plot coordinates to canvas coords
      proc world2canvascoor {wx wy} {
         variable xplotmin
         variable yplotmin
         variable scalex
         variable scaley
         variable xmin
         variable ymin
         set x [expr {$xplotmin + ($scalex*($wx-$xmin))}]
         set y [expr {$yplotmin + ($scaley*($wy-$ymin))}]
         return [list $x $y]
      }                    
      
      proc redraw_objects {} {
        variable objectlist
        foreach object $objectlist {
           draw_object $object
        }
      }
      
      proc draw_object {object} {
        variable c
        set oname [lindex $object 0]
        set optpos [lsearch -regexp $object {^-[^[:digit:]]}]
        set options {}
        if {$optpos<0} {
          set optpos end
        } else {
          set options [lrange $object $optpos end]
          incr optpos -1
        }
        set coords [join [lrange $object 1 $optpos]]
        foreach {wx wy} $coords {
          lappend plotcoords [world2canvascoor $wx $wy]
        }
        if {$oname=="circle" || $oname=="square"} {
          set rad 1.0
          set pos [lsearch $options "-radius"]
          if {$pos>=0} {
             if {$pos+1<[llength $options]} {
                set rad [lindex $options [expr {$pos+1}]]
             }
             set options [lreplace $options $pos [expr {$pos+1}]]
          }
          foreach {x y} [join $plotcoords] {break}
          set plotcoords  [list [expr {$x-$rad}] [expr {$y-$rad}] [expr {$x+$rad}] [expr {$y+$rad}]]
          if {$oname=="circle"} { 
             set oname "oval"
          } else { set oname "rectangle" }
        }
        
        set evalstr "$c create $oname [join $plotcoords] $options"
        set item [eval $evalstr]
        $c addtag objects withtag $item
      }
      
      # grab_legend --
      # This procedure is invoked when the mouse is pressed over one of the
      # legend items.  It sets up state to allow the legend to be dragged.
      #
      # Arguments:
      # w -             The canvas window.
      # x, y -  The coordinates of the mouse press.

      proc grab_legend {w x y} {
         variable legendpos
         #$w dtag selected
         #$w addtag selected withtag current
         $w raise legend
         set legendpos(lastX) $x
         set legendpos(lastY) $y
      }

      # move_legend --
      # This procedure is invoked during mouse motion events.  It drags the
      # legend.
      #
      # Arguments:
      # w -             The canvas window.
      # x, y -  The coordinates of the mouse.

      proc move_legend {w x y} {
         variable legendpos
         $w move legend [expr {$x-$legendpos(lastX)}] [expr {$y-$legendpos(lastY)}]
         set legendpos(lastX) $x
         set legendpos(lastY) $y
      }

      proc calculate_range {} {
         # Get min/max values

         variable predefRange
         variable givenXmin
         variable givenYmin
         variable givenXmax
         variable givenYmax

         variable datasets
         set lxmin {}
         set lxmax {}
         set lymin {}
         set lymax {}
         foreach dsxmin $datasets(xmin) dsxmax $datasets(xmax) \
                 dsymin $datasets(ymin) dsymax $datasets(ymax) \
                 dsxmin_y $datasets(xmin_y) dsxmax_y $datasets(xmax_y) \
                 dsymin_x $datasets(ymin_x) dsymax_x $datasets(ymax_x) {
               lappend lxmin [list $dsxmin $dsxmin_y]
               lappend lymin [list $dsymin $dsymin_x]
               lappend lxmax [list $dsxmax $dsxmax_y]
               lappend lymax [list $dsymax $dsymax_x]
            }

         if { $predefRange } {
            if { $givenXmin == "auto" || $givenXmin == "Auto" } {
               set givenXmin [lindex [lsort -real -index 0 $lxmin] 0 0]
            }
            if { $givenXmax == "auto" || $givenXmax == "Auto" } {
               set givenXmax [lindex [lsort -real -index 0 $lxmax] end 0]
            }
            if { $givenYmin == "auto" || $givenYmin == "Auto" } {
               set givenYmin [lindex [lsort -real -index 0 $lymin] 0 0]
            }
            if { $givenYmax == "auto" || $givenYmax == "Auto" } {
               set givenYmax [lindex [lsort -real -index 0 $lymax] end 0]
            }
            if { $givenXmin < $givenXmax && $givenYmin < $givenYmax } {
               set tmpxmin $givenXmin
               set tmpymin $givenYmin
               set tmpxmax $givenXmax
               set tmpymax $givenYmax
            } else {
               variable predefRange 0
               set givenXmin auto
               set givenXmax auto
               set givenYmin auto
               set givenYmax auto
            }
         } 

         if { !$predefRange } {
            set tmpxmin [lindex [lsort -real -index 0 $lxmin] 0 0]
            set tmpymin [lindex [lsort -real -index 0 $lymin] 0 0]
            set tmpxmax [lindex [lsort -real -index 0 $lxmax] end 0]
            set tmpymax [lindex [lsort -real -index 0 $lymax] end 0]
         }

         variable xmin     
         variable ymin     
         variable xmax     
         variable ymax     
         if {$tmpxmin<$xmin || $tmpxmax>$xmax || $tmpymin<$ymin || $tmpymax>$ymax} {
            variable resize 1
         }

         variable xmin     [format "%10g" $tmpxmin]
         variable ymin     [format "%10g" $tmpymin]
         variable xmax     [format "%10g" $tmpxmax]
         variable ymax     [format "%10g" $tmpymax]
         variable spanx    [expr $xmax-$xmin]
         variable spany    [expr $ymax-$ymin]

         if { $predefRange } {
            variable xmin_y $ymin 
            variable ymin_x $xmin
            variable xmax_y $ymax
            variable ymax_x $xmax
         } else {
            variable xmin_y   [format "%10g" [lindex [lsort -real -index 0 $lxmin] 0 1]]
            variable ymin_x   [format "%10g" [lindex [lsort -real -index 0 $lymin] 0 1]]
            variable xmax_y   [format "%10g" [lindex [lsort -real -index 0 $lxmax] end 1]]
            variable ymax_x   [format "%10g" [lindex [lsort -real -index 0 $lymax] end 1]]
         }

         # Order of magnitude of value range
         if {$spanx==0.0} { variable spanx 1 } 
         if {$spany==0.0} { variable spany 1 } 
         variable dimx [expr 0.5*pow(10,floor(log10($spanx)))]
         variable dimy [expr 0.5*pow(10,floor(log10($spany)))]
      }
         
      proc calculate_ticspacing {} {
         variable spanx
         variable spany
         variable dimx
         variable dimy

         # Total number of tics between two major tics
         variable minorticx 5
         if {[expr $spanx/$dimx]>5} {
            variable minorticx 2
         }
         
         variable minorticy 5
         if {[expr $spany/$dimy]>5} {
            variable minorticy 2
         }

         variable anglescalex
         variable anglescaley
         if {$anglescalex} {
            set dimx 90
            set minorticx 3
         }
         if {$anglescaley} {
            set dimy 90
            set minorticy 3
         }

         variable xmajortics
         variable ymajortics
         variable xminortics
         variable yminortics
         if {[llength $xmajortics]} { set dimx $xmajortics }
         if {[llength $ymajortics]} { set dimy $ymajortics }
         if {[llength $xminortics]} { set minorticx $xminortics }
         if {[llength $yminortics]} { set minorticy $yminortics }
#        set i 0
#        while {1} {
#           variable loticx [expr $i*$minorticx]
#           if {$loticx<$xmin} { return [expr $i*$minorticx]}
#           incr i
#        }
         variable xmin
         variable ymin
         variable xmax
         variable ymax
         if {${::MultiPlot::verbose}} {
           vmdcon -info "dimx=$dimx xmin=$xmin xmax=$xmax ceil=[expr ceil($xmin/$dimx*$minorticx)]"
           vmdcon -info "dimy=$dimy ymin=$ymin ymax=$ymax ceil=[expr ceil($ymin/$dimy*$minorticy)]"
         }
         variable loticx [expr $dimx*ceil($xmin/$dimx*$minorticx)/$minorticx]
         variable loticy [expr $dimy*ceil($ymin/$dimy*$minorticy)/$minorticy]
         variable hiticx [expr $dimx*floor($xmax/$dimx*$minorticx)/$minorticx]
         variable hiticy [expr $dimy*floor($ymax/$dimy*$minorticy)/$minorticy]
      }

      proc calculate_labelsize {} {
         # Measure y-axis label size
         variable c
         variable ylabelfont
         variable ylabelfontsize
         variable xlabelfont
         variable xlabelfontsize
         variable xlabeltext
         variable ylabeltext
         variable ylabelorient
         if {[llength $ylabeltext]} {
            set item [$c create text 0 0 -text $ylabeltext -font "$ylabelfont $ylabelfontsize" -anchor nw]
            set bbox [$c bbox $item]
            variable ylabelheight [expr [lindex $bbox 3]-[lindex $bbox 1]]
            variable ylabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0] + $ylabelheight]
            $c delete $item
         } else {
            variable ylabelwidth 0.0
         }

         # Measure x-axis label height
         if {[llength $xlabeltext]} {
            set item [$c create text 0 0 -text $xlabeltext -font "$xlabelfont $xlabelfontsize" -anchor nw]
            set bbox [$c bbox $item]
            $c delete $item
            variable xlabelheight  [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         } else {
            variable xlabelheight 0.0
         }

         
         ## Measure x-axis ticlabel size
         variable loticx
         variable hiticx
         variable ticfont
         # Compare smallest and biggest tics
         set absxmax [lindex [lsort -real [list [expr abs($loticx)] [expr abs($hiticx)]]] end]
         set item [$c create text 0 0 -text [format "-%g" $absxmax] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable ticlabelheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         variable xticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]

         ## Measure y-axis ticlabel size
         variable dimx
         variable loticy
         variable hiticy
         variable ticfont
         # Compare smallest and biggest tics
         set absymax [lindex [lsort -real [list [expr abs($loticy)] [expr abs($hiticy)]]] end]
         set item [$c create text 0 0 -text [format "-%g" $absymax] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable ticlabelheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
         variable yticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]
         # Check if the neighboring ticlabel is wider since it could involve more decimal places
         set item [$c create text 0 0 -text [format "-%g" [expr $absymax+$dimx]] -font $ticfont -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         if {[expr 1.5*[lindex $bbox 2]-[lindex $bbox 0]]>$yticlabelwidth} {
            variable yticlabelwidth [expr [lindex $bbox 2]-[lindex $bbox 0]]
         }

         # Measure title height
         variable title
         variable titlefont
         variable titlefontsize
         if {![llength $title]} { set title [namespace current] }
         set item [$c create text 0 0 -text $title -font "$titlefont $titlefontsize" -anchor nw]
         set bbox [$c bbox $item]
         $c delete $item
         variable titleheight [expr 1.5*[lindex $bbox 3]-[lindex $bbox 1]]
      } 
      
      proc calculate_plot_geometry {} {
         # Compute legend height
         variable legend         
         variable ticfontsize
         variable legendheight 0.0
         foreach legitem $legend {
            if {[llength $legitem]} {
               set legendheight [expr $legendheight+1.8*$ticfontsize]
            }
         }
         
         ## Plot geometry
         variable rim
         variable canh
         variable canw
         variable ticlen
         variable xlabelheight
         variable ylabelwidth
         variable xticlabelwidth
         variable yticlabelwidth
         variable ticlabelheight
         variable titleheight
         variable xplotmin [expr $rim+$ylabelwidth+$yticlabelwidth+$ticlen]
         variable yplotmin [expr $canh-($rim+$xlabelheight+$ticlabelheight+$ticlen)]
         variable xplotmax [expr $canw-$rim-0.5*$xticlabelwidth]
         variable yplotmax [expr $rim+$titleheight]

         # Scaling factor to convert world coordinates into plot coordinates
         variable spanx
         variable spany
         variable scalex [expr ($xplotmax-$xplotmin)/(1.0*$spanx)]      
         variable scaley [expr ($yplotmax-$yplotmin)/(1.0*$spany)]      

         variable dimx
         variable scalex
         if {[expr $xticlabelwidth]>[expr $dimx*$scalex*0.7]} {
            set dimx [expr 2.0*$dimx]
            calculate_ticspacing
            calculate_labelsize
            calculate_plot_geometry
         }
      }

      proc draw_periphery {} {
         # Draw title
         variable c
         variable rim
         variable canw
         variable canh
         variable title
         variable titlefont
         variable titlefontsize
         variable bkgcolor
         $c create text [expr $canw/2] $rim -anchor n -text $title -font "$titlefont $titlefontsize" -fill brown

         # Draw bounding box
         variable xplotmin
         variable yplotmin
         variable xplotmax
         variable yplotmax
         $c create line $xplotmin $yplotmin $xplotmax $yplotmin -width 2
         $c create line $xplotmin $yplotmin $xplotmin $yplotmax -width 2
         $c create line $xplotmax $yplotmin $xplotmax $yplotmax -width 2
         $c create line $xplotmin $yplotmax $xplotmax $yplotmax -width 2

         # Scaling factor to convert plot coordinates into canvas coordinates
         variable spanx
         variable spany
         variable scalex        
         variable scaley
         variable xmin  
         variable ymin  
         variable xmax  
         variable ymax  

         # x-axis, y=0
         if {$ymin<0 && $ymax>0} {
            set zero [expr $yplotmin-($scaley*$ymin)]
            $c create line $xplotmin $zero $xplotmax $zero -width 1 -dash -
         }
         # y-axis, x=0
         if {$xmin<0 && $xmax>0} {
            set zero [expr $xplotmin-($scalex*$xmin)]
            $c create line $zero $yplotmin $zero $yplotmax -width 1 -dash -
         }

         # x-label
         variable ticlen
         variable xlabelfont
         variable xlabelfontsize
         variable xlabeltext
         variable xlabelheight
         if {[llength $xlabeltext]} {
            set labelposx [expr $xplotmin+($xplotmax-$xplotmin)*0.5]
            #set labelposy [expr $yplotmin+$ticlen+$ticlabelheight+0.2*$xlabelheight]
            set labelposy [expr $canh-$rim-0.2*$xlabelheight]
            $c create text $labelposx $labelposy -text $xlabeltext -font "$xlabelfont $xlabelfontsize" -anchor s
         }

         # y-label
         variable ylabelfont
         variable ylabelfontsize
         variable ylabeltext
         if {[llength $ylabeltext]} {
            set labelposy [expr $yplotmin+($yplotmax-$yplotmin)*0.5]
            $c create text $rim $labelposy -text $ylabeltext -font "$ylabelfont $ylabelfontsize" -anchor w
         }

         # Draw x-tics
         variable ticfont
         set i 0 
         set ticval $xmin
         variable dimx
         variable hiticx
         variable loticx
         variable minorticx
         variable ticlabelheight
         set firstmajor [expr abs(int($loticx-int($loticx/$minorticx)*$minorticx))]
         #set firstmajor [expr $loticx%$minorticx]
         set ticlabelposy [expr $yplotmin+$ticlen+0.2*$ticlabelheight]
         while {$ticval<$hiticx} {
            set ticval [expr $loticx+$i*$dimx/$minorticx]
            set x [expr $xplotmin + ($ticval-$xmin)*$scalex]
            if {![expr ($i-$firstmajor)%$minorticx]} {
               $c create line $x $yplotmin $x [expr $yplotmin+$ticlen] -width 2
               $c create text $x $ticlabelposy -text [format "%g" $ticval] -anchor n -font $ticfont
            } else {
               $c create line $x $yplotmin $x [expr $yplotmin+0.5*$ticlen] -width 2
            }
            incr i
         }

         # Draw y-tics
         set i 0
         variable dimy
         variable hiticy
         variable loticy
         variable minorticy
         set firstmajor [expr abs(int($loticy-int($loticx/$minorticy)*$minorticy))]
         set ticlabelposx [expr $xplotmin-$ticlen-0.2*$ticlabelheight]
         set ticval $ymin
         while {$ticval<$hiticy} {
            set ticval [expr $loticy+$i*$dimy/$minorticy]
            set y [expr $yplotmin + ($ticval-$ymin)*$scaley]
            if {![expr ($i-$firstmajor)%$minorticy]} {
               $c create line $xplotmin $y [expr $xplotmin-$ticlen] $y -width 2
               $c create text $ticlabelposx $y -text [format "%g" $ticval] -anchor e -font $ticfont
            } else {
               $c create line $xplotmin $y [expr $xplotmin-0.5*$ticlen] $y -width 2
            }
            incr i
         }

         # Draw user specified horizontal lines
         variable hline
         foreach line $hline {
            set y [lindex $line 0]
            set zero [expr $yplotmin-($scaley*$ymin)]
            set opt [lrange $line 1 end]
            set ypos [expr $yplotmin+($scaley*($y-$ymin))]
            if {${::MultiPlot::verbose}} {
              vmdcon -info "$c create line $xplotmin $ypos $xplotmax $ypos $opt"
            }
            eval $c create line $xplotmin $ypos $xplotmax $ypos $opt
         }

         # Draw user specified vertical lines
         variable vline
         foreach line $vline {
            set x [lindex $line 0]
            set opt [lrange $line 1 end]
            set xpos [expr $xplotmin+($scalex*($x-$xmin))]
            eval $c create line $xpos $yplotmin $xpos $yplotmax $opt
         }
      }

      proc print_datapoint {x y} {
         variable xplotmin
         variable yplotmin
         variable scalex
         variable scaley
         variable xmin
         variable ymin
         set coords [format "%8g %8g" [expr ($x-$xplotmin)/$scalex+$xmin] [expr ($y-$yplotmin)/$scaley+$ymin]]
         puts $coords
      }

      proc draw_infobox {} {
         variable w
         variable infoFont
         labelframe $w.info -text "Info"
         label $w.info.headx -text [format "%10s" "x"] -font $infoFont
         label $w.info.heady -text [format "%10s" "y"] -font $infoFont
         grid $w.info.headx -row 1 -column 2
         grid $w.info.heady -row 1 -column 3
         label $w.info.xmint -text "X min: "       -font $infoFont
         label $w.info.xmin  -textvariable [namespace current]::xmin   -font $infoFont
         label $w.info.xminy -textvariable [namespace current]::xmin_y -font $infoFont
         grid $w.info.xmint  -row 2 -column 1
         grid $w.info.xmin   -row 2 -column 2
         grid $w.info.xminy  -row 2 -column 3
         label $w.info.xmaxt -text "X max: "       -font $infoFont
         label $w.info.xmax  -textvariable [namespace current]::xmax   -font $infoFont
         label $w.info.xmaxy -textvariable [namespace current]::xmax_y -font $infoFont
         grid $w.info.xmaxt  -row 3 -column 1
         grid $w.info.xmax   -row 3 -column 2
         grid $w.info.xmaxy  -row 3 -column 3
         label $w.info.ymint -text "Y min: "       -font $infoFont
         label $w.info.ymin  -textvariable [namespace current]::ymin   -font $infoFont
         label $w.info.yminx -textvariable [namespace current]::ymin_x -font $infoFont
         grid $w.info.ymint  -row 4 -column 1
         grid $w.info.ymin   -row 4 -column 2
         grid $w.info.yminx  -row 4 -column 3
         label $w.info.ymaxt -text "Y max: "       -font $infoFont
         label $w.info.ymax  -textvariable [namespace current]::ymax   -font $infoFont
         label $w.info.ymaxx -textvariable [namespace current]::ymax_x -font $infoFont
         grid $w.info.ymaxt  -row 5 -column 1
         grid $w.info.ymax   -row 5 -column 2
         grid $w.info.ymaxx  -row 5 -column 3
         #pack  $w.info -side top -pady 2m -ipadx 5m -ipady 2m
         grid  $w.info
      }

      proc savedialog {} {
         variable w
         set types {
            {{Postscript files} {.ps}}
            {{All files}        *   }
         }
         variable postscript
         set newfile [tk_getSaveFile \
                             -title "Choose file name" -parent $w \
                             -initialdir [pwd] -filetypes $types -initialfile $postscript]
         variable c
         if {[llength $newfile]} {
            $c postscript -file $newfile
            vmdcon -info "Wrote plot postscript file to $newfile."
         }
      }

      proc xmgracedialog {} {
         variable w
         set types {
            {{xmgrace files} {.agr}}
            {{All files}        *  }
         }
         set newfile [tk_getSaveFile \
                      -title "Choose file name" -parent $w \
                     -initialdir [pwd] -filetypes $types -initialfile "multiplot.agr"]
         if {[llength $newfile]} {
           [namespace current]::plothandle export xmgrace $newfile
           vmdcon -info "Wrote plot to $newfile, running xmgrace..."
           ::ExecTool::exec xmgrace $newfile &
         }
      }

      proc savevectors {} {   # TONI
         variable w
         set types {
            {{ASCII data files} {.dat}}
            {{All files}        *  }
         }
         set newfile [tk_getSaveFile \
                      -title "Choose file name" -parent $w \
                      -initialdir [pwd] -filetypes $types -initialfile "multiplot.dat"]
         if {[llength $newfile]} {
           [namespace current]::plothandle export vectors $newfile
         }
      }

      # return 1 if all datasets have the same X vector
      proc cansavematrix {} {
	  variable datasets
	  variable nsets
	  set x0 [lindex $datasets(X) 0]; # X of set 0
	  set x0l [llength $x0]
	  for {set s 0} {$s<$nsets} {incr s} {
	      if {[llength [lindex $datasets(X) $s]] != $x0l} { return 0 }
	      for {set i 0} {$i<$x0l} {incr i} {
		  if { [lindex $x0 $i] != [lindex $datasets(X) $s $i] } {
		      return 0
		  }
	      }
	  }
	  return 1
      }

      proc savematrix {} {	# TONI
         variable w
	  if {![cansavematrix] } {
	      tk_messageBox -icon error -message "All of the data sets should have the same X vector to be saved in matrix form"
	      return
	  }
         set types {
            {{ASCII data files} {.dat}}
            {{All files}        *  }
         }
         set newfile [tk_getSaveFile \
                      -title "Choose file name" -parent $w \
                      -initialdir [pwd] -filetypes $types -initialfile "multiplot.dat"]
         if {[llength $newfile]} {
           [namespace current]::plothandle export matrix $newfile
         }
      }

      proc editplot {} {
         puts "sending [namespace current] to editPlot"
         ::MultiPlot::editPlot [namespace current]
      }

      proc managestyle {mode} {
         ::MultiPlot::manageStyle $mode [namespace current] 
      }

   } ; # END namespace $ns

   return "::MultiPlot::Plot${::MultiPlot::plotcount}::plothandle"
}

proc multiplot { args } {
   set keyword [lindex $args 0]
   if {![llength $keyword]} { return }
   if {$keyword=="list"} {
      set plist {}
      foreach plot [namespace children ::MultiPlot "Plot*"] { 
         lappend plist [subst $plot]::plothandle
      }
      return $plist
   } elseif {$keyword=="reset"} {
      # on reset we may delete only toplevel widgets
      foreach ploth [namespace children ::MultiPlot "Plot*"] {
         eval "set tl $[subst $ploth]::istoplevel"
         if { $tl } {
            if {$::MultiPlot::verbose} {
              vmdcon -info "deleting toplevel widget of $ploth"
            }
            eval "destroy $[subst $ploth]::w"
            namespace delete $ploth
         }
      }
      return
   }

   variable plothandle
   if {$keyword=="embed"} {
     set parent [lindex $args 1]
     if {$parent == ""} {
        return
     } else {
        set plothandle [::MultiPlot::init_plot $parent]
     }
   } else {
     set plothandle [::MultiPlot::init_plot ]
   } 
   #puts "$plothandle configure $args"
   eval $plothandle configure $args

   return $plothandle
}

################################################################################
### Testing Adding QMOrbital plot from QwikMD here to use in QMtool
################################################################################
################################################################################
### PLot orbitals energies for the current frame
################################################################################
proc ::MultiPlot::plotOrbitals {plugin} {
    
    puts "In plotOrbitals. plugin = $plugin"
    array set orb_frame ""
    #array set orb_xvals ""
    #array set orb_evals ""
    #array set orb_frame ""
    
    if {$plugin == "QWIKMD"} {
      set qmID [expr $QWIKMD::curframe +1]
      set phandle QWIKMD::advGui(analyze,advance,plot,qmorbital,$qmID)
      set plabel QWIKMD::advGui(analyze,advance,plot,qmorbital,$qmID,label)
      set ptable QWIKMD::advGui(analyze,qmorb,table)
      set qmorbmol QWIKMD::qmorbmol
      set qmtslist QWIKMD::qmtslist
      #set orb_frame $QWIKMD::advGui(analyze,advance,plot,qmorbital,$qmID,frame)
      set orbselected QWIKMD::advGui(analyze,advance,plot,qmorbital,$qmID,orbselected)
    } elseif {$plugin == "QMtool"} {    
      set qmID $::QMtool::molid
      set phandle ::QMtool::plothandle2
      set qmorbmol ::QMtool::molid
      set ptable ::QMtool::qmorbtable
      set qmtslist ::QMtool::qmtslist
      set plabel ::QMtool::qmplabel
      set orbselected ::QMtool::qmorbselected
    }
    puts "Using new plotorbitals"
    
    set orb_frame($qmID) [list]
    
    
    
    if {[info exists $phandle] == 0 ||
        [set $phandle] == ""} {
        if {[QWIKMD::callQMOrbitals] < 0 } {
            return
        }
    }

    set orblist [[set $ptable] getcolumns 0]
    set descrlist [[set $ptable] getcolumns 1]
    set energylist [[set $ptable] getcolumns 2]

    set curframe [molinfo [set $qmorbmol] get frame]
    set tmstep [lindex [lindex [set $qmtslist] 0 ] 1]
    set i 0
    while {$curframe > [lindex [lindex [set $qmtslist] $i] 0] && $i < [llength [set $qmtslist]]} {
        set tmstep [lindex [lindex [set $qmtslist] $i] 1]
        incr i
    } 

    if {$i > [llength [set $qmtslist]]} {
        set tmstep [lindex [lindex [set $qmtslist] [expr $i -1] ] 1]
    }
    set xval [expr $curframe * $tmstep]
    if {[lsearch $orb_frame($qmID) $curframe] == -1} {
        foreach energy $energylist orb $orblist descr $descrlist {
            lappend orb_xvals($qmID) $xval
            lappend orb_evals($qmID) $energy
            lappend orb_orb($qmID) $orb
            lappend orb_frame($qmID) $curframe
            if {$descr == "HOMO"} {
                lappend orb_homoframelist($qmID) [list $xval $orb $energy]
            }
        }
    }

    set lframes [lsort -unique -real $orb_xvals($qmID)]
    set framdist 0.25
    if {[llength $lframes] > 1} {
        ### Set frame dist to evaluate the biggest width for energy lines
        ### framdist is so big in the beginning because we want to find the smallest different 
        ### and use 25% of the difference as width
        
        set framdist 999999
        for {set i 1} {$i < [llength $lframes]} {incr i} {
            set dist [expr [lindex $lframes $i] - [lindex $lframes [expr $i -1]] ]
            if {$dist < $framdist} {
                set framdist $dist
            }
        }
       
        set framdist [expr abs($framdist * 0.25)]
    }
    set tbind [[set $ptable] curselection]
    if {[llength $tbind] == 2 && [lsearch [set $orbselected] "$xval [lindex $energylist [lindex $tbind 0]] [lindex $orblist [lindex $tbind 1]] * * *"] == -1 && [lsearch [set $orbselected] "$xval [lindex $energylist [lindex $tbind 1]] [lindex $orblist [lindex $tbind 0]] * * *"]} {
        ### orbselected stores {xmin ymin xmax ymax} of the selected orbitals
        set orb1 [lindex $orblist [lindex $tbind 0]]
        set orb2 [lindex $orblist [lindex $tbind 1]]
        set ene1 [lindex $energylist [lindex $tbind 0]]
        set ene2 [lindex $energylist [lindex $tbind 1]]
        if {$orb2 < $orb1} {
            set aux $ene1
            set ene1 $ene2
            set ene1 $aux
        }
        lappend $orbselected [list $xval $ene1 $ene2 $tbind $curframe $framdist]
        
    }
    


    set xmin [expr [lindex $lframes 0] - 2 * $framdist ]
    set xmax [expr [lindex $lframes end] + 2 * $framdist]

    set ymin [expr [QWIKMD::mincalc $orb_evals($qmID)] - 0.1]
    set ymax [expr [QWIKMD::maxcalc $orb_evals($qmID)] + 0.1]

    ### Configure plot axis
    [set $phandle] configure -xmin $xmin -xmax $xmax -ymin $ymin -ymax $ymax
    
    ### Clean up lines and data points
    [set $phandle] clear
    ::MultiPlot::deletePlotLines $phandle $plabel

    set preframe [lindex $orb_frame($qmID) 0]
    set xframe [list]
    set yenergy [list]
    
    foreach energy $orb_evals($qmID) \
    x $orb_xvals($qmID) \
    orb $orb_orb($qmID) \
    frame $orb_frame($qmID) {
        set xminm [expr $x - $framdist]
        set xmaxi [expr $x + $framdist]

        if {$frame != $preframe } {
            [set $phandle] add $xframe $yenergy -marker none -nolines
            set index [lsearch -all -index 0 [set $orbselected] [lindex $xframe 0] ]
            foreach ind $index {
                ::MultiPlot::addQMOrbDeltaLine [lindex [set $orbselected] $ind] $framdist
            }
            
            set xframe $x
            set yenergy $energy
            set preframe $frame
        } else {
            lappend xframe $x
            lappend yenergy $energy
        }
        set fill black
        if {[lsearch $orb_homoframelist($qmID) "* $orb $energy" ] != -1} {
            set fill blue
        }
        [set $phandle] draw line $xminm $energy $xmaxi $energy -width 2 -fill $fill -tags ${orb}orb$frame
        
        ### Add the binding command to the lines when pressed with the left click mouse
        lappend $plabel ${orb}orb$frame
        
        set text "animate goto $frame;::${plugin}::updateOrbitalsTable;\[set $ptable\] selection clear 0 end;\[set $ptable\] selection set \[lsearch \[\[set $ptable\] getcolumns 0\] $orb\];::${plugin}::showOrbitals"
        eval "$[[set $phandle] namespace]::c bind ${orb}orb$frame <Any-ButtonPress> \"$text\""
        
        
    }

    ### connect the homo orbitals energy values
    if {[llength $orb_homoframelist($qmID)] > 0} {

        set homoframelist [lsort -index 0 $orb_homoframelist($qmID)]
        for {set i 1} {$i < [llength $homoframelist]} {incr i} {
            set orbframe1 [lindex $homoframelist [expr $i -1]]
            set orbframe2 [lindex $homoframelist $i]

            set frame1 [lindex $orbframe1 0]
            set orb1 [lindex $orbframe1 1]
            set energy1 [lindex $orbframe1 2]

            set frame2 [lindex $orbframe2 0]
            set orb2 [lindex $orbframe2 1]
            set energy2 [lindex $orbframe2 2]
           
            set xmin [expr $frame1 + $framdist]
            set ymin $energy1

            set xmax [expr $frame2 - $framdist]
            set ymax $energy2

            [set $phandle] draw line $xmin $ymin $xmax $ymax -width 2 -fill blue -dash - -tags ${frame1}homo$frame2

            lappend $plabel ${frame1}homo$frame2
        }
    }

    ### Add points and text to the plot
    if {[llength $xframe] >0} {
        [set $phandle] add $xframe $yenergy -marker none -nolines
        set index [lsearch -all -index 0 [set $orbselected] [lindex $xframe 0] ]
        foreach ind $index {
            ::MultiPlot::addQMOrbDeltaLine [lindex [set $orbselected] $ind] $framdist
        }
        set xframe [list]
        set yenergy [list]
    }

    [set $phandle] replot
}
################################################################################
### Add delta lines to the plot
################################################################################
proc ::MultiPlot::addQMOrbDeltaLine {orbselected framdistt} {
    variable plugin
    variable ptable
    if {$plugin == "QWIKMD"} {
       set phandle QWIKMD::advGui(analyze,advance,plot,qmorbital,[expr $QWIKMD::curframe +1])
       set plabel QWIKMD::advGui(analyze,advance,plot,qmorbital,[expr $QWIKMD::curframe +1],label)
    } elseif {$plugin == "QMtool"} {
       set phandle ::QMtool::plothandle2
       set ptable ::QMtool::qmorbtable
    }    
    set x [lindex $orbselected 0]
    set ymin [lindex $orbselected 2]
    set ymax [lindex $orbselected 1]
    set y [expr $ymax - $ymin]
    set tbind [lindex $orbselected 3]
    set frame [lindex $orbselected 4]
    set framdist [lindex $orbselected 5]

    [set $phandle] draw line $x $ymin $x $ymax -width 4 -fill orange -tags deltaline$frame

    lappend $plabel deltaline$frame
    set command "animate goto $frame;::${plugin}::updateOrbitalsTable;$ptable selection clear 0 end;\$ptable selection set [list $tbind];::${plugin}::showOrbitals"
    eval "$[[set $phandle] namespace]::c bind deltaline$frame <Any-ButtonPress> \"$command\""

    set text [QWIKMD::format5Dec [expr ([lindex $orbselected 2] - [lindex $orbselected 1])]]
    set position [expr ([lindex $orbselected 2] + [lindex $orbselected 1])/2]
    [set $phandle] draw text $x $position -fill black -tags ${ymin}y${ymax}text$frame -text $text -justify center -anchor w

    lappend $plabel ${ymin}y${ymax}text$frame

    eval "$[[set $phandle] namespace]::c bind ${ymin}y${ymax}text$frame <Any-ButtonPress> \"$command\""
    
}

################################################################################
### Delete lines from plot of the orbitals
################################################################################
proc ::MultiPlot::deletePlotLines {phandle plabel} {
    if {[llength [set $plabel]] > 0} {
        foreach label [set $plabel]) {
            [set $phandle] undraw $label
            eval "$[[set $phandle] namespace]::c delete $label"
            set [[set $phandle] namespace]::objectlist {}
        }
        set $plabel [list]
    }
}


##################################################################################
### Idea: Add GUI menu to change plot atributes
### independet of namespaces/plothandles, a single instance
### that can select the plothandle to work on
### OBS: need to add code to fetch title/labels/properties of plots
### Author: Mariano Spivak
##################################################################################


proc ::MultiPlot::resizeGUI {w} {
  # taken from fftk
  update idletasks
  regexp {([0-9]+)x[0-9]+[\+\-]+[0-9]+[\+\-]+[0-9]+} [wm geometry $w] all dimW
  set dimH [winfo reqheight $w]
  set dimH [expr {$dimH + 10}]
  set dimW [winfo reqwidth $w]
  set dimW [expr {$dimW + 5}]
  wm geometry $w [format "%ix%i" $dimW $dimH]
  update idletasks

}


proc ::MultiPlot::editPlot {ns} { 
   variable editlist
   variable editselected
   variable editlastselected
   variable editdataset
   variable editdatasetlist
   variable editplotysize
   variable editplotxsize

   variable edittitle
   variable edittitlefont
   variable edittitlefontsize
   variable editlinedash
   variable editlinewidth
   variable editlinecolor
   variable editmarkersize
   variable editmarkercolor
   variable editmarkertype

   variable editylabel
   variable editylabelfont
   variable editylabelfontsize
   variable editxlabel
   variable editxlabelfont
   variable editxlabelfontsize

   variable editgridy
   variable editgridx
   
   variable editymin
   variable editymax
   variable editxmin
   variable editxmax
   
   # check window was present, if already initialized just turn on
   set v .editmultiplot
   if { [winfo exists $v] } {
      wm deiconify $v
      raise $v
      return
   }

   toplevel $v
   wm title $v "Edit Plot"
   wm resizable $v 1 1
   set width 240 ;# in pixels
   set height 358 ;# in pixels
   wm geometry $v ${width}x${height}+68+125
   grid columnconfigure $v 0 -weight 1
   grid columnconfigure $v 1 -weight 0
   grid rowconfigure $v 0 -weight 0
   grid rowconfigure $v 1 -weight 1

   grid [frame $v.options] -row 0 -column 0 -sticky news
   #grid columnconfigure $w.options 0 -weight 1
   #grid rowconfigure $w.options 0 -weight 1

   set editlist {}
   foreach s [multiplot list] {
      lappend editlist "::[lindex [split $s ":"] 2]::[lindex [split $s ":"] 4]"
   }

   set editselected $ns
   set editlastselected $editselected

   # find the number of datasets in plot
   set nsets [${editselected}::plothandle nsets]
   # create a list from 0 to nsets-1
   set editdatasetlist {}
   set editdataset 0
   for {set i 0} {$i < $nsets} {incr i} {lappend editdatasetlist $i}
   
   grid [labelframe $v.options.handle -text "Select Plot"] -row 0 -column 0 -sticky news -padx 4 -pady 2
   grid [ttk::combobox $v.options.handle.sel -values $::MultiPlot::editlist -width 34 -state readonly -justify left -textvariable ::MultiPlot::editselected] -row 0 -column 0 -columnspan 4 -sticky ew -padx 2 -pady 2
   bind .editmultiplot.options.handle.sel <<ComboboxSelected>> {
      ::MultiPlot::editSelect
      return
   }
   grid [label $v.options.handle.ysizelbl -text "Heigth" -width 10] -row 1 -column 0 -sticky news -padx 1 -pady 2
   grid [spinbox $v.options.handle.ysize -textvariable ::MultiPlot::editplotysize -width 4 -from 100 -to 2000 -increment 100 -command {::MultiPlot::editPeriphery size}] -row 1 -column 1 -sticky news -padx 1 -pady 2
   grid [label $v.options.handle.xsizelbl -text "Width" -width 10] -row 1 -column 2 -sticky news -padx 1 -pady 2
   grid [spinbox $v.options.handle.xsize -textvariable ::MultiPlot::editplotxsize -width 4 -from 100 -to 2000 -increment 100 -command {::MultiPlot::editPeriphery size}] -row 1 -column 3 -sticky news -padx 1 -pady 2

   set wd2 10
   set px 2
   grid [labelframe $v.options.l -text "Title Options"] -row 2 -column 0 -sticky news -padx 4 -pady 2
   grid columnconfigure $v.options.l 1 -weight 1
   grid [entry $v.options.l.tentry -textvariable ::MultiPlot::edittitle -width 18] -row 0 -column 0 -sticky news -padx 1 -pady 2
   bind .editmultiplot.options.l.tentry <Return> {
      ::MultiPlot::editPeriphery title
      return
   }
   set fonttype "Helvetica Courier Times"
   grid [ttk::combobox $v.options.l.font -values $fonttype -width 8 -state readonly -justify left -textvariable ::MultiPlot::edittitlefont] -row 0 -column 1 -sticky ns -padx 1 -pady 2
   bind .editmultiplot.options.l.font <<ComboboxSelected>> {
      ::MultiPlot::editPeriphery title
      return
   }
   grid [spinbox $v.options.l.size -textvariable ::MultiPlot::edittitlefontsize -width 4 -from 10 -to 50 -increment 5 -command {::MultiPlot::editPeriphery title}] -row 0 -column 2 -sticky news -padx 1 -pady 2

   set ew 5
   set wd 12
   grid [labelframe $v.options.a -text "Axis Options"] -row 3 -column 0 -sticky news -padx 4 -pady 2
   grid [label $v.options.a.ylbl -text "Y label" -width 5] -row 0 -column 0 -sticky news -padx 2 -pady 2
   grid [entry $v.options.a.yentry -textvariable ::MultiPlot::editylabel -width $wd] -row 0 -column 1 -sticky news -padx 1 -pady 2
   bind .editmultiplot.options.a.yentry <Return> {
      ::MultiPlot::editPeriphery ylabel
      return
   }
   grid [ttk::combobox $v.options.a.ylabelfont -values $fonttype -width 8 -state readonly -justify left -textvariable ::MultiPlot::editylabelfont] -row 0 -column 2 -sticky ns -padx 1 -pady 2
   bind .editmultiplot.options.a.ylabelfont <<ComboboxSelected>> {
      ::MultiPlot::editPeriphery ylabel
      return
   }
   grid [spinbox $v.options.a.ylabelsize -textvariable ::MultiPlot::editylabelfontsize -width 4 -from 10 -to 50 -increment 5 -command {::MultiPlot::editPeriphery ylabel}] -row 0 -column 3 -sticky news -padx 1 -pady 2
   
   grid [label $v.options.a.xlbl -text "X label" -width 5] -row 1 -column 0 -sticky news -padx 2 -pady 2
   grid [entry $v.options.a.xentry -textvariable ::MultiPlot::editxlabel -width $wd] -row 1 -column 1 -sticky news -padx 1 -pady 2
   bind .editmultiplot.options.a.xentry <Return> {
      ::MultiPlot::editPeriphery xlabel
      return
   }
   grid [ttk::combobox $v.options.a.xlabelfont -values $fonttype -width 8 -state readonly -justify left -textvariable ::MultiPlot::editxlabelfont] -row 1 -column 2 -sticky ns -padx 1 -pady 2
   bind .editmultiplot.options.a.xlabelfont <<ComboboxSelected>> {
      ::MultiPlot::editPeriphery xlabel
      return
   }
   grid [spinbox $v.options.a.xlabelsize -textvariable ::MultiPlot::editxlabelfontsize -width 4 -from 10 -to 50 -increment 5 -command {::MultiPlot::editPeriphery xlabel}] -row 1 -column 3 -sticky news -padx 1 -pady 2
   
   grid [ttk::separator $v.options.a.s -orient horizontal] -row 2 -column 0 -columnspan 5 -sticky ew -padx 1 -pady 1

   grid [frame $v.options.a.minmax] -row 3 -column 0 -columnspan 5 -padx 1 -pady 1
   grid [label $v.options.a.minmax.yminlbl -text "Y min" -width 5] -row 2 -column 0 -sticky news
   grid [entry $v.options.a.minmax.ymin -textvariable ::MultiPlot::editymin -width $ew] -row 2 -column 1 -sticky news
   bind .editmultiplot.options.a.minmax.ymin <Return> {
      ::MultiPlot::editPeriphery ymin
      return
   }
   grid [label $v.options.a.minmax.ymaxlbl -text "Y max" -width 5] -row 2 -column 2 -sticky news
   grid [entry $v.options.a.minmax.ymax -textvariable ::MultiPlot::editymax -width $ew] -row 2 -column 3 -sticky news
   bind .editmultiplot.options.a.minmax.ymax <Return> {
      ::MultiPlot::editPeriphery ymax
      return
   }
   grid [button $v.options.a.minmax.yauto -text "Auto Y" -width 6 -command {::MultiPlot::editPeriphery yauto; ::MultiPlot::editSelect}] -row 2 -column 4 -sticky news
   
   grid [label $v.options.a.minmax.xminlbl -text "X min" -width 5] -row 3 -column 0 -sticky news
   grid [entry $v.options.a.minmax.xmin -textvariable ::MultiPlot::editxmin -width $ew] -row 3 -column 1 -sticky news
   bind .editmultiplot.options.a.minmax.xmin <Return> {
      ::MultiPlot::editPeriphery xmin
      return
   }
   grid [label $v.options.a.minmax.xmaxlbl -text "X max" -width 5] -row 3 -column 2 -sticky news
   grid [entry $v.options.a.minmax.xmax -textvariable ::MultiPlot::editxmax -width $ew] -row 3 -column 3 -sticky news
   bind .editmultiplot.options.a.minmax.xmax <Return> {
      ::MultiPlot::editPeriphery xmax
      return
   }
   grid [button $v.options.a.minmax.xauto -text "Auto X" -width 6 -command {::MultiPlot::editPeriphery xauto; ::MultiPlot::editSelect}] -row 3 -column 4 -sticky news
   grid [ttk::separator $v.options.a.s2 -orient horizontal] -row 4 -column 0 -columnspan 5 -sticky ew -padx 1 -pady 1
   grid [frame $v.options.a.grid] -row 5 -column 0 -columnspan 5 -padx 1 -pady 1
   set gridoptions "none major all"
   grid [label $v.options.a.grid.ylbl -text "Grid Y axis" -width 15 -anchor w] -row 0 -column 0 -sticky news  -padx 2 -pady 0
   grid [ttk::combobox $v.options.a.grid.ycomb -values $gridoptions -width 10 -state readonly -justify left -textvariable ::MultiPlot::editgridy] -row 0 -column 1 -sticky ns -padx 2 -pady 0
   bind .editmultiplot.options.a.grid.ycomb <<ComboboxSelected>> {
      set text [%W get]
      ::MultiPlot::editGrid y $text
      return
   }
   grid [label $v.options.a.grid.xlbl -text "Grid X axis" -width 15 -anchor w] -row 1 -column 0 -sticky news -padx 2 -pady 0
   grid [ttk::combobox $v.options.a.grid.xcomb -values $gridoptions -width 10 -state readonly -justify left -textvariable ::MultiPlot::editgridx] -row 1 -column 1 -sticky ns -padx 2 -pady 0
   bind .editmultiplot.options.a.grid.xcomb <<ComboboxSelected>> {
      set text [%W get]
      ::MultiPlot::editGrid x $text
      return
   }
   grid [labelframe $v.options.b -text "Plot Options"] -row 4 -column 0 -sticky news -padx 4 -pady 2
   grid [label $v.options.b.labelset -text "Select Dataset" -width 15 -anchor w] -row 0 -column 0 -columnspan 2 -sticky news -padx 2 -pady 2
   grid [ttk::combobox $v.options.b.setsel -values $::MultiPlot::editdatasetlist -width 6 -state readonly -justify left -textvariable ::MultiPlot::editdataset] -row 0 -column 2 -columnspan 1 -sticky ew -padx 2 -pady 2
   bind .editmultiplot.options.b.setsel <<ComboboxSelected>> {
      ::MultiPlot::editSelect
      return
   }

   grid [label $v.options.b.llbl -text "Lines" -width 10 -anchor w] -row 1 -column 0 -sticky news -padx 2 -pady 2
   grid [spinbox $v.options.b.lsize -textvariable ::MultiPlot::editlinewidth -width 4 -from 1 -to 10 -increment 1 -command {::MultiPlot::editLine width}] -row 1 -column 1 -sticky news -padx 1 -pady 2
   set linetype {"full" "dash -" "dash ," "dash ." "dash _"}
   grid [ttk::combobox $v.options.b.combline1 -values $linetype -width 6 -state readonly -justify left -textvariable ::MultiPlot::editlinedash] -row 1 -column 2 -pady 0 -sticky ns -padx $px
   bind .editmultiplot.options.b.combline1 <<ComboboxSelected>> {
      ::MultiPlot::editLine dash
      return
   }
   grid [button $v.options.b.lcolor -text "Color" -command {
      set auxcolor [tk_chooseColor -initialcolor $::MultiPlot::editlinecolor -title "Choose line color"]
      if {$auxcolor != ""} {set ::MultiPlot::editlinecolor $auxcolor ;::MultiPlot::editLine color}
      }] -row 1 -column 3 -sticky news

   grid [label $v.options.b.mlbl -text "Marker" -width 10 -anchor w] -row 2 -column 0 -sticky news -padx 2 -pady 2
   grid [spinbox $v.options.b.msize -textvariable ::MultiPlot::editmarkersize -width 4 -from 1 -to 10 -increment 1 -command {::MultiPlot::editMarker size}] -row 2 -column 1 -sticky news -padx 1 -pady 2
   set marktype {"none" "point" "circle" "square"}
   grid [ttk::combobox $v.options.b.combmark1 -values $marktype -width 6 -state readonly -justify left -textvariable ::MultiPlot::editmarkertype] -row 2 -column 2 -pady 0 -sticky ns -padx $px
   bind .editmultiplot.options.b.combmark1 <<ComboboxSelected>> {
      ::MultiPlot::editMarker type
      return
   }
   grid [button $v.options.b.mcolor -text "Color" -command {
      set auxcolor [tk_chooseColor -initialcolor $::MultiPlot::editmarkercolor -title "Choose marker color"]
      if {$auxcolor != ""} {set ::MultiPlot::editmarkercolor $auxcolor ;::MultiPlot::editMarker color}
      }] -row 2 -column 3 -sticky news
   
   ::MultiPlot::editSelect

   ::MultiPlot::resizeGUI $v

}

proc ::MultiPlot::editLine {args} {
   # this proc tranlastes the selected line type to the -dash configure option of plothandle
   variable editlinewidth
   variable editlinecolor
   variable editlinedash
   variable editselected
   variable editdataset

   set type [lindex $args 0]
   puts "editLine: $type, editlinedash: $editlinedash"
   switch $type {
      width {
         if {[string is double $editlinewidth]} {
            ${::MultiPlot::editselected}::plothandle configure -set $editdataset -linewidth $editlinewidth -plot
         } else {
            puts "Not a number" ;# ADD ERROR MESSAGE
         }
      }
      color {
         set ::MultiPlot::editlinecolor [lindex $::MultiPlot::editlinecolor 0]
         ${::MultiPlot::editselected}::plothandle configure -set $editdataset -linecolor $editlinecolor -plot
         #switch $editlinecolor { black   {[set ::MultiPlot::editselected]::plothandle configure -linecolor "black" -plot} }
      }
      dash {
         switch $editlinedash {
            "full"   {[set ::MultiPlot::editselected]::plothandle configure -set $editdataset -dash "" -plot}
            "dash -" {[set ::MultiPlot::editselected]::plothandle configure -set $editdataset -dash "-" -plot}
            "dash ," {[set ::MultiPlot::editselected]::plothandle configure -set $editdataset -dash "," -plot} 
            "dash ." {[set ::MultiPlot::editselected]::plothandle configure -set $editdataset -dash "." -plot}
            "dash _" {[set ::MultiPlot::editselected]::plothandle configure -set $editdataset -dash "_" -plot}
         }
      }
   }
}

proc ::MultiPlot::editMarker {args} {
   variable editmarkercolor
   variable editmarkersize
   variable editmarkertype
   variable editselected
   variable editdataset

   set aux [lindex $args 0]
   puts "editMarker: $aux"
   switch $aux {
      size {
         if {[string is double $editmarkersize]} {
            ${::MultiPlot::editselected}::plothandle configure -set $editdataset -radius $editmarkersize -plot
         } else {
            puts "Not a number" ;# ADD ERROR MESSAGE
         }
      }
      color {
         set ::MultiPlot::editmarkercolor [lindex $::MultiPlot::editmarkercolor 0]
         ${::MultiPlot::editselected}::plothandle configure -set $editdataset -fillcolor $editmarkercolor -plot
      }
      type {
         if {$editmarkertype == "none" || $editmarkertype == "point" || $editmarkertype == "circle" || $editmarkertype == "square"} {
            ${::MultiPlot::editselected}::plothandle configure -set $editdataset -marker $editmarkertype -plot
         } else {
            puts "Not a marker type" ;# ADD ERROR MESSAGE
         }
      }
   }

}

proc ::MultiPlot::editPeriphery {args} {
   variable editplotysize
   variable editplotxsize
   variable edittitle
   variable edittitlefont
   variable edittitlefontsize
   variable editylabel
   variable editylabelfont
   variable editylabelfontsize
   variable editxlabel
   variable editxlabelfont
   variable editxlabelfontsize
   variable editymin
   variable editymax
   variable editxmin
   variable editxmax
   variable editselected
   
   set aux [lindex $args 0]
   switch $aux {
      size   {${editselected}::plothandle configure -ysize $editplotysize -xsize $editplotxsize -plot}
      title  {${editselected}::plothandle configure -title $edittitle -titlefont $edittitlefont -titlefontsize $edittitlefontsize -plot}
      ylabel {${editselected}::plothandle configure -ylabel $editylabel -ylabelfont $editylabelfont -ylabelfontsize $editylabelfontsize -plot}
      xlabel {${editselected}::plothandle configure -xlabel $editxlabel -xlabelfont $editxlabelfont -xlabelfontsize $editxlabelfontsize -plot}
      ymin   {if {[regexp {[0-9]} $editymin]} {${editselected}::plothandle configure -ymin $editymin -plot}}
      ymax   {if {[regexp {[0-9]} $editymax]} {${editselected}::plothandle configure -ymax $editymax -plot}}
      xmin   {if {[regexp {[0-9]} $editxmin]} {${editselected}::plothandle configure -xmin $editxmin -plot}}
      xmax   {if {[regexp {[0-9]} $editxmax]} {${editselected}::plothandle configure -xmax $editxmax -plot}}
      yauto  {${editselected}::plothandle configure -ymin auto -ymax auto -plot}
      xauto  {${editselected}::plothandle configure -xmin auto -xmax auto -plot}
   }

}

proc ::MultiPlot::editSelect {args} {
   variable editlist
   variable editselected
   variable editdataset
   variable editauxvalues
   variable editplotysize
   variable editplotxsize
   variable edittitle
   variable edittitlefont
   variable edittitlefontsize

   variable editlinewidth
   variable editlinecolor
   variable editlinedash

   variable editmarkersize
   variable editmarkercolor
   variable editmarkertype
   
   variable editylabel
   variable editylabelfont
   variable editylabelfontsize
   variable editxlabel
   variable editxlabelfont
   variable editxlabelfontsize

   variable editgridy
   variable editgridx
   
   variable editylabel
   variable editxlabel
   variable editymin
   variable editymax
   variable editxmin
   variable editxmax

   set ::MultiPlot::editplotysize [set ${editselected}::canh]
   set ::MultiPlot::editplotxsize [set ${editselected}::canw]

   set ::MultiPlot::edittitle [set ${editselected}::title]
   set ::MultiPlot::edittitlefont [set ${editselected}::titlefont]
   set ::MultiPlot::edittitlefontsize [set ${editselected}::titlefontsize]
   set ::MultiPlot::editylabelfont [set ${editselected}::ylabelfont]
   set ::MultiPlot::editylabelfontsize [set ${editselected}::ylabelfontsize]
   set ::MultiPlot::editxlabelfont [set ${editselected}::xlabelfont]
   set ::MultiPlot::editxlabelfontsize [set ${editselected}::xlabelfontsize]
   set ::MultiPlot::editylabel [set ${editselected}::ylabeltext]
   set ::MultiPlot::editxlabel [set ${editselected}::xlabeltext]
   set ::MultiPlot::editymin [string trim [set ${editselected}::ymin]]
   set ::MultiPlot::editymax [string trim [set ${editselected}::ymax]]
   set ::MultiPlot::editxmin [string trim [set ${editselected}::xmin]]
   set ::MultiPlot::editxmax [string trim [set ${editselected}::xmax]]
   set ::MultiPlot::editgridy [set ${editselected}::gridy]
   set ::MultiPlot::editgridx [set ${editselected}::gridx]
   
   # Variables that depend on dataset
   set ::MultiPlot::editlinewidth [lindex [set ${editselected}::linewidth] $editdataset]
   set ::MultiPlot::editlinecolor [lindex [set ${editselected}::linecolor] $editdataset]
   set auxlinedash [lindex [set ${editselected}::dashed] $editdataset]
   switch $auxlinedash {
      {} {set ::MultiPlot::editlinedash "full"}
      -  {set ::MultiPlot::editlinedash "dash -"} 
      ,  {set ::MultiPlot::editlinedash "dash ,"}
      .  {set ::MultiPlot::editlinedash "dash ."} 
      _  {set ::MultiPlot::editlinedash "dash _"} 
   }
   set ::MultiPlot::editmarkersize  [lindex [set ${editselected}::radius] $editdataset]
   set ::MultiPlot::editmarkercolor [lindex [set ${editselected}::fillcolor] $editdataset]
   set ::MultiPlot::editmarkertype  [lindex [set ${editselected}::marker] $editdataset]
   
   # update editlist
   set editlist {}
   foreach s [multiplot list] {
      lappend editlist "::[lindex [split $s ":"] 2]::[lindex [split $s ":"] 4]"
   }

   set editdatasetlist {}
   set nsets [${editselected}::plothandle nsets]
   for {set i 0} {$i < $nsets} {incr i} {lappend editdatasetlist $i}

}

proc ::MultiPlot::getVars {plot var} {
   if {[catch {set ::MultiPlot::Plot${plot}::title}]} {
      puts "Not a valid plot number"
      return
   }
   # check the args is a valid variable
   if {[catch {set ::MultiPlot::Plot${plot}::$var}]} {
      puts "Not a valid variable"
      return
   }
   return [set ::MultiPlot::Plot${plot}::$var]
}

proc ::MultiPlot::manageStyle {mode ns} {
   variable edittitle 
   variable edittitlefont 
   variable edittitlefontsize
   variable editylabel
   variable editylabelfont
   variable editylabelfontsize 
   variable editxlabel 
   variable editxlabelfont
   variable editxlabelfontsize
   variable editymin 
   variable editymax 
   variable editxmin 
   variable editxmax 
   variable editlinedash 
   variable editlinewidth 
   variable editlinecolor 
   variable editmarkersize
   variable editmarkercolor
   variable editmarkertype

   set varlist [list title titlefont titlefontsize linewidth linecolor \
      radius fillcolor marker ylabelfont ylabelfontsize xlabelfont xlabelfontsize \
      ylabeltext xlabeltext ymin ymax xmin xmax canh canw]

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
         set fid [open $newpathfile w]
         puts $fid "# MultiPlot options saved by style manager"
         foreach var $varlist {
               puts $fid "set ${ns}::${var} \"[::MultiPlot::getVars [string range $ns end end] $var]\""
         }
         puts $fid "${ns}::plothandle replot"
         close $fid  
      }
      load {
         set newpathfile [tk_getOpenFile \
			  -title "Choose file name" \
			  -initialdir [pwd] -filetypes $types]
         set fid [open $newpathfile r]
         #check this is a multiplot options file
         set line [gets $fid]
         if {[regexp {# MultiPlot} $line] == 1} {
               source $newpathfile   
         } else {puts "MultiPlot) Not a valid style file"}
      }
   }   

}

proc ::MultiPlot::editGrid {val mode} {
  # mode refers to which tics to grid can be either "major, all, none"
  variable editselected
  set a ${editselected}::plothandle
  set xmin [string trim [set ${editselected}::xmin]]
  set xmax [string trim [set ${editselected}::xmax]]
  set ymin [string trim [set ${editselected}::ymin]]
  set ymax [string trim [set ${editselected}::ymax]]

  set hiticx [set ${editselected}::hiticx]
  set loticx [set ${editselected}::loticx]
  set dimx [set ${editselected}::dimx]
  set minorticx [set ${editselected}::minorticx]
  set hiticy [set ${editselected}::hiticy]
  set loticy [set ${editselected}::loticy]
  set dimy [set ${editselected}::dimy]
  set minorticy [set ${editselected}::minorticy]

  set i 1
  
  # Lets try the following idea:
  # first define the x values to draw based on the tick mode
  switch $val {
    "x" {
      set ticval [expr $loticx+$i*$dimx/$minorticx]
      set firstmajor [expr abs(int($loticx-int($loticx/$minorticx)*$minorticx))]
      set xtogrid {}
      if {$mode == "all"} {
        # add minor tics
        while {$ticval<$hiticx} {
          lappend xtogrid $ticval       
          incr i
          set ticval [expr $loticx+$i*$dimx/$minorticx]
        }
        set ${editselected}::gridx all
      } elseif {$mode == "major"} {
        while {$ticval<$hiticx} {
          if {![expr ($i-$firstmajor)%$minorticx]} {
            # This means I am at a major tic
            lappend xtogrid $ticval
          }       
          incr i
          set ticval [expr $loticx+$i*$dimx/$minorticx]
        }
        set ${editselected}::gridx major
      } elseif {$mode == "none"} {
         set ${editselected}::gridx none
      }
      # After populating xtogrid, draw the lines
      $a undraw gridx
      foreach x $xtogrid {
        $a draw line $x $ymin $x $ymax -dash - -tag gridx
      }
    }
    "y" {
      set ticval [expr $loticy+$i*$dimy/$minorticy]
      set firstmajor [expr abs(int($loticy-int($loticy/$minorticy)*$minorticy))]
      set ytogrid {}
      if {$mode == "all"} {
        # add minor tics
        while {$ticval<$hiticy} {
          lappend ytogrid $ticval       
          incr i
          set ticval [expr $loticy+$i*$dimy/$minorticy]
        }
        set ${editselected}::gridy all
      } elseif {$mode == "major"} {
        while {$ticval<$hiticy} {
          if {![expr ($i-$firstmajor)%$minorticy]} {
            # This means I am at a major tic
            lappend ytogrid $ticval
          }       
          incr i
          set ticval [expr $loticy+$i*$dimy/$minorticy]
        }
        set ${editselected}::gridy major
      } elseif {$mode == "none"} {
         set ${editselected}::gridy none   
      }
      $a undraw gridy
      # After populating xtogrid, draw the lines
      foreach y $ytogrid {
        $a draw line $xmin $y $xmax $y -dash - -tag gridy
      }
    }
  }
}
