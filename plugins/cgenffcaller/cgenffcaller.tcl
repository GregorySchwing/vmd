#
# $Id: cgenffcaller.tcl,v 1.2 2020/02/04 19:25:33 jribeiro Exp $
#
#==============================================================================
# CgenFFCaller plugin is intended to be an interface between VMD and the 
# atom typing and CHARMM force field toolkit webserver CgenFF 
# (https://cgenff.umaryland.edu/). The user is required to register at the 
# the website, and obtain a username and password. 
#
# The plugin always expects a mol2 or pdb file, username. 
# 
# NOTE: The password is not a CGenFF requirement. CGenFFCALLER is keeping
#       the use of the password in case they ever change this. The GUI for
#       submitting the password will be hidden, but the hooks are here.
#
# TODO: finish the return documentation
# 
# Authors:
#   João V. Ribeiro
#   Beckman Institute for Advanced Science and Technology
#   University of Illinois, Urbana-Champaign
#   jribeiro@ks.uiuc.edu
#   http://www.ks.uiuc.edu/~jribeiro
#=============================================================================


package provide cgenffcaller 1.0
package require exectool
package require http
package require json
package require platform

namespace eval ::CGENFFCALLER:: {

  set streamfiles [list]; # the stream files returned by the server
  set totstreamfiles 0; # number of stream files stored
  set count -1; # the number of submissions triggered by the user with 
                # the same account
}

### Global call for the the ::CGENFFCALLER::callCGENFF

proc cgenffcaller { args } {

  if {[llength $args ] != 0} {
    ### Global call of the main proc
    if {[lindex $args 0] == "callcgenff"} {
      return [eval CGENFFCALLER::callCGENFF [lrange $args 1 end]]
    } elseif {[lindex $args 0] == "savemol"} {
      ### Global call of the proc to export CGenFF compatible molecule files
      return [eval CGENFFCALLER::saveMol [lrange $args 1 end]]
    } elseif {[lindex $args 0] == "numstr"} {
      return [eval CGENFFCALLER::getNumStreamFiles]
    } elseif {[lindex $args 0] == "getstr"} {
      return [eval CGENFFCALLER::getStream [lindex $args end]]
    } elseif {[lindex $args 0] == "getsubcount"} {
      return [eval CGENFFCALLER::getCount]
    }
  }

  ::CGENFFCALLER::callCGENFF_usage
  ::CGENFFCALLER::saveMol_usage
  ::CGENFFCALLER::getNumStreamFiles_usage
  ::CGENFFCALLER::getStream_usage
  ::CGENFFCALLER::getCount_usage
  return -1

}

proc ::CGENFFCALLER::callCGENFF_usage {} {

  puts "cgenffcaller callcgenff -molfiles <List of Files> -username <CGenFF Server Username>\
          (optional) -method <Web Server Fetching Tool - default value = \"auto\">: "
  #puts "callCGENFF -molfile <List of Files> -username <CGenFF Server Username>\
          -password <CGenFF Server Password> -method <Web Server Fetching Tool - default value = \"auto\">"
  puts "Fetch the topology and parameters for the molecule files in the -molfiles\
         list, and store them in a stream file (str)."
  puts "The stream files are stored internally, and can be accessed executing the\
         command \"callcgenff getstr <stream index>\".\n"
  puts "List of Options:"
  puts "    -molfiles <List of Files> = List of files to be submitted to CGenFF Server. Mol2 and PDB formats supported"
  puts "    -username <CGenFF Server Username> = CGenFF Server Account Username"
  #puts "    -password <CGenFF Server Password> = CGenFF Server Account Password"
  puts "    -method <Web Server Fetching Tool> = Tool to use to fetch the information from the CGenFF Server."
  puts "     The tool is automatically select whenever is available. Linux and Mac use \"wget\" and \"curl\";"
  puts "     Windows machine use \"curl\" executed inside the powershell (Windows 10 only). This flag enforces" 
  puts "     the use the specified tool (for testing and debugging).\n"

  return -1
}


proc ::CGENFFCALLER::saveMol_usage {} {

  puts "cgenffcaller savemol -sel <atom selection object> -output <output file name>\
          -format <file format: mol2, xyz, or pdb>\n"
  return -1
}

proc ::CGENFFCALLER::getNumStreamFiles_usage {} {

  puts "cgenffcaller numstr : Returns the number of stream files currently stored\
         in the plugin.\n"

  return -1
}

proc ::CGENFFCALLER::getStream_usage {} {

  puts "cgenffcaller getstr <index of the stream> : Returns the stream file\
         stored currently stored in the <index>th position.\n"
  return -1
}

proc ::CGENFFCALLER::getCount_usage {} {

  puts "cgenffcaller getsubcount: Returns the number of submissions performed\
         by the user. This command only works after the one submission\
         had been performed.\n"
  return -1
}
###############################################################################
# Fetch topology file from the CGenFF server:
#   INPUT: list of files with molecule information (pdb or mol2)
#   OUTPUT: Stored in the global variables streamfiles and cout and can be
#           retrieved using the commands ::CGENFFCALLER::getStream and 
#           ::CGENFFCALLER::getCount'
#
#   The "method" defines which command will be executed to fetch info
#   from the server. This is useful to test the different mechanism. Default
#   value is auto
#
#   ::CGENFFCALLER::saveMol save a molecule as mol2, pdb and xyz and improves the
#   compatibility with the CGenFF (see comments in the proc definition)
#
###############################################################################

proc ::CGENFFCALLER::callCGENFF { args } {

  set molfiles [list]
  set username ""

  set password "auto"
  ### Whenever the CGenFF requirements change, activate the cleaning of the 
  ### password variable
  #set password ""
  set method "auto"
  set i 0
  set j 0
  foreach {i j} $args {

    if { $i == "-molfiles" } {
      set molfiles $j
      continue  
    }
    if { $i == "-username" } {
      set username $j
      continue
    }
    if { $i== "-password" } {
      set password $j
      continue
    }
    if { $i== "-method" } {
      set method $j
      continue
    }
  }

  if {[llength $molfiles] == 0} {
    puts "Please provide at least one file to be submitted to the CGenFF Server"
    ::CGENFFCALLER::callCGENFF_usage
    return -1
  }
  if {[llength $username] == 0} {
    puts "Please provide a valid CGenFF Server account username"
    ::CGENFFCALLER::callCGENFF_usage
    return -1
  }
  ### Uncomment whenever CGenFF requires password
  # if {[llength $password] == 0} {
  #   puts "Please provide a valid CGenFF Server account password"
  #   ::CGENFFCALLER::callCGENFF_usage
  #   return -1
  # }

  ### Reset global variables

  set CGENFFCALLER::streamfiles [list]
  set CGENFFCALLER::totstreamfiles 0
  set CGENFFCALLER::count -1

  ### Check if the arguments were defined correctly

  if {[llength $molfiles] == 0} {
    puts "Please provide one or more files (pdb or mol2) with molecule information. \
    Different molecules should be defined in different files. All files must \
    define a list."
    return -1
  }

  set username [string trim $username]

  if { $username == "" || [llength $username] == 0 } {
    puts "Please provide a valid username."
    return -1
  }

  if { $password == "" || [llength $password] == 0 } {
    puts "Please provide a valid password."
    return -1
  }


  set topfile -1

  foreach mol $molfiles {
    ### Set the parameters for the CGenFF server

    set filename $mol

    ### List of commands that will be used based on the platform and tools
    ### installed in the user's machine
    if {$method != "auto"} {
      set commandslist $method
      puts "DEBUG Command $commandslist"
    } else {
      set commandslist [list curl wget powershell python]
    }
    

    set command ""
    set cmdindex 0

    ### Check which commands are available to be used
    foreach cmd $commandslist {
      set command [auto_execok $cmd]
      if {[string length $command] != 0} {
        break
      }
      incr cmdindex
    }
    set https -1
    set str ""

    
    #set url "https://cgenff.umaryland.edu/rest/?username=${username}&?password:${password}&filename=$filename&all=true"

    set url "https://cgenff.umaryland.edu/rest/?username=${username}&filename=$filename&all=true"

    ### The command list will be used to avoid the use of tls package that, in
    ### in most cases, will not be available. In case nothing works, tls and 
    ### http will packages will be used.

    if {[string length $command] > 0} {
      switch [lindex $commandslist $cmdindex] {
        curl {
           if {[catch {eval exec "$command -s -X POST \"$url\" -H \
                \"Content-Type: text/xml\" --data-binary @$filename"} str]} {
               puts "Error: $str"   
           }
          set str [list [lindex $str 0]]
        }
        wget {

          if {[catch {eval exec $command "$url --post-file=$filename -qO-"} str]} {
            puts "Error: $str"
          }

          set str [list [lindex $str 0]]
        }
        powershell {
          ### On windows, the powershell and python should work.
        ### I noticed that to use the powershell, one needs to open the INTERNET
        ### EXPLORER first to set some security settings.
          if {[catch {exec $command "curl -Method POST -Uri \"$url\" \
            -ContentType: \"text/xml\" -InFile \"$filename\"  | Select-Object -Expand Content"} str]} {
            puts "Error: $str"
          }
        }
        python {
          ### In case none of the other methods work, use python.
          ### Some issues might raise because of the settings not being properly

          global env
          set pexec [open "pyurl.py" w+]

          set pystr "import os
print(\"here is the environ\")
print(os.environ)
print(\"end of environ\")
import requests
data = open(\'$filename\','rb').read()
resp = requests.post(\'$url\', data=data)

strfile = open ('pystream.str','w')
strfile.write(resp.text)
strfile.close()

"
          puts $pexec $pystr
          close $pexec

          package require platform
          set sourcefile ""
          set error ""
          set commandline "$command pyurl.py"

          if {[catch {eval "exec $commandline"} error]} {
            puts "ERROR $error"
          } else {

            if {[file exists "pystream.str"]} {
              set outfile [open pystream.str r]
              set str [read -nonewline $outfile]
              close $outfile
            }
          }
        }
      }
    } else {

      set https [catch {package require tls}]

      if {$https == 1} {
        puts "CGenFFCaller uses many mechanisms to access CGenFF server and \
        non of them worked. Please submit a mol2 or pdb file file directly to \
        the CGenFF server at https://cgenff.umaryland.edu/."
        return -2
      } elseif {$https == 0} {
        http::register https 443 tls::socket
      }
      set file [open $filename r]
      set data [read -nonewline $file]
      close $file
      set httpHandle [::http::geturl \"$url\" -query "$data"]
      set str [::http::data $httpHandle]
    }


    ### Parse the json object returned by the server
    set serveroutput [::json::json2dict $str ]
    set message ""

    ### Format: 
    ### flag: false - Connection successful. True - connection refused
    ### count: [num of submission] - Number of submissions made by the user
    ### output: <stream file> - stream file containing the topology and the 
    ###         parameters

    ### Check if connection failed
    if {[dict get $serveroutput flag] == true} {
      set mtype [dict get $serveroutput type]
      set message [dict get $serveroutput message]

      puts "ERROR. Connection to CGenFF refused. $message"
      return -2
    }
    
    set count [dict get $serveroutput count]
    if {$count != ""} {
      set CGENFFCALLER::count $count
    }

    set stream [dict get $serveroutput output]
    if {$stream == ""} {
      set stream [dict get $serveroutput message]
    }

    lappend CGENFFCALLER::streamfiles $stream
    incr CGENFFCALLER::totstreamfiles

    if {$https == 0} {
      http::cleanup $httpHandle
      http::unregister https
    }
  }

  return 0
}

###############################################################################
# Return the stream file stored at the position index. 
###############################################################################

proc ::CGENFFCALLER::getStream { {index -1} } {

  if {$index > -1 && $index < $CGENFFCALLER::totstreamfiles} {
    return [lindex $CGENFFCALLER::streamfiles $index]
  } else {
    puts "Please provide a valid index between 0 and [expr $CGENFFCALLER::totstreamfiles -1]"
    return -1
  }

}

###############################################################################
# Return the number of submissions triggered by the user with the same account
###############################################################################
proc ::CGENFFCALLER::getCount {} {
  
  return $CGENFFCALLER::count

}

###############################################################################
# Return the number of stream files stored
###############################################################################
proc ::CGENFFCALLER::getNumStreamFiles {} {
  
  return $CGENFFCALLER::totstreamfiles

}

###############################################################################
# Export molecule to file (XYZ, MOL2 or PDB), and edit to improve compatibility
# with CGenFF by setting the atoms' element as atoms' type and put the 
# line "generated by VMD" at the end of mol2 files 
# The command gets a atomselection, name for the output file and the format
# Although CGenFF cannot use XYZ, This file format is very useful to model small
# molecules, that I added in here, while keeping the same treatment as for
# mol2 and pdb
###############################################################################
proc ::CGENFFCALLER::saveMol { args } {
  set tmpsel ""
  set filename ""
  set format ""

  foreach {i j} $args {

    if { $i == "-sel" } {
      set tmpsel $j
      continue  
    }
    if { $i == "-output" } {
      set filename $j
      continue
    }
    if { $i== "-format" } {
      set format $j
      continue
    }
  }

  if {[llength $tmpsel] == 0} {
    puts "Please provide a valid atom selection object."
    ::CGENFFCALLER::saveMol_usage
    return -1
  }
  if {[llength $filename] == 0} {
    puts "Please provide a valid output file name"
    ::CGENFFCALLER::saveMol_usage
    return -1
  }
  if {[llength $format] == 0} {
    puts "Please provide a valid file format: mol2, xyz, or pdb."
    ::CGENFFCALLER::saveMol_usage
    return -1
  }

  ### Check if the molecule has the residue name defined. If not, abort
  set resname [lsort -unique [$tmpsel get resname]]
  if {[llength $resname] == 0} {
    puts "ERROR. Please defined an unique residue name."
    $tmpsel delete
    return -1
  }
  ### set the element as the atom's type. This tends to produce
  ### better results with the cgenff, as it eliminates confusion
  ### with atom types
  set typelist [$tmpsel get type]
  set elementlist [$tmpsel get element]
  $tmpsel set type $elementlist
  if {$format == "mol2"} {
    $tmpsel writemol2 $filename
  } elseif {$format == "xyz"} {
    $tmpsel writexyz $filename
  } elseif {$format == "pdb"} {
    $tmpsel writepdb $filename
  }

  $tmpsel set type $typelist

  if {[file exists $filename] == 1 && $format == "mol2"} {
    file copy -force ${filename} ${filename}_bkup
    set filein [open ${filename}_bkup r]
    set fileout [open ${filename} w+]
    while {[eof $filein] != 1} {
       set line [gets $filein]
       ### Send the "generated by VMD" line to the end of the file
       ### as a comment to ensure compatibility with CGenFF
       if {$line == "generated by VMD"} {
          puts $fileout "$resname"
          seek $filein [tell $filein]
          set lines [read $filein]
          puts $fileout $lines
          puts $fileout "\#$line"
          break
       } else {
          puts $fileout $line
       }
    }
    close $fileout
    close $filein
    file delete -force ${filename}_bkup
  }

}


