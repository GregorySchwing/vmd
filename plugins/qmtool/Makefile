.SILENT:

VMFILES = qmtool.tcl qmtool_analysis.tcl qmtool_atomedit.tcl \
          qmtool_aux.tcl qmtool_charges.tcl qmtool_intcoor.tcl \
          qmtool_readwrite.tcl qmtool_setup.tcl qmtool_gui.tcl \
	  qmtool_gaussian.tcl qmtool_orca.tcl qmtool_gamess.tcl \
	  pkgIndex.tcl

VMVERSION = 2.0
DIR = $(PLUGINDIR)/noarch/tcl/qmtool$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying qmtool $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp $(VMFILES) $(DIR) 

	
