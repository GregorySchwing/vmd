.SILENT:

VMFILES = pkgIndex.tcl qwikmd_ballon.tcl qwikmd_func.tcl qwikmd_info.tcl qwikmd.tcl qwikmd_logText.tcl arrow_right.gif arrow_down.gif Annealing.conf Equilibration.conf MD.conf Minimization.conf SMD.conf QMMM-Min.conf QMMM-Ann.conf QMMM-Equi.conf QMMM.conf

VMVERSION = 1.3
DIR = $(PLUGINDIR)/noarch/tcl/qwikmd$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying qwikmd $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp $(VMFILES) $(DIR) 

	
