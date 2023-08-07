.SILENT:

VMFILES := $(wildcard *.tcl ) 
VMVERSION = 0.5
DIR = $(PLUGINDIR)/noarch/tcl/vnd$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying vnd $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp -r $(VMFILES) $(DIR) 

