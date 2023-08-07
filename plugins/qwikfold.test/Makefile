.SILENT:

VMFILES := $(wildcard *.tcl *.py ) scripts examples 
VMVERSION = 0.5
DIR = $(PLUGINDIR)/noarch/tcl/qwikfold$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying qwikfold $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR) 
	cp -r $(VMFILES) $(DIR) 

