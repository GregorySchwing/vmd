.SILENT:

VMFILES := $(wildcard *.tcl *.py ) examples patches python
VMVERSION = 0.7
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

