.SILENT:

VMFILES = CHANGES.txt COPYRIGHT.txt README.txt pkgIndex.tcl \
	tablelistPublic.tcl tablelist.tcl tablelist_tile.tcl \
	tablelistCommon.tcl

VMFILES2 = scripts/mwutil.tcl scripts/pencil.cur scripts/repair.tcl \
	scripts/tablelistBind.tcl scripts/tablelistConfig.tcl \
	scripts/tablelistEdit.tcl scripts/tablelistImages.tcl \
	scripts/tablelistMove.tcl scripts/tablelistSort.tcl \
	scripts/tablelistThemes.tcl scripts/tablelistUtil.tcl \
	scripts/tablelistWidget.tcl scripts/tclIndex

VMFILES3 = scripts/utils/mwutil.tcl scripts/utils/pkgIndex.tcl \
	scripts/utils/scaleutilMisc.tcl scripts/utils/scaleutil.tcl 

VMVERSION = 6.11
DIR = $(PLUGINDIR)/noarch/tcl/tablelist$(VMVERSION)

bins:
win32bins:
dynlibs:
staticlibs:
win32staticlibs:

distrib:
	@echo "Copying tablelist $(VMVERSION) files to $(DIR)"
	mkdir -p $(DIR)/scripts
	mkdir -p $(DIR)/scripts/utils
	cp $(VMFILES) $(DIR) 
	cp $(VMFILES2) $(DIR)/scripts
	cp $(VMFILES3) $(DIR)/scripts/utils

	
