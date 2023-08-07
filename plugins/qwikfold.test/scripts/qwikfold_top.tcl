#~ proc QWIKFOLD::top_frame {} {
#~ variable FoldingEngine
#~ variable OutputModels
#~ variable DOWNLOAD_DIR
########################################################################
# Top Frame
########################################################################
grid [ ttk::frame $QWIKFOLD::topGui.top -borderwidth 10 ] -column 0  -row 0

# Folding Engine -------------------------------------------------------
grid [ ttk::label 	$QWIKFOLD::topGui.top.label0 -text "Folding Engine:" ] -column 0 -row 1 -padx 5 -sticky w
tk_optionMenu 		$QWIKFOLD::topGui.top.value0 FoldingEngine "AlphaFold 2" "RoseTTAfold" "GSAFold"
grid 				$QWIKFOLD::topGui.top.value0 -column 0 -row 2 -padx 5  -sticky w

#grid [ ttk::label	$QWIKFOLD::topGui.top.label1 -text "Models to create:" ]
#tk_optionMenu 		$QWIKFOLD::topGui.top.value1 OutputModels "1" "2" "3" "4" "5"
#grid 							$QWIKFOLD::topGui.top.label1 -column 1 -row 1
#grid 							$QWIKFOLD::topGui.top.value1 -column 1 -row 2

grid [ ttk::label		$QWIKFOLD::topGui.top.label1 -text "AlphaFold Models" ] -column 1 -columnspan 2 -row 1 
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check1 -text "model_1" -variable model_1 ]   -column 1 -row 2 -padx 2 -sticky w
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check2 -text "model_2" -variable model_2 ]   -column 1 -row 3 -padx 2 -sticky w
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check3 -text "model_3" -variable model_3 ]   -column 1 -row 4 -padx 2 -sticky w
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check4 -text "model_4" -variable model_4 ]   -column 2 -row 2 -padx 2 -sticky w
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check5 -text "model_5" -variable model_5 ]   -column 2 -row 3 -padx 2 -sticky w
grid [ ttk::checkbutton $QWIKFOLD::topGui.top.check6 -text "all"     -variable model_all ] -column 2 -row 4 -padx 2 -sticky w
#~ }
