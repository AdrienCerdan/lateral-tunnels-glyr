# Permeation pathway for ions in GlyR
# Adrien Cerdan, Univerty of strasbourg, France



# Track the pathway and the direction of all Chloride.

# Usage:
#vmd -dispdev text prot_ions.psf -e pathway-v5.tcl -args ${prodlist[*]}
# or
#vmd -dispdev text prot_ions.psf -e pathway-v5.tcl -args prod01.dcd



#exec rm comz.dat
exec rm portal-trans.dat
exec rm tmp

source bigdcd-skip.tcl



#applied at each frame
proc COMZ {frame} {
	#puts "hello"
	#global molID molList atomID1 atomID2 framerate
	global molID molList ringIN1 ringOUT1 ringOUT1BOT ringIN2 ringOUT2 ringOUT2BOT ringIN3 ringOUT3 ringOUT3BOT ringIN4 ringOUT4 ringOUT4BOT ringIN5 ringOUT5 ringOUT5BOT framerate flagOUT flagIN OUT IN flag radius listCLA Dflag Din Dout Dflagout Dflagin bulk vestibule centerRING transition
	#puts "bye"
	#get postion of GM, and position of the CL-
	set posIN1 [measure center $ringIN1]
	set posIN2 [measure center $ringIN2]
	set posIN3 [measure center $ringIN3]
	set posIN4 [measure center $ringIN4]
	set posIN5 [measure center $ringIN5]
	set posOUT1 [measure center $ringOUT1]
	set posOUT2 [measure center $ringOUT2]
	set posOUT3 [measure center $ringOUT3]
	set posOUT4 [measure center $ringOUT4]
	set posOUT5 [measure center $ringOUT5]
	set posOUT1BOT [measure center $ringOUT1BOT]
	set posOUT2BOT [measure center $ringOUT2BOT]
	set posOUT3BOT [measure center $ringOUT3BOT]
	set posOUT4BOT [measure center $ringOUT4BOT]
	set posOUT5BOT [measure center $ringOUT5BOT]

	set posCENTER [measure center $centerRING]

#loop over all selected molecules
	foreach mol $listCLA {
	#selections
		set CLA [atomselect 0 "index $mol"]

	
  	#set posCLA [lindex [$CLA get {x y z}] 0]
		set posCLA [measure center $CLA]
	
	
		# distance with the ringOUT
		set distanceOUT1 [veclength [vecsub $posOUT1 $posCLA]]
		set distanceOUT2 [veclength [vecsub $posOUT2 $posCLA]]
		set distanceOUT3 [veclength [vecsub $posOUT3 $posCLA]]
		set distanceOUT4 [veclength [vecsub $posOUT4 $posCLA]]
		set distanceOUT5 [veclength [vecsub $posOUT5 $posCLA]]

		# distance with the ringOUT
		set distanceOUT1BOT [veclength [vecsub $posOUT1BOT $posCLA]]
		set distanceOUT2BOT [veclength [vecsub $posOUT2BOT $posCLA]]
		set distanceOUT3BOT [veclength [vecsub $posOUT3BOT $posCLA]]
		set distanceOUT4BOT [veclength [vecsub $posOUT4BOT $posCLA]]
		set distanceOUT5BOT [veclength [vecsub $posOUT5BOT $posCLA]]

		# distance with the ringIN
		set distanceIN1 [veclength [vecsub $posIN1 $posCLA]]
		set distanceIN2 [veclength [vecsub $posIN2 $posCLA]]
		set distanceIN3 [veclength [vecsub $posIN3 $posCLA]]
		set distanceIN4 [veclength [vecsub $posIN4 $posCLA]]
		set distanceIN5 [veclength [vecsub $posIN5 $posCLA]]

		# distance bulk vs vestibule
		set distanceCENTER [veclength [vecsub $posCENTER $posCLA]]

		
		# [OPTI] we could continue all ion either compartement and track only the ones in the middle.
		if { [expr {$distanceCENTER < 10}] } {
			if { [expr {$frame != 1 } ] && [expr {[dict get $vestibule $mol] == 0 } ] } {
				dict set transition $mol 1
				dict set vestibule $mol 1
				dict set bulk $mol 0
				#puts "transition $mol , vestibule = [dict get $vestibule $mol] , frame = $frame " 
			} else {
				dict set transition $mol 0
				dict set vestibule $mol 1
				dict set bulk $mol 0
				continue 
			} 
		} elseif { [expr {$distanceCENTER > 40}] } {
			if { [expr {$frame != 1 } ] && [expr {[dict get $bulk $mol] == 0 } ] } {
				dict set transition $mol 1
				dict set vestibule $mol 0
				dict set bulk $mol 1
				#puts "transition $mol , vestibule = [dict get $vestibule $mol] , frame = $frame "
			} else {
				dict set transition $mol 0
				dict set vestibule $mol 0
				dict set bulk $mol 1
				continue
			}
		} else {
			dict set transition $mol 0
		}

## Record the OUT entry
		if { [expr {$distanceOUT1 < $radius}] } {
			if { [expr {$distanceOUT1 < $distanceIN1}] } {
				#set flagOUT 1
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 1
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT2 < $radius}] } {
			if { [expr {$distanceOUT2 < $distanceIN2}] } {
				#set flagOUT 2
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 2
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT3 < $radius}] } {
			if { [expr {$distanceOUT3 < $distanceIN3}] } {
				#set flagOUT 3
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 3
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT4 < $radius}] } {
			if { [expr {$distanceOUT4 < $distanceIN4}] } {
				#set flagOUT 4
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 4
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT5 < $radius}] } {
			if { [expr {$distanceOUT5 < $distanceIN5}] } {
				#set flagOUT 5
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 5
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT1BOT < $radius}] } {
			if { [expr {$distanceOUT1BOT < $distanceIN1}] } {
				#set flagOUT 1
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 1
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT2BOT < $radius}] } {
			if { [expr {$distanceOUT2BOT < $distanceIN2}] } {
				#set flagOUT 2
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 2
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT3BOT < $radius}] } {
			if { [expr {$distanceOUT3BOT < $distanceIN3}] } {
				#set flagOUT 3
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 3
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT4BOT < $radius}] } {
			if { [expr {$distanceOUT4BOT < $distanceIN4}] } {
				#set flagOUT 4
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 4
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceOUT5BOT < $radius}] } {
			if { [expr {$distanceOUT5BOT < $distanceIN5}] } {
				#set flagOUT 5
				#set flag 1
				#set OUT $frame
				dict set Dflag $mol 0
				dict set Dflagout $mol 5
				dict set Dout $mol $frame
				#puts "frame = [dict get $Dout $mol]"
				#puts "OUT = [dict get $Dflagout $mol]"
			}
		} elseif { [expr {$distanceIN1 < $radius}] } {
			if { [expr {$distanceIN1 < $distanceOUT1}] } {
				#set flagIN 1
				#set flag 1
				#set IN $frame
				dict set Dflag $mol 0
				dict set Dflagin $mol 1
				dict set Din $mol $frame
				#puts "frame = [dict get $Din $mol]"
				#puts "IN = [dict get $Dflagin $mol]"
			}
		} elseif { [expr {$distanceIN2 < $radius}] } {
			if { [expr {$distanceIN2 < $distanceOUT2}] } {
				#set flagIN 2
				#set flag 1
				#set IN $frame
				dict set Dflag $mol 0
				dict set Dflagin $mol 2
				dict set Din $mol $frame
				#puts "frame = [dict get $Din $mol]"
				#puts "IN = [dict get $Dflagin $mol]"
			}
		} elseif { [expr {$distanceIN3 < $radius}] } {
			if { [expr {$distanceIN3 < $distanceOUT3}] } {
				#set flagIN 3
				#set flag 1
				#set IN $frame
				dict set Dflag $mol 0
				dict set Dflagin $mol 3
				dict set Din $mol $frame
				#puts "frame = [dict get $Din $mol]"
				#puts "IN = [dict get $Dflagin $mol]"
			}
		} elseif { [expr {$distanceIN4 < $radius}] } {
			if { [expr {$distanceIN4 < $distanceOUT4}] } {
				#set flagIN 4
				#set flag 1
				#set IN $frame
				dict set Dflag $mol 0
				dict set Dflagin $mol 4
				dict set Din $mol $frame
				#puts "frame = [dict get $Din $mol]"
				#puts "IN = [dict get $Dflagin $mol]"
			}
		} elseif { [expr {$distanceIN5 < $radius}] } {
			if { [expr {$distanceIN5 < $distanceOUT5}] } {
				#set flagIN 5
				#set flag 1
				#set IN $frame
				dict set Dflag $mol 0
				dict set Dflagin $mol 5
				dict set Din $mol $frame
				#puts "frame = [dict get $Din $mol]"
				#puts "IN = [dict get $Dflagin $mol]"
			}
		} else {
			dict set Dflag $mol [expr {[dict get $Dflag $mol] + 1 } ]
		}




		# test is permeation happend [NEW]
		if { [expr {[dict get $transition $mol] == 1 } ] && [expr {[dict get $Dflagin $mol] == [dict get $Dflagout $mol]} ] && [expr {[dict get $Dflagin $mol] != 0} ] } {
			puts "permeation: ion = $mol in gate = [dict get $Dflagin $mol] at frame = $frame"
			dict set permeation $mol 0
			set interface [dict get $Dflagin $mol]
			dict set Dflagin $mol 0
			dict set Dflagout $mol 0
			if { [expr {[dict get $bulk $mol] == 0 } ] } {
				puts "outside -> inside"
				set direction "inward"	
			} else {
				puts "inside -> outside"
				set direction "outward"			
			}
			exec echo "$mol,$frame,$interface,$direction" >> portal-trans.dat
		}


		if { [expr {[dict get $Dflag $mol] > 200 } ] } {
			dict set Dflagin $mol 0
			dict set Dflagout $mol 0		
		}
		#set nf [expr {$frame * $framerate}]
	}
	return
}

# General Parameters
set framerate 1

set selectCLA [atomselect 0 "name CLA"]
# can be used for specific ion only:
#set selectCLA [atomselect 0 "name CLA and resid 248"]

set listCLA [$selectCLA list]

# to initiate all the dicts
foreach mol $listCLA {
	dict set Dflag $mol 0
	dict set Dflagin $mol 0
	dict set Dflagout $mol 0
	dict set Din $mol 0
	dict set Dout $mol 0
	dict set bulk $mol 0
	dict set vestibule $mol 0
	dict set transition $mol 0
	
}
puts "number of tracked ions: [dict size $Dflag]"


# set ring atoms
set ringIN1 [atomselect 0 "protein and ((segname PROA and resid 120 71) or (segname PROB and resid 66 75))"]
set ringIN2 [atomselect 0 "protein and ((segname PROB and resid 120 71) or (segname PROC and resid 66 75))"]
set ringIN3 [atomselect 0 "protein and ((segname PROC and resid 120 71) or (segname PROD and resid 66 75))"]
set ringIN4 [atomselect 0 "protein and ((segname PROD and resid 120 71) or (segname PROE and resid 66 75))"]
set ringIN5 [atomselect 0 "protein and ((segname PROE and resid 120 71) or (segname PROA and resid 66 75))"]

#OLD RINGout
#set ringOUT1 [atomselect 0 "protein and ((segname PROA and resid 153 155 158 291) or (segname PROB and resid 198 199 202))"]
#set ringOUT2 [atomselect 0 "protein and ((segname PROB and resid 153 155 158 291) or (segname PROC and resid 198 199 202))"]
#set ringOUT3 [atomselect 0 "protein and ((segname PROC and resid 153 155 158 291) or (segname PROD and resid 198 199 202))"]
#set ringOUT4 [atomselect 0 "protein and ((segname PROD and resid 153 155 158 291) or (segname PROE and resid 198 199 202))"]
#set ringOUT5 [atomselect 0 "protein and ((segname PROE and resid 153 155 158 291) or (segname PROA and resid 198 199 202))"]

#new RINGout1
#set ringOUT1 [atomselect 0 "protein and ((segname PROA and resid 153 155 118) or (segname PROB and resid 198 199 63))"]
#set ringOUT2 [atomselect 0 "protein and ((segname PROB and resid 153 155 118) or (segname PROC and resid 198 199 63))"]
#set ringOUT3 [atomselect 0 "protein and ((segname PROC and resid 153 155 118) or (segname PROD and resid 198 199 63))"]
#set ringOUT4 [atomselect 0 "protein and ((segname PROD and resid 153 155 118) or (segname PROE and resid 198 199 63))"]
#set ringOUT5 [atomselect 0 "protein and ((segname PROE and resid 153 155 118) or (segname PROA and resid 198 199 63))"]

#new RINGout2TOP
set ringOUT1 [atomselect 0 "protein and ((segname PROA and resid 216 229) or (segname PROB and resid 196))"]
set ringOUT2 [atomselect 0 "protein and ((segname PROB and resid 216 229) or (segname PROC and resid 196))"]
set ringOUT3 [atomselect 0 "protein and ((segname PROC and resid 216 229) or (segname PROD and resid 196))"]
set ringOUT4 [atomselect 0 "protein and ((segname PROD and resid 216 229) or (segname PROE and resid 196))"]
set ringOUT5 [atomselect 0 "protein and ((segname PROE and resid 216 229) or (segname PROA and resid 196))"]

#new RINGout2BOT
set ringOUT1BOT [atomselect 0 "protein and ((segname PROA and resid 155) or (segname PROB and resid 199))"]
set ringOUT2BOT [atomselect 0 "protein and ((segname PROB and resid 155) or (segname PROC and resid 199))"]
set ringOUT3BOT [atomselect 0 "protein and ((segname PROC and resid 155) or (segname PROD and resid 199))"]
set ringOUT4BOT [atomselect 0 "protein and ((segname PROD and resid 155) or (segname PROE and resid 199))"]
set ringOUT5BOT [atomselect 0 "protein and ((segname PROE and resid 155) or (segname PROA and resid 199))"]

# set ring atoms
set centerRING [atomselect 0 "protein and (segname PROA PROB PROC PROD PROE and resid 120)"]

set flagOUT 0
set flagIN 0

set IN 0
set OUT 0

set flag 0

# Most important parameter to change:
set radius 6




# Do the job

# only every $framerate frames
eval bigdcd COMZ $framerate $argv 

bigdcd_wait

exit
