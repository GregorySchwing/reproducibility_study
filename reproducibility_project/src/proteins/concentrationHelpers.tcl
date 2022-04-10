
proc get_total_charge {{molid $solute}} {
	eval "vecadd [[atomselect $molid all] get charge]"
}

proc getNumIons {psfSolvated pdbSolvated {cation SOD} {cationCharge 1} {anion CLA} {anionCharge -1} {saltConcentration 0.110} {mode "sc"}} {

  # load solute molecule
  set solute [mol new $psfSolvated waitfor all]
  mol addfile $pdbSolvated mol $solute waitfor all

  # Get pbc info for later
  set xdim [molinfo $solute get a]
  set ydim [molinfo $solute get b]
  set zdim [molinfo $solute get c]

  # Compute net charge of the system
  set sel [atomselect $solute all]
  set netCharge [eval "vecadd [$sel get charge]"]
  set roundNetCharge [expr round($netCharge)]
  $sel delete

  set done 0 ;# flag to tell if there are ions to be placed
  set nCation 0
  set nAnion 0

  set nIonsList {}
  # For each ion placement mode, calculate the number of each ion
  # and set the nIonsList accordingly.

  ###
  ### Ion placement mode 'neutralize'. Also called in combination 
  ### with other modes.
  ###
    # XXX - The following implementation will work only for ions with 
    #       charge -2, -1, +1, and +2. 
  if {$mode != "nions"} {
    set errflag 0
    if {$roundNetCharge > 0} {
      if {$anionCharge == -1} {
        set nAnion $roundNetCharge
        set nCation 0
      } elseif {$anionCharge == -2} {
        if {[expr {$roundNetCharge % 2}] == 0} { ;# total charge is even
          set nAnion [expr {$roundNetCharge / 2}]
          set nCation 0
        } else { ;# total charge is odd
          if {$cationCharge == 1} {
            set nAnion [expr {1+int($roundNetCharge/2)}]
            set nCation 1
          } else {
            set errflag 1
          }
        }
      }
    } elseif {$roundNetCharge < 0} {
      set roundNetChargeAbs [expr {-1*$roundNetCharge}]
      if {$cationCharge == 1} {
        set nAnion 0
        set nCation $roundNetChargeAbs
      } elseif {$cationCharge == 2} {
        if {[expr {$roundNetChargeAbs % 2}] == 0} { ;# total charge is even
          set nAnion 0
          set nCation [expr {$roundNetChargeAbs / 2}]
        } else { ;# total charge is odd
          if {$anionCharge == -1} {
            set nAnion 1
            set nCation [expr {1+int($roundNetChargeAbs/2)}]
          } else {
            set errflag 1
          }
        }
      }
    }
  }

  ###
  ### Ion placement modes 'sc' = salt concentration
  ###                     'is' = ionic strength
  ###
  if {$mode == "sc"} {
    puts "Autoionize) Desired salt concentration: ${saltConcentration} mol/L."
  } elseif {$mode == "is"} {
    puts "Autoionize) Desired ionic strength: ${ionicStrength} mol/L."
  }
  
  if {$mode == "sc" || $mode == "is"} {

    set sel [atomselect $solute "water and noh"]
    set nWater [$sel num]
    $sel delete

    if {$nWater == 0} {
      error "Autoionize) ERROR: Cannot add ions to unsolvated system."
    }

    # Guess chemical formula. 
    # XXX - We currently only support ions with charge -2, -1, +1, and +2.
    if {$cationCharge == 1 && $anionCharge == -1} { ;# e.g., NaCl, KCl, ...
      set cationStoich 1
      set anionStoich 1
    } elseif {$cationCharge == 2 && $anionCharge == -1} { ;# e.g., MgCl2
      set cationStoich 1
      set anionStoich 2
    } elseif {$cationCharge == 1 && $anionCharge == -2} {
      set cationStoich 2
      set anionStoich 1
    } elseif {$cationCharge == 2 && $anionCharge == -2} {
      set cationStoich 1
      set anionStoich 1
    } else {
      error "Autoionize) ERROR: Unsupported ion charge; cannot guess chemical formula."
    }

    if {$mode == "is"} { ;# convert ionic strength to salt concentration
      set cationConcentration [expr {2 * $ionicStrength / ( sqrt($cationCharge * $cationCharge * $anionCharge * $anionCharge) + $cationCharge * $cationCharge)}]
      set saltConcentration [expr {$cationConcentration * $cationStoich}]
    }

    # As long as saltConcentration and ionicStrength are non-negative,
    # no error checking is needed here...
    set num [expr {int(0.5 + 0.0187 * $saltConcentration * $nWater)}]
    set nCation [expr {$nCation + $cationStoich * $num}]
    set nAnion [expr {$nAnion + $anionStoich * $num}]

  }

  if {$mode != "nions"} {
    if {$nCation > 0} {
      lappend nIonsList [list $cation $cationCharge $nCation]
    } 
    if {$nAnion > 0} {
      lappend nIonsList [list $anion $anionCharge $nAnion]
    }
    if {$nCation == 0 && $nAnion == 0} {
      # Just in case the system is neutral and the requested ionic
      # strength or salt concentration is zero, in which case the
      # behavior is the same as -neutralize, i.e., copy files and exit
      # normally.
      set done 1
    }
    if {$nCation < 0 || $nAnion < 0} {
      error "Autoionize) ERROR: Internal error; negative number of ions."
    }
  }

  return [list $nCation $nAnion]
}


