* Note: for now, assumes reference is always the first category
* and categories will be indexed 1, 2, ... regardless of actual numeric values

program define ipdfloat, eclass
	
	_on_colon_parse `0'
	local 0       `"`s(before)'"'
	local command `"`s(after)'"'
	
	syntax, STUDY(varname) POOLvars(varlist numeric min=2) [EPS(real 1e-7) FORWards BACKWards SAVING(string) SGLABel(string) CLEAR]
	gettoken trtvar intvars : poolvars
	local k : word count `poolvars'
	local k1 = `k' - 1
	tokenize `intvars'
	
	// QUICK AND DIRTY!!
	// remove `intvars' from `command', for later offset command
	local command_offset = subinstr("`command'", "`intvars'", "", 1)

	// QUICK AND DIRTY!!
	// labels -- note all subgroups are assumed to be 1, 2, ...  so labels may need to be re-defined before use with -ipdfloat-
	if "`sglabel'"!="" {
		tempfile sglabfile
		qui label save `sglabel' using `sglabfile'
	}
	
	// Interaction effects (contrasts)
	tempfile mvmeta
	qui mvmeta_make `command', usevars(`poolvars') by(`study') ppfix(check) hard saving(`mvmeta')
	preserve
		use `mvmeta', clear
		foreach v of local intvars {
			local vlist `vlist' y`v'
		}
		mvmeta y S, vars(`vlist') fixed
		test `vlist'
		
		tempname GammaHat VarGammaHat
		mat `GammaHat' = e(b)
		mat `VarGammaHat' = e(V)
		
		if "`saving'"!="" | "`clear'"!="" {
			rename y`trtvar' _y_ES1
			qui gen double _y_seES1 = sqrt(S`trtvar'`trtvar')
			forvalues i = 1 / `k1' {
				local ii = `i' + 1
				qui gen double _y_ES`ii' = _y_ES1 + y``i''
				qui gen double _y_seES`ii' = sqrt(S`trtvar'`trtvar' + S``i''``i'' + 2*S`trtvar'``i'')
				
				rename y``i'' _yInt_ES`ii'
				qui gen double _yInt_seES`ii' = sqrt(S``i''``i'')
			}			
			keep Trial _y_ES* _y_seES* _yInt_ES* _yInt_seES*

			// local sglab : value label `subgroup'
			// tempfile sglabfile
			// qui label save `sglab' using `sglabfile'
			
			tempfile subgroups
			qui save `subgroups'
		}
	restore

	// Treatment effects in reference subgroup
	tempvar offset offset_h
	qui gen double `offset_h' = .
	qui gen double `offset' = 0
	forvalues i = 1 / `k1' {
		qui replace `offset' = `offset' + `GammaHat'[1,`i']*``i''
	}
	qui ipdmetan, study(`study') poolvar(`trtvar') nogr : `command_offset', offset(`offset')
	tempname ThetaHat VarThetaHat
	scalar `ThetaHat' = r(eff)
	scalar `VarThetaHat' = r(se_eff)^2

	// Obtain "delta method" gradients, in order to form the final variance-covariance matrix
	tempname Grad h
	matrix define `Grad' = J(1, `k1', .)

	forvalues i = 1 / `k1' {
		scalar `h' = max(abs(`GammaHat'[1,`i']), 1)*sqrt(`eps')

		local sign "+"
		if "`backwards'"!="" local sign "-"
		else if "`forwards'"=="" scalar `h' = `h'/2
		
		qui replace `offset_h' = `offset' `sign' `h'*``i''
		tab ``i''
		summ `offset_h'
		qui ipdmetan, study(`study') poolvar(`trtvar') nogr : `command_offset', offset(`offset_h')
		
		if `"`forwards'`backwards'"'!=`""' {
			matrix `Grad'[1, `i'] = `sign'(r(eff) - `ThetaHat')/`h'
			continue
		}
		tempname eff_pos eff_neg
		scalar `eff_pos' = r(eff)

		qui replace `offset_h' = `offset' - `h'*``i''
		qui ipdmetan, study(`study') poolvar(`trtvar') nogr : `command_offset', offset(`offset_h')
		scalar `eff_neg' = r(eff)

		matrix `Grad'[1, `i'] = (`eff_pos' - `eff_neg')/(2*`h')
	}
		
	// Assemble matrices
	tempname ones Z A BetaHat VarBetaHat
	
	matrix define `ones' = J(`k', 1, 1)
	matrix define `Z' = J(`k', `k1', 0)
	matrix define `Z'[2, 1] = I(`k1')	
	matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat''

	// matrix define grad = grad2, grad3
	matrix define `A' = `ones'*`Grad' + `Z'
	matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	

	*** Construct saved dataset
	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {
		preserve
		
		// store pooled estimates
		tempfile pooled_ests
		tempname BetaTable GammaTable
		mat `GammaTable' = `GammaHat'', vecdiag(`VarGammaHat')'
		
		// TEMP!!
		local r = rowsof(`GammaTable')
		local ++r
		forvalues i=2/`r' {
			local rn `rn' r`i'
		}
		matrix rownames `GammaTable' = `rn'
		
		qui metani `GammaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
			saving(`pooled_ests', stacklabel) prefix(_yInt)

		mat `BetaTable' = `BetaHat', vecdiag(`VarBetaHat')'
		// matrix rownames `BetaTable' = `sgvals'
		qui metani `BetaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
			clearstack prefix(_y)

		rename _y_LABELS _yInt_LABELS
		tempvar sortorder
		qui gen int `sortorder' = _n		// -merge- will sort the data on _yInt_LABELS, so need to re-sort afterwards
		qui merge 1:1 _yInt_LABELS using `pooled_ests', assert(match master)
		
		/*
		cap {
			assert _merge==1 if _yInt_LABELS=="`ref'"
			assert _merge==3 if _yInt_LABELS!="`ref'"
		}
		if _rc {
			nois disp as err "Error in colnames [NEED TO IMPROVE THIS MESSAGE]"
			exit _rc
		}
		*/
		drop _merge *_EFFECT _yInt_USE _yInt_STUDY
		rename _y_STUDY _STUDY
		rename _yInt_LABELS _LABELS		
		qui replace _y_USE = 5 if inlist(_y_USE, 1, 2)
		sort `sortorder'
		qui save `pooled_ests', replace
		
		
		// merge subgroup estimates and interaction estimates into a single dataset
		qui use `subgroups', clear
		qui reshape long _y_ES _y_seES _yInt_ES _yInt_seES, i(Trial) j(Subgroup)
		
		if "`sglabel'"!="" {
			qui do `sglabfile'
			label values Subgroup `sglabel'
		}
		
		qui metan _y_ES _y_seES, nograph nohet nosubgroup nooverall keeporder ///
			study(Subgroup) by(Trial) rcols(_yInt_ES _yInt_seES) ///
			clear prefix(_y)
		qui drop *_EFFECT
		qui rename _y_BY _BY
		qui rename _y_STUDY _STUDY
		qui rename _y_LABELS _LABELS
		qui order _y_USE, before(_y_ES)

		qui gen double _yInt_LCI = _yInt_ES - invnormal(.975)*_yInt_seES
		qui gen double _yInt_UCI = _yInt_ES + invnormal(.975)*_yInt_seES
		qui gen double _yInt_WT = 1/_yInt_seES^2 if _y_USE!=5
		summ _yInt_WT, meanonly
		qui replace _yInt_WT = 100*_yInt_WT / r(sum)
		label variable _yInt_WT "% Weight"
		format _yInt_WT %6.2f
		
		// append subgroup estimates
		// ... and "correct" the contents of _STUDY (see above)
		tempvar pooled obs
		qui append using `pooled_ests', gen(`pooled')
		qui gen int `obs' = _n

		qui levelsof _STUDY, local(sgvals)
		forvalues i = 1 / `k' {
			// local x : word `i' of `sgvals`string''
			local x : word `i' of `sgvals'
			summ `obs' if _STUDY==`x' & !`pooled', meanonly
			local labi = _LABELS[`r(min)']
			qui replace _LABELS = `"`labi'"' if `pooled' & _STUDY==`i'	// REVISIT!!  may need to edit -metani-
			qui replace _STUDY = `x' if `pooled' & _STUDY==`i'			// REVISIT!!  may need to edit -metani-
		}
		drop `pooled' `obs'
		
		qui gen byte _yInt_USE = _y_USE
		qui replace _yInt_USE = 9 if _STUDY==1		// reference
		qui replace _yInt_USE = 2 if _yInt_USE==1 & missing(_yInt_seES)
		order _yInt_USE, before(_yInt_ES)
		
		qui compress
				
		if `"`saving'"'!=`""' {
			_prefix_saving `saving'
			qui save `s(filename)', `s(replace)'
		}
		if `"`clear'"'!=`""' {
			restore, not
		}
	}	
		
		
	*** Finally, present the results
	// mat list `BetaHat'
	// mat list `VarBetaHat'
	mat `BetaHat' = `BetaHat''
	
	ereturn post `BetaHat' `VarBetaHat', depname(Subgroup)
	// ereturn matrix GammaHat    = `GammaHat'
	// ereturn matrix VarGammaHat = `VarGammaHat'
	// ereturn matrix ThetaHat    = `ThetaHat'
	// ereturn matrix VarThetaHat = `VarThetaHat'
	ereturn matrix GammaHat    = `GammaHat'
	ereturn matrix VarGammaHat = `VarGammaHat'
	ereturn scalar ThetaHat    = `ThetaHat'
	ereturn scalar VarThetaHat = `VarThetaHat'
	ereturn matrix Grad        = `Grad'
	// ereturn matrix A = `A'
	// ereturn matrix SigmaGamma  = `SigmaGamma'
	// ereturn matrix SigmaBeta   = `SigmaBeta'
	// ereturn scalar n = `t'						// number of studies
	ereturn scalar k = `k'						// number of subgroups
	
	nois disp as text _n "Floating subgroups:"
	ereturn display, `eformopt'	
	
end

