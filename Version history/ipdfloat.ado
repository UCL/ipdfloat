*** IPDFLOAT ***
* Routine for fitting IPD-based "floating subgroup" models
// by David Fisher, Peter Godolphin, Ian White
// (MRC Clinical Trials Unit at UCL)

* Note: ultimately we'll try to redesign this command so it handles modern factor-variable notation
* and hence `command' is the actual interaction-based command you'd want to fit

* For now, though, things are a little more clunky:
* The `command' needs to be in the following form:
*   [command] [depvar] trtvar intvar [adjustment vars] [if] [in] , [cmd_opts]

* (where, as usual with Stata, [depvar] is ONLY optional IF it's e.g. a survival model; otherwise NOT optional)

* ...and where the two other non-optional variables are:
* 1) trtvar = the treatment variable
* 2) intvar = the covariate *TO BE* interacted with treatment (either continuous or categorical)

* For now, need to specify the option `factor' to treat as factor; otherwise will treat as continuous by default

* Reference value or factor can be specified using ref(#)
* If left blank, automatically uses lowest observed factor level, or sample mean of continuous covariate

* v0.1  David Fisher 17sep2021
* v0.2  David Fisher 28sep2021 


program define ipdfloat, eclass
	
	version 11.0
	local version : di "version " string(_caller()) ":"	
	
	_on_colon_parse `0'
	local 0       `"`s(before)'"'
	local command `"`s(after)'"'
	
	syntax, STUDY(varname) POOLvars(varlist numeric min=2 max=2) ///
		[FACTOR SEPGroups REFerence(string) EPS(real 1e-7) FORWards BACKWards Level(cilevel) SAVING(string) CLEAR]

	`version' _prefix_command ipdfloat2, level(`level') : `command'
	tempvar touse
	mark `touse' `s(if)' `s(in)'
	markout `touse' `study' `poolvars'
	
	cap drop _Trt*
	cap drop _Fac*
	
	gettoken trtvar intvar : poolvars
	qui tab `trtvar'
	if `r(r)'!=2 {
		nois disp as err "Treatment variable must have exactly two values"
		exit 198
	}
	summ `trtvar', meanonly
	qui gen byte _Trt = (`trtvar'==r(max))
	
	local k1 = 1
	if "`factor'"!="" {
		cap tab `intvar'
		if _rc | `r(r)'>10 {
			nois disp as err `"Too many factor levels in variable {bf:`intvar'}"'
			exit 198
		}
		local k = r(r)
		local k1 = `k' - 1
		
		if `"`reference'"'==`""' {
			summ `intvar', meanonly
			local ref = r(min)
		}
		else {
			cap confirm integer number `reference'
			if _rc {
				nois disp as err `"With {bf:factor} option, {bf:reference(#)} must be a single integer value"'
				exit 198
			}
			qui count if `intvar'==`reference'
			if !`r(N)' {
				nois disp as err `"With {bf:factor} option, {bf:reference(#)} must be a value taken on by variable {bf:`intvar'}"'
				exit 198
			}
			local ref : copy local reference
		}
					
		qui levelsof `intvar', local(levs)
		foreach i of local levs {
			cap confirm integer number `i'
			if _rc {
				nois disp as err `"With {bf:factor} option, variable {bf:`intvar'} must only take on integer values"'
				exit 198
			}
			if `i'!=`ref' {
				qui gen byte _Fac`i' = (`intvar'==`i')
				qui gen byte _TrtFac`i' = _Trt * _Fac`i'
				
				local fac_coefs `fac_coefs' _Fac`i'
				local int_coefs `int_coefs' _TrtFac`i'
			}
			else {
				qui gen byte _TrtFac`i' = _Trt * (`intvar'==`i')
				local ref_coef _TrtFac`i'
			}
		}
		
		local sglabel : value label `intvar'
		if "`sglabel'"!="" {
			tempfile sglabfile
			qui label save `sglabel' using `sglabfile'
		}
	}
	
	else {
		if `"`saving'"'!=`""' {
			nois disp as err "Cannot use option {bf:saving()} with continuous interactions"
			exit 198
		}
		if `"`clear'"'!=`""' {
			nois disp as err "Cannot use option {bf:clear} with continuous interactions"
			exit 198
		}
		if `"`reference'"'==`""' {
			summ `intvar' if `touse', meanonly
			qui gen double _Fac = `intvar' - r(mean)		// centered at mean
		}
		else {
			cap confirm number `reference'
			if _rc {
				nois disp as err `"Option {bf:reference(#)} must contain a single numeric value"'
				exit 198
			}
			summ `intvar' if `touse', meanonly
			qui gen double _Fac = `intvar' - `reference'	// centered at specified reference
		}
		qui gen double _TrtFac = _Trt * _Fac
		local fac_coefs _Fac
		local int_coefs _TrtFac
	}
	
	
	*** Find interaction effect(s)  ( = contrasts if `factor')
	
	// QUICK AND DIRTY!!
	// remove `trtvar' and `intvar' from `command', to be replaced with the variables defined above
	local command_new = subinstr(`"`command'"', `"`poolvars'"', `"_Trt `int_coefs' `fac_coefs'"', 1)	
	
	preserve
		qui mvmeta_make `command_new' if `touse', usevars(_Trt `int_coefs') names(_yInt _VInt) by(`study') ppfix(check) hard clear

		local vlist				
		foreach v of local int_coefs {
			local vlist `vlist' _yInt`v'
		}
		
		mvmeta _yInt _VInt, vars(`vlist') fixed
		test `vlist'
		
		tempname GammaHat VarGammaHat
		matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix define `GammaHat'    = e(b)[1,       1..`k1']'
		
		if "`saving'"!="" | "`clear'"!="" {
			
			if "`sepgroups'"=="" {
				// qui levelsof `intvar', local(levs)
				foreach i of local levs {
					if `i'==`ref' {
						rename _yInt_Trt _y_ES`i'
						qui gen double _y_seES`i' = sqrt(_VInt_Trt_Trt)
					}
					else {
						qui gen double _y_ES`i' = _y_ES`ref' + _yInt_TrtFac`i'
						qui gen double _y_seES`i' = sqrt(_VInt_Trt_Trt + _VInt_TrtFac`i'_TrtFac`i' + 2*_VInt_Trt_TrtFac`i')
					}
				}
				local extra _y_ES* _y_seES*
				local extra_reshape _y_ES _y_seES
			}
						
			foreach i of local levs {
				if `i'!=`ref' {
					rename _yInt_TrtFac`i' _yInt_ES`i'
					qui gen double _yInt_seES`i' = sqrt(_VInt_TrtFac`i'_TrtFac`i')
				}
			}
			
			keep Trial _yInt_ES* _yInt_seES* `extra'
			qui reshape long _yInt_ES _yInt_seES `extra_reshape', i(Trial) j(Subgroup)
			
			tempfile interactions
			qui save `interactions'
		}
	restore	

	// Treatment effect at reference value / in reference subgroup
	tempvar offset offset_h
	qui gen double `offset_h' = .
	qui gen double `offset' = 0

	tokenize `int_coefs'
	forvalues i = 1 / `k1' {
		qui replace `offset' = `offset' + `GammaHat'[`i',1]*``i''
	}
	
	if `"`factor'"'!=`""' {
		
		// fixed-effects: do it using -ipdmetan- instead of -mvmeta-
		if "`sepgroups'"=="" {
			
			// QUICK AND DIRTY!!
			// remove `trtvar' and `intvar' from `command', to be replaced with the variables defined above
			local command_new = subinstr(`"`command'"', `"`poolvars'"', `"_Trt `fac_coefs'"', 1)
			
			qui ipdmetan_new, study(`study') poolvar(_Trt) nograph norsample : `command_new' if `touse', offset(`offset')
			tempname ThetaHat VarThetaHat
			matrix define `VarThetaHat' = r(se_eff)^2
			matrix define `ThetaHat'    = r(eff)
		}
		
		else {
			// QUICK AND DIRTY!!
			// remove `trtvar' and `intvar' from `command', to be replaced with the variables defined above
			local command_new = subinstr(`"`command'"', `"`poolvars'"', `"`ref_coef' `int_coefs' `fac_coefs'"', 1)	
			
			preserve
				qui mvmeta_make `command_new' if `touse', offset(`offset') usevars(`ref_coef' `int_coefs') names(_y _S) by(`study') ppfix(check) hard clear
				
				local vlist
				local use_coefs `ref_coef' `int_coefs'
				foreach v of local use_coefs {
					local vlist `vlist' _y`v'
				}
				qui mvmeta _y _S, vars(`vlist') fixed commonparm
			
				tempname ThetaHat VarThetaHat
				matrix define `VarThetaHat' = e(V)[1,1]
				matrix define `ThetaHat'    = e(b)[1,1]
				
				if "`saving'"!="" | "`clear'"!="" {
					keep Trial _y* _S*
					qui reshape long _y_ES _y_seES _yInt_ES _yInt_seES, i(Trial) j(Subgroup)					
					tempfile subgroups
					qui save `subgroups'
				}
			restore
		}
	}
	
	else {
		/* STILL TO DO !!! */
		preserve
			qui ipdmetan_new, study(`study') poolvar(_Trt) nograph norsample clear : `command_new' if `touse', offset(`offset')
			tempname ThetaHat VarThetaHat
			matrix define `VarThetaHat' = r(se_eff)^2
			matrix define `ThetaHat'    = r(eff)
			
			rename _STUDY Trial
			keep Trial _ES _seES
			tempfile subgroups
			qui save `subgroups'
		restore
	}

	
	*** Obtain "delta method" gradients, in order to form the final variance-covariance matrix
	tempname Grad h eff_pos eff_neg
	matrix define `Grad' = J(1, `k1', .)

	forvalues i = 1 / `k1' {
		scalar `h' = max(abs(`GammaHat'[1,`i']), 1)*sqrt(`eps')

		local sign "+"
		if "`backwards'"!="" local sign "-"
		else if "`forwards'"=="" scalar `h' = `h'/2
		
		qui replace `offset_h' = `offset' `sign' `h'*``i''
		
		if ("`factor'"!="" & "`sepgroups'"=="") | "`factor'"=="" {
			qui ipdmetan_new, study(`study') poolvar(_Trt) nograph norsample : `command_new' if `touse', offset(`offset_h')
			scalar `eff_pos' = r(eff)
		}
		else {
			preserve
				qui mvmeta_make `command_new' if `touse', offset(`offset_h') usevars(`ref_coef' `int_coefs') names(_y _S) by(`study') ppfix(check) hard clear
				local vlist
				local use_coefs `ref_coef' `int_coefs'
				foreach v of local use_coefs {
					local vlist `vlist' _y`v'
				}
				qui mvmeta _y _S, vars(`vlist') fixed commonparm
				scalar `eff_pos' = e(b)[1,1]
			restore
		}
		
		if `"`forwards'`backwards'"'!=`""' {
			matrix `Grad'[1, `i'] = `sign'(`eff_pos' - `ThetaHat'[1,1])/`h'
			continue
		}

		qui replace `offset_h' = `offset' - `h'*``i''
		
		if ("`factor'"!="" & "`sepgroups'"=="") | "`factor'"=="" {
			qui ipdmetan_new, study(`study') poolvar(_Trt) nograph norsample : `command_new' if `touse', offset(`offset_h')
			scalar `eff_neg' = r(eff)
		}
		else {
			preserve
				qui mvmeta_make `command_new' if `touse', offset(`offset_h') usevars(`ref_coef' `int_coefs') names(_y _S) by(`study') ppfix(check) hard clear
				local vlist
				local use_coefs `ref_coef' `int_coefs'
				foreach v of local use_coefs {
					local vlist `vlist' _y`v'
				}
				qui mvmeta _y _S, vars(`vlist') fixed commonparm
				scalar `eff_neg' = e(b)[1,1]
			restore
		}

		matrix `Grad'[1, `i'] = (`eff_pos' - `eff_neg')/(2*`h')
	}
		

	*** Assemble matrices
	if "`factor'"!="" {
		tempname ones Z A BetaHat VarBetaHat
		
		matrix define `ones' = J(`k', 1, 1)
		matrix define `Z' = J(`k', `k1', 0)
		matrix define `Z'[2, 1] = I(`k1')	
		matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat'

		// matrix define grad = grad2, grad3
		matrix define `A' = `ones'*`Grad' + `Z'
		matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	}
	
	// continuous
	else {
		cap drop _Int_xb
		cap drop _Int_stdp
		cap drop _Int_lci
		cap drop _Int_uci
		qui gen double _Int_xb = `ThetaHat'[1,1]*_Trt + `GammaHat'[1,1]*`1'
		qui gen double _Int_stdp = sqrt(`VarThetaHat'[1,1]*_Trt + `VarGammaHat'[1,1]*((`1' + `Grad'[1,1])^2))
		qui gen double _Int_lci = _Int_xb - invnormal(.5 + `level'/200)*_Int_stdp
		qui gen double _Int_uci = _Int_xb + invnormal(.5 + `level'/200)*_Int_stdp
		
		ereturn matrix GammaHat    = `GammaHat'
		ereturn matrix VarGammaHat = `VarGammaHat'
		ereturn scalar ThetaHat    = `ThetaHat'[1,1]
		ereturn scalar VarThetaHat = `VarThetaHat'[1,1]
		ereturn matrix Grad        = `Grad'		
	}
	
	
	*** Construct saved dataset
	if `"`saving'"'!=`""' | `"`clear'"'!=`""' {
		preserve
		
		// store pooled estimates
		tempfile pooled_ests
		tempname BetaTable GammaTable
		mat `GammaTable' = `GammaHat', vecdiag(`VarGammaHat')'
		local levs_not_ref : list levs - ref
		matrix rownames `GammaTable' = `levs_not_ref'
		qui metani `GammaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
			saving(`pooled_ests', stacklabel) prefix(_yInt)

		mat `BetaTable' = `BetaHat', vecdiag(`VarBetaHat')'
		matrix rownames `BetaTable' = `levs'
		qui metani `BetaTable', variances nograph nooverall rownames rowtitle(Pooled estimates) ///
			clearstack prefix(_y)

		rename _y_LABELS _yInt_LABELS
		
		tempvar sortorder
		qui gen int `sortorder' = _n		// -merge- will sort the data on _yInt_LABELS, so need to re-sort afterwards
		qui merge 1:1 _yInt_LABELS using `pooled_ests', assert(match master)		
		cap {
			assert _merge==1 if _yInt_LABELS=="`ref'"
			assert _merge==3 if _yInt_LABELS!="`ref'"
		}
		if _rc {
			nois disp as err "Error in colnames [NEED TO IMPROVE THIS MESSAGE]"
			exit _rc
		}
		
		drop _merge *_EFFECT _yInt_USE _yInt_STUDY
		rename _y_STUDY _STUDY
		rename _yInt_LABELS _LABELS		
		qui replace _y_USE = 5 if inlist(_y_USE, 1, 2)
		sort `sortorder'
		drop `sortorder'
		qui save `pooled_ests', replace
		
		
		// merge subgroup estimates and interaction estimates into a single dataset
		// use -metan- to restructure
		qui use `interactions', clear
		
		if `"`subgroups'"'!=`""' {
			* stuff
		}
		if "`sglabel'"!="" {
			qui do `sglabfile'
			label values Subgroup `sglabel'
		}		
		
		qui metan _y_ES _y_seES, nograph nohet nosubgroup nooverall keeporder ///
			study(Subgroup) by(`study') rcols(_yInt_ES _yInt_seES) ///
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
		qui replace _yInt_USE = 9 if _STUDY==`ref'		// reference
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
	if "`factor'"!="" {
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
		ereturn scalar ThetaHat    = `ThetaHat'[1,1]
		ereturn scalar VarThetaHat = `VarThetaHat'[1,1]
		ereturn matrix Grad        = `Grad'
		// ereturn matrix A = `A'
		// ereturn matrix SigmaGamma  = `SigmaGamma'
		// ereturn matrix SigmaBeta   = `SigmaBeta'
		// ereturn scalar n = `t'						// number of studies
		ereturn scalar k = `k'						// number of subgroups
		
		nois disp as text _n "Floating subgroups:"
		ereturn display, `eformopt'
	}
	
end

