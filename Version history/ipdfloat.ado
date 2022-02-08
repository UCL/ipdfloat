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
* v0.3  David Fisher 17dec2021
*! v0.3.1  David Fisher 05jan2022 -- intermediate test version with some options/routines removed

* Note DF 17th Dec 2021: needs to use -ipdmetan- version 4.03  David Fisher  29nov2021
* (as yet unreleased)

program define ipdfloat, eclass
	
	version 11.0
	local version : di "version " string(_caller()) ":"	
	
	// Check that -metan- v4.05+ and -ipdmetan- v4.03+ are installed
	cap metan
	if `"`r(metan_version)'"'==`""' {
		nois disp as err `"{bf:ipdfloat} requires version 4.05 or higher of {bf:metan}"'
		exit 499
	}
	else {
		local reqd_version = 4.05
		if `r(metan_version)' < `reqd_version' {
			nois disp as err `"{bf:ipdfloat} requires version 4.05 or higher of {bf:metan}"'
			exit 499
		}
	}
	cap ipdmetan
	if `"`r(ipdmetan_version)'"'==`""' {
		nois disp as err `"{bf:ipdfloat} requires version 4.03 or higher of {bf:ipdmetan}"'
		exit 499
	}
	else {
		local reqd_version = 4.03
		if `r(ipdmetan_version)' < `reqd_version' {
			nois disp as err `"{bf:ipdfloat} requires version 4.03 or higher of {bf:ipdmetan}"'
			exit 499
		}
	}
	
	_on_colon_parse `0'
	local 0       `"`s(before)'"'
	local command `"`s(after)'"'
	
	syntax, STUDY(varname) ///
		[ TRTVAR(varname) INTVARs(varlist) ///
		  SHOWmodels HARD SEPGroups EPS(real 1e-7) FORWards BACKWards Level(cilevel) SAVING(string) CLEAR ]

	`version' _prefix_command ipdfloat, level(`level') : `command'
	local cmdname `"`s(cmdname)'"'
	local cmdopts `"`s(options)'"'
	local rest `"`s(anything)'"'	
	
	tempvar touse
	mark `touse' `s(if)' `s(in)'
	markout `touse' `study', strok
	
	
	*** Find treatment and interaction effect(s)  ( = contrasts if `factor')
	// Assume first element of `command' is of the form: i.`trtvar'*`i'`intvar'
	if substr(`"`cmdname'"', 1, 2)!="st" {		// survival commands don't have a dependent variable
		gettoken depvar rest : rest
		// local depvar : copy local cmdname
		// gettoken cmdname : rest
	}
	
	tokenize `"`rest'"', parse(" .*")
	local abbtrtvar : copy local 3
	unab fulltrtvar : `3'
	capture {
		assert `"`1'`2'"'==`"i."'
		assert `"`4'"'==`"*"'
	}
	if _rc {
		nois disp as err  "Could not identify treatment and subgroup variables"
		nois disp as err `"First element of {it:command} should be of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
		exit 198
	}
	
	if `"`5'`6'"'==`"i."' {
		local factor factor
		unab unab7 : `7'
		fvunab other : `8'
		
		if `"`8'"'==`"*"' {
			nois disp as err `"{bf:ipdfloat} currently only supports two-way treatment-covariate interaction terms"'
			exit 198
		}
	}
	else {
		// if continuous, `5' should contain the (first) covariate interaction variable
		unab unab7 : `5'
	
		// also need to search for any other terms beginning "i.trt*[covariate]"
		// as these will not be included in `VarsToDrop' but will be needed to form `intvars'
		local abbtrtlen = length(`"`abbtrtvar'"')
		gettoken el newrest : rest
		local k = 1						// for treatment
		while `"`el'"'!=`""' {
			tokenize `"`el'"', parse(".*")
			if `"`1'`2'"'==`"i."' & `"`4'"'==`"*"' {
				local len3 = length(`"`3'"')
				if substr(`"`3'"', 1, min(`len3', `abbtrtlen'))==substr(`"`abbtrtvar'"', 1, min(`len3', `abbtrtlen')) {
					unab foo : `5'
					local covars `covars' `foo'
					local ++k
				}
			}
			else local other `other' `el'
			gettoken el newrest : newrest
		}		
		fvunab other : `other'
	}
	
	qui xi, prefix(_J) `rest'
	local OldVarsPrefix : char _dta[__xi__Vars__Prefix__]
	local OldVarsToDrop: char _dta[__xi__Vars__To__Drop__]
	local i = 0
	foreach el of local OldVarsPrefix {
		local ++i
		if `"`el'"'==`"_J"' {
			local eli : word `i' of `OldVarsToDrop'
			local VarsToDrop `VarsToDrop' `eli'
		}
	}
	
	// Process treatment variable
	gettoken trtvar rest : VarsToDrop			// first element should be of the form _J[trtvar]_1
	local trtlen = length(`"`trtvar'"')
	local trtvar7 = substr(`"`trtvar'"', 3, min(7, `=`trtlen'-4'))	// isolates [trtvar] truncated to 7 characters (so that _J[trtvar]_1 is 11 characters max)
	
	cap	assert `"`trtvar7'"'==substr(`"`fulltrtvar'"', 1, `trtlen')
	if _rc {
		nois disp as err "Error in isolating treatment variable"
		exit _rc
	}	
	cap tab `trtvar7' if `touse'									// this will also trap varname ambiguities
	if _rc==111 {
		nois disp as err "Treatment variable {bf:`trtvar7'} must be uniquely identifiable by its first 7 characters"
	}
	if _rc exit _rc
	cap assert `r(r)'==2
	if _rc {
		nois disp as err "Treatment variable {bf:`trtvar7'} must take exactly two values"
		exit _rc
	}	
	
	local _Jtrt = "_J" + substr(`"`trtvar7'"', 1, 3)	// isolates the string "_J[trt]" where [trt] = [trtvar] truncated to 3 characters
	tokenize `VarsToDrop'
	while `"`1'"'!=`""' {
		if substr(`"`1'"', 1, 5)==`"`_Jtrt'"' {
			local usevars `usevars' `1'			// `usevars' is a subset of `VarsToDrop', beginning "_J[trt]X..."  (or _J[trtvar] itself)
		}
		else {									// `facvars' contain the remainder...
			local facvars `facvars' `1'			// ...should be the main-effects of the factor-variable covariate (if applicable)
		}
		macro shift
	}
	if `"`usevars'"'==`""' {
		nois disp as err  "Could not identify treatment and subgroup variables"
		nois disp as err `"First element of {it:command} should be of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
		exit 198
	}
	

	// Process interaction variable
	gettoken intvar : rest
	cap assert `"`intvar'"'!=`""'
	if _rc {
		nois disp as err "Subgroup variable could not be identified"
		exit 198
	}

	// if factor, second element should be of the form _J[intvar]_1
	// else if continuous, "xi:" will not create a new variable for `intvar', so...
	//    ...second element will be the interaction itself, of the form _J[trtvar3]X[intvar]
	local intlen = length(`"`intvar'"')
	local _Jtrtlen = length(`"`_Jtrt'"')
	local intvar7 = substr(`"`intvar'"', 3, min(7, `=`intlen'-4'))	// isolates [intvar] truncated to 7 characters (so that _J[intvar]_1 is 11 characters max)
	cap assert (`"`factor'"'==`""') == (substr(`"`intvar'"', 1, min(5, `_Jtrtlen'))==`"`_Jtrt'"')
	if _rc {
		nois disp as err "Subgroup variable could not be identified"
		exit 198
	}
	
	if `"`factor'"'!=`""' {
		cap tab `intvar7' if `touse'	// this will trap varname ambiguities
		local k = r(r)
		local k1 = `k' - 1
		if _rc==111 {
			nois disp as err "Subgroup variable {bf:`intvar7'} must be uniquely identifiable by its first 7 characters"
		}
		else if `k'>10 {
			nois disp as err `"Too many factor levels in variable {bf:`intvar7'} (max of 10)"'
			exit 198
		}
		if _rc exit _rc
				
		// find reference value
		unab intvar : `intvar7'
		local ref : char `intvar'[omit]
		if `"`ref'"'==`""' {
			summ `intvar', meanonly
			local ref = r(min)
		}
		else {
			confirm integer number `ref'
			qui count if `intvar'==`ref'
			if `r(N)' {
				summ `intvar', meanonly
				nois disp `"{error}Reference value `ref' is not taken on by variable {bf:`intvar'} in any observations currently in memory"'
				nois disp `"{error}Reference value will be re-set to `r(min)'"'
				local ref = r(min)
			}
		}
	}
	
	// continuous `intvar'
	else {
		cap assert `"`facvars'"'==`""'
		if _rc {
			nois disp as err "Error in isolating covariate variable(s)"
			exit 198
		}
		
		if `"`saving'"'!=`""' {
			nois disp as err "Cannot use option {bf:saving()} with continuous interactions"
			exit 198
		}
		if `"`clear'"'!=`""' {
			nois disp as err "Cannot use option {bf:clear} with continuous interactions"
			exit 198
		}	
		if `"`sepgroups'"'!=`""' {
			nois disp as err "Cannot use option {bf:sepgroups} with continuous interactions"
			exit 198
		}	
	
		local intvarlab : variable label `intvar'
		tokenize `intvarlab', parse("*")
		cap {
			assert `"`2'"'==`"*"'
			assert `"`4'"'==`""'
			summ `3', meanonly
		}
		if _rc {
			nois disp as err "Error in identification of treatment-by-subgroup interaction terms"
			exit 198
		}
		local intvar : copy local 3
		nois disp `"Note: Interaction variable {bf:`intvar'} will be treated as continuous"'
		
		// find reference value
		local ref : char `intvar'[omit]
		if `"`ref'"'==`""' {
			summ `intvar', meanonly
			local ref = r(mean)
		}
		else {
			confirm number `ref'
			summ `intvar', meanonly
			if !inrange(`ref', `r(min)', `r(max)') {
				summ `intvar', meanonly
				nois disp `"{error}Reference value `ref' does not lie within the range of variable {bf:`intvar'} among the observations currently in memory"'
				nois disp `"{error}Reference value will be re-set to `r(mean)'"'
				local ref = r(mean)
			}
		}		
	}

	cap assert `: word count `usevars'' == `k'
	if _rc {
		nois disp as err "Error in identification of treatment-by-subgroup interaction terms"
		exit 198
	}
	
	
	// Note:
	// `command' should be the same as:
	// `cmdname' `depvar' `trtvar' `intvars' `facvars'/`covars' `other'
	// (where `"`usevars'"' == `"`trtvar' `intvars'"')
	if "`showmodels'"=="" local quietly quietly
	
	// Interaction effects (contrasts)
	preserve
		qui xi, prefix(_J) : mvmeta_make `command' if `touse', `cmdopts' usevars(`usevars') /*names(_yInt _SInt)*/ by(`study') ppfix(check) `hard' names(y S) `collapse_opt' nodetails clear
		local yvars : colnames e(b)
		
		local allvars `"`usevars' `facvars'"'
		local others : list yvars - allvars
		
		gettoken trtvar intvars : usevars
		local vlist
		local k_check = 1
		foreach v of local intvars {
			local vlist `vlist' y`v'
			local ++k_check
		}
		assert `k_check'==`k'
		
		nois disp as text _n "Test of interaction(s):"
		`quietly' mvmeta y S, vars(`vlist') fixed
		nois test `vlist'
		local chi2_int = r(chi2)
		local df_int = r(df)
		local p_int = r(p)

		local k1 = `k' - 1
		tempname GammaHat VarGammaHat		
		matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix define `GammaHat'    = e(b)[1,       1..`k1']

		if `"`saving'"'!=`""' | `"`clear'"'!=`""' {
			
			tokenize intvars
			if `"`sepgroups'"'==`""' {
				qui levelsof `intvar', local(levs)
				rename y_J`trtvar' _y_ES`ref'
				qui gen double _y_seES`ref' = sqrt(S_J`trtvar'_J`trtvar')
				
				foreach i of local levs {
					if `i'!=`ref' {
						qui gen double _y_ES`i' = _y_ES`ref' + y_J``i''
						qui gen double _y_seES`i' = sqrt(S_J`trtvar'_J`trtvar' + S_J``i''_J``i'' + 2*S_J`trtvar'_J``i'')
					}
				}
				local extra _y_ES* _y_seES*
				local extra_reshape _y_ES _y_seES
			}
						
			foreach i of local levs {
				if `i'!=`ref' {
					rename y_J``i'' _yInt_ES`i'
					qui gen double _yInt_seES`i' = sqrt(S_J``i''_J``i'')
				}
			}
			
			keep Trial _yInt_ES* _yInt_seES* `extra'
			qui reshape long _yInt_ES _yInt_seES `extra_reshape', i(Trial) j(Subgroup)
			
			tempfile interactions
			qui save `interactions'
		}

	restore
	
	// Treatment effect at reference value (if continuous), or in reference subgroup (if factor)
	tempvar offset offset_h
	qui gen double `offset_h' = .
	qui gen double `offset' = 0

	tokenize `intvars'
	forvalues i = 1 / `k1' {
		qui replace `offset' = `offset' + `GammaHat'[1, `i']*``i''
	}

	// setup offset
	unabcmd `cmdname'
	if inlist(`"`r(cmd)'"', "regress", "mixed") {
		tempvar depvar_offset
		qui gen double `depvar_offset' = `depvar' - `offset'
		local depvar2 : copy local depvar_offset
	}
	else {
		local offset_opt offset(`offset')
		local depvar2 : copy local depvar
	}

	// NOTE: commands allowing offset() option:
	// binreg, biprobit, clogit, cloglog, cpoisson, expoisson, fracreg, glm, heckman, heckoprobit, heckpoisson, heckprobit, hetoprobit, hetprobit,
	// intreg, ivpoisson, logistic, logit, nbreg, gnbreg, ologit, oprobit, poisson, probit, scobit, tnbreg, tobit, tpoisson, truncreg, zinb, ziologit, zioprobit, zip
	
	// for now, just assume everything except -regress-, -xtreg- and -mixed-
	
	if `"`sepgroups'"'==`""' | `"`factor'"'==`""' {
		
		// remove `intvars' from `command'		
		local command_new `cmdname' `depvar2' `trtvar' `facvars' `covars' `other'
		local usevars `trtvar'
		
		// fixed-effects: do it using -ipdmetan- instead of -mvmeta-
		/*
		qui ipdmetan, study(`study') poolvar(`trtvar') nograph norsample : `command_new' if `touse', `offset_opt' `cmdopts'
		tempname ThetaHat VarThetaHat
		matrix define `VarThetaHat' = r(se_eff)^2
		matrix define `ThetaHat'    = r(eff)
		*/
		
		/*
		if `"`factor'"'==`""' {
			rename _STUDY Trial
			keep Trial _ES _seES
			tempfile subgroups
			qui save `subgroups'
		}
		*/
	}
		
	else {
		// replace `trtvar' with `refintvar'
		tempvar refintvar
		qui gen byte `refintvar' = `trtvar' * (`intvar'==`ref')
		
		local command_new `cmdname' `depvar2' `refintvar' `intvars' `facvars' `covars' `other'
		local usevars `refintvar' `intvars'
		
		/*
		preserve
			qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`useintvars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
			mvmeta y S, fixed commonparm
		
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
		*/
	}
	
	preserve
		qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`trtvar') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
		mvmeta y S, fixed commonparm

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
	

	*** Obtain "delta method" gradients, in order to form the final variance-covariance matrix
	tempname Grad h eff_pos eff_neg
	matrix define `Grad' = J(1, `k1', .)

	forvalues i = 1 / `k1' {
		scalar `h' = max(abs(`GammaHat'[1, `i']), 1)*sqrt(`eps')

		local sign "+"
		if `"`backwards'"'!=`""' local sign "-"
		else if `"`forwards'"'==`""' scalar `h' = `h'/2
	
		// update offset
		qui replace `offset_h' = `offset' `sign' `h'*``i''
		unabcmd `cmdname'
		if inlist(`"`r(cmd)'"', "regress", "xtreg", "mixed") {
			qui replace `depvar_offset' = `depvar' - `offset_h'
		}
		else {
			local offset_opt offset(`offset_h')
		}
		
		/*
		if `"`sepgroups'"'==`""' | `"`factor'"'==`""' {
			qui ipdmetan, study(`study') poolvar(`trtvar') nograph norsample : `command_new' if `touse', `offset_opts' `cmdopts'
			scalar `eff_pos' = r(eff)
		}
		else {
		*/
			preserve
				qui mvmeta_make `command_new' if `touse', `offset_opts' `cmdopts' usevars(`usevars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
				qui mvmeta y S, fixed commonparm
				scalar `eff_pos' = e(b)[1,1]
			restore
		// }
		
		if `"`forwards'`backwards'"'!=`""' {
			matrix `Grad'[1, `i'] = `sign'(`eff_pos' - `ThetaHat'[1,1])/`h'
			continue
		}

		// update offset
		qui replace `offset_h' = `offset' - `h'*``i''
		unabcmd `cmdname'
		if inlist(`"`r(cmd)'"', "regress", "mixed") {
			qui replace `depvar_offset' = `depvar' - `offset_h'
		}
		
		/*
		if "`sepgroups'"=="" | "`factor'"=="" {
			qui ipdmetan, study(`study') poolvar(`trtvar') nograph norsample : `command_new' if `touse', `offset_opts' `cmdopts'
			scalar `eff_neg' = r(eff)
		}
		else {
		*/
			preserve
				qui mvmeta_make `command_new' if `touse', `offset_opts' `cmdopts' usevars(`usevars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
				qui mvmeta y S, fixed commonparm
				scalar `eff_neg' = e(b)[1,1]
			restore
		// }

		matrix `Grad'[1, `i'] = (`eff_pos' - `eff_neg')/(2*`h')
	}
		

	*** Assemble matrices
	if "`factor'"!="" {
		tempname ones Z A BetaHat VarBetaHat
		
		matrix define `ones' = J(`k', 1, 1)
		matrix define `Z' = J(`k', `k1', 0)
		matrix define `Z'[2, 1] = I(`k1')
		matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat''

		matrix define `A' = `ones'*`Grad' + `Z'
		matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
	}
	
	// continuous
	else {
		cap drop _Int_xb
		cap drop _Int_stdp
		cap drop _Int_lci
		cap drop _Int_uci
		
		qui gen double _Int_xb = `ThetaHat'[1,1]*`trtvar'
		qui gen double _Int_stdp = .
		tokenize `intvars'
		forvalues i = 1 / `k1' {
			qui replace _Int_xb = _Int_xb + `GammaHat'[1, `i']*``i''
		}
		
		mata: st_view(X=., ., "`trtvar'", "`touse'")
		mata: st_view(Z=., ., "`intvars'", "`touse'")
		mata: VarThetaHat = st_matrix("`VarThetaHat'")
		mata: VarGammaHat = st_matrix("`VarGammaHat'")
		mata: GradZ = st_matrix("`Grad'")
		mata: GradZ = Z + J(rows(Z), 1, 1) * GradZ
		mata: Var_Xb = diag(VarThetaHat*X) + GradZ*VarGammaHat*GradZ'
		mata: st_store(st_viewobs(X), "_Int_stdp", sqrt(diagonal(Var_Xb)))

		// qui gen double _Int_xb = `ThetaHat'[1,1]*`trtvar' + `GammaHat'[1,1]*`intvar'
		// qui gen double _Int_stdp = sqrt(`VarThetaHat'[1,1]*`trtvar' + `VarGammaHat'[1,1]*((`intvar' + `Grad'[1,1])^2))
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
	
	ereturn local command     `"`command'"'
	ereturn local command_new `"`command_new'"'

	// qui drop `VarsToDrop'
	// char define _dta[__xi__Vars__Prefix__] `OldVarsPrefix'
	// char define _dta[__xi__Vars__To__Drop__] `OldVarsToDrop'	
	
end

