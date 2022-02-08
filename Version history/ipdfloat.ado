*** IPDFLOAT ***
* Routine for fitting IPD-based "floating subgroup" models
// by David Fisher, Peter Godolphin, Ian White
// (MRC Clinical Trials Unit at UCL)

* Syntax:
*   ipdfloat, study(studyname) [options] : cmdname [depvar] i.trtvar*[i.]intvar [adjustment vars] [if] [in] [wgt] , [cmd_opts]

* where:
*  [depvar] = dependent variable (ONLY optional IF it's e.g. a survival model; otherwise NOT optional)
*  trtvar = treatment variable, taking on precisely two values (e.g. 0, 1)
*  intvar = covariate to be interacted with treatment
*    - if factor, write i.trt*i.intvar and only a single intvar is allowed
*    - if continuous, write i.trt*intvar1 [i.trt*intvar2 ...] and multiple intvars are allowed (e.g. to specify a FP or spline)

* Note: -ipdfloat- uses old-style "xi"-based factor notation, although "xi" does not need to be explicitly written within the command.
* (ultimately we may try to redesign this command so it handles modern factor-variable notation, but no immediate need)

* Reference level (for factor variable) or reference value (for first continuous variable) can be specified in advance using "char intvar[omit] #"; see -help char-
* Otherwise, -ipdfloat- will automatically use the lowest observed factor level, or the sample mean of (first) continuous covariate

* v0.1  David Fisher 17sep2021
* v0.2  David Fisher 28sep2021 
* v0.3  David Fisher 17dec2021
*! v0.3.2  David Fisher 07jan2022 -- intermediate test version with some options/routines removed

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
		[ TRTVAR(varname numeric) INTVARs(varlist numeric) ///
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
	
	if `"`trtvar'"'==`""' {
		tokenize `"`rest'"', parse(" .*")
		unab trtvar : `3'
		capture {
			assert `"`1'`2'"'==`"i."'
			assert `"`4'"'==`"*"'
		}
		if _rc {
			nois disp as err  "Could not identify treatment variable. If option {bf:trtvar()} is not specified,"
			nois disp as err `" then first element of {it:command} should be of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
			exit 198
		}
	}
	cap tab `trtvar' if `touse'
	if _rc {
		tab `trtvar' if `touse'
	}
	cap assert `r(r)'==2
	if _rc {
		nois disp as err "Treatment variable {bf:`trtvar'} must take exactly two values"
		exit _rc
	}
	
	// find reference value for treatment variable
	local trtref : char `trtvar'[omit]
	if `"`trtref'"'==`""' {
		summ `trtvar', meanonly
		local trtref = r(min)
	}
	else {
		confirm integer number `trtref'
		qui count if `trtvar'==`trtref'
		if !`r(N)' {
			summ `trtvar', meanonly
			nois disp `"{error}Reference value `trtref' is not taken on by variable {bf:`trtvar'} in any observations currently in memory"'
			nois disp `"{error}Reference value will be re-set to `r(min)'"'
			local trtref = r(min)
		}
	}	
	
	// Now search `rest' for terms beginning "i.trt* ... "
	local specintvars : copy local intvars
	local remain : copy local intvars
	gettoken el newrest : rest
	local k = 1						// for treatment
	while `"`el'"'!=`""' {
		tokenize `"`el'"', parse(".*|")
		
		if lower(`"`1'`2'"')==`"i."' {
			unab newvar : `3'
			local facvars `facvars' `el'
			if `"`factor'"'!=`""' {
				gettoken el newrest : newrest
				continue
			}
		
			if inlist(`"`4'"', `"*"', `"|"') & `"`newvar'"'==`"`trtvar'"' {
				
				// factor
				if lower(`"`5'`6'"')==`"i."' {
					unab newvar : `7'
					if (`: list newvar in intvars') & (`: word count `intvars'' > 1) {
						nois disp as err `"If option {bf:intvars()} is specified, and contains a factor variable, then no other variables may be specified to this option"'
						exit 198
					}
					if `"`intvars'"'==`""' | `"`newvar'"'==`"`intvars'"' {
						local remain : list remain - newvar
						local intvars : copy local newvar
						local factor factor
						// local adjvars `adjvars' `newrest'
						gettoken el newrest : newrest
						continue
					}
					// else local adjvars `adjvars' `el'
				}
				
				// continuous
				else {
					unab newvar : `5'
					if `"`intvars'"'==`""' | `: list newvar in intvars' {
						local remain : list remain - newvar
						local covars `covars' `newvar'
						local ++k
					}
					// else local adjvars `adjvars' `el'
				}
			}
			// else local adjvars `adjvars' `el'
		}
		else local adjvars `adjvars' `el'
		gettoken el newrest : newrest
	}
	
	if `"`remain'"'!=`""' {
		nois disp as err `"Failed to identify interaction(s) between treatment variable {bf:`trtvar'} and covariate(s) {bf:`remain'}"'
		nois disp as err `"Please check {it:command} and {bf:intvars()} option, and try again"'
		exit 198
	}
	
	if `"`factor'"'!=`""' {				// factor-variable; recall that `intvars' must contain only a single varname
		cap tab `intvars' if `touse'
		if _rc {
			tab `intvars' if `touse'
		}
		local k = r(r)		// can overwrite previously-defined local k
		if `k'>10 {
			nois disp as err `"Too many factor levels in variable {bf:`intvars'} (max of 10)"'
			exit 198
		}
				
		// find reference value
		local ref : char `intvars'[omit]
		if `"`ref'"'==`""' {
			summ `intvars', meanonly
			local ref = r(min)
		}
		else {
			confirm integer number `ref'
			qui count if `intvars'==`ref'
			if !`r(N)' {
				summ `intvars', meanonly
				nois disp `"{error}Reference value `ref' is not taken on by variable {bf:`intvars'} in any observations currently in memory"'
				nois disp `"{error}Reference value will be re-set to `r(min)'"'
				local ref = r(min)
			}
		}
	}
	
	else {		// continuous
		if `"`intvars'"'!=`""' {
			local problem : list intvars - covars
			if `"`problem'"'!=`""' {
				local p1 : word 1 of `problem'
				nois disp as err `"Could not find term {bf:i.`trtvar'*`p1'} in {it:command}"'
				exit 198
			}
		}
		else {
			local intvars : copy local covars	// use all of the identified `covars'
		}
		
		// find reference value
		local intvar1 : word 1 of `intvars'
		local ref : char `intvar1'[omit]
		if `"`ref'"'==`""' {
			summ `intvar1', meanonly
			local ref = r(mean)
		}
		else {
			confirm number `ref'
			summ `intvar1', meanonly
			if !inrange(`ref', `r(min)', `r(max)') {
				summ `intvar1', meanonly
				nois disp `"{error}Reference value `ref' does not lie within the range of variable {bf:`intvar1'} among the observations currently in memory"'
				nois disp `"{error}Reference value will be re-set to `r(mean)'"'
				local ref = r(mean)
			}
		}
		
		// options not allowed with continuous interactions
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
	}
	local k1 = `k' - 1

	// so far, macros defined are:
	// `covars' (only if continuous; if factor, defined below)
	// `trtvar', `intvars'
	// `facvars', `adjvars'

	// `covars' (if factor)
	// `adjfacvars'
	// `usetrtvar', `useintvars'
	
	// Now run -xi- on the supplied command
	// and use the results to check again for consistency and for correct variables to use
	qui xi, prefix(_J) `rest'
	local OldVarsPrefix : char _dta[__xi__Vars__Prefix__]
	local OldVarsToDrop: char _dta[__xi__Vars__To__Drop__]
	local i = 0
	foreach el of local OldVarsPrefix {
		local ++i
		if `"`el'"'==`"_J"' {
			local eli : word `i' of `OldVarsToDrop'
			local VarsToDrop `VarsToDrop' `eli'
			
			local intvarlab : variable label `eli'
			tokenize `"`intvarlab'"', parse("*&")
			args trtpart symbol intpart rest
			
			// parse first part (treatment)
			local open
			local close
			tokenize `"`trtpart'"', parse("(=)")
			if `"`1'"'!=`"("' {
				cap {
					assert `"`factor'"'!=`""' | `"`symbol'"'==`""'
					assert `"`2'"'==`"=="'
					assert `"`4'"'==`""'
				}
				args trtname eq trtval					// format:  trtname==trtval [& intname==intval]
			}
			else {
				cap {
					assert `"`1'`3'`5'"'==`"(==)"'
					assert `"`6'"'==`""'
				}
				args open trtname eq trtval close		// format:  (trtname==trtval)*intname
			}
			if _rc ExitErr
			
			// parse second part (interaction)
			tokenize `"`intpart'"', parse("=")
			args intname eq intval
			cap {
				assert `"`rest'"'==`""'
				assert inlist(`"`symbol'"', `""', `"*"', `"&"')
				assert `"`symbol'"'==`""' if `"`intname'"'==`""'
				assert `"`symbol'"'==`"*"' & `"`intname'"'!=`""' & `"`eq'"'==`""' if `"`open'"'!=`""'
				assert `"`intval'"'!=`""' & `"`eq'"'==`"=="' if `"`open'"'==`""' & `"`symbol'"'==`"&"'
				assert `"`4'"'==`""'
			}
			if _rc ExitErr
			// if not treatment (and also not interaction with treatment), should be either:
			//  - a main effect of a factor-variable covariate, or
			//  - part of another interaction which is being adjusted-for but is not relevant to the procedure.
			if `"`trtname'"'!=`"`trtvar'"' {
				cap assert `"`intname'"'!=`"`trtvar'"'
				if _rc {
					nois disp as err "Cannot specify an interaction term where the treatment variable is in second position"
					exit 198
				}
				if `"`factor'"'!=`""' & `"`trtname'"'==`"`intvars'"' {
					cap assert `trtval'!=`ref'
					if _rc ExitErr
					if `"`intpart'"'==`""' local covars `covars' `eli'
					else {
						// (non-treatment) interaction with continuous covariate; need to save main-effect of this covariate
						if `"`symbol'"'==`"*"' local adjfacvars `adjfacvars' `intname'
						local adjfacvars `adjfacvars' `eli'
					}
				}
				else {
					// (non-treatment) interaction with continuous covariate; need to save main-effect of this covariate
					if `"`symbol'"'==`"*"' local adjfacvars `adjfacvars' `intname'
					local adjfacvars `adjfacvars' `eli'
				}
				continue
			}
			cap assert `trtval' != `trtref'
			if _rc ExitErr
					
			// if treatment variable alone (i.e. *not* interaction), then continue to process next element
			if `"`symbol'"'==`""' {
				local usetrtvar : copy local eli
				continue
			}
	
			// at this point we definitely have a treatment-covariate interaction
			if `"`factor'"'!=`""' {
				
				// if treatment interacted with a continuous covariate
				if `"`open'"'!=`""' & `"`specintvars'"'==`""' {
					nois disp as err "Cannot interact treatment with both a factor and a continuous variable"
					exit 198
				}
				
				cap assert `"`intname'"'==`"`intvars'"'
				if _rc {
					// if user-specified `intvars', then ignore it and add to `adjfacvars'
					if `"`specintvars'"'!=`""' {
						local adjfacvars `adjfacvars' `eli'
						continue
					}
					
					// otherwise, exit with error
					else {
						nois disp as err "Cannot interact treatment with more than one factor variable"
						exit 198
					}
				}
				cap {
					assert `"`eq'"'==`"=="'
					assert `intval' != `ref'
				}
			}
			else {
				cap confirm variable `intpart'
				local intname : copy local intpart
			}
			if _rc ExitErr
			
			// if successful, add to `useintvars' (unless excluded from user-specified `intvars')
			if `"`specintvars'"'==`""' | `: list intname in specintvars' {
				local useintvars `useintvars' `eli'
				local intnames `intnames' `intname'
			}		
			else local adjfacvars `adjfacvars' `eli'
		}	
	}

	cap {
		confirm variable `trtvar'
		confirm variable `intvars'
		confirm numeric variable `usetrtvar'
		confirm numeric variable `useintvars'
		assert `: word count `useintvars'' == `k1'
		
		if `"`factor'"'==`""' {
			assert `: list intvars === intnames'		// this is trivial; remove??
			assert `: word count `intvars''==`: word count `useintvars''
		}
	}
	if _rc ExitErr
	
	// to form `command':
	// take `adjvars' and replace `facvars' with `adjfacvars'
	
	local adjfacvars : list uniq adjfacvars
	
	disp `"trtvar: `trtvar'"'
	disp `"intvars: `intvars'"'
	disp `"facvars: `facvars'"'
	disp `"adjvars: `adjvars'"'
	disp `"adjfacvars: `adjfacvars'"'
	disp `"covars: `covars'"'						// main effects of interacted covars
	disp `"usevars: `usetrtvar' `useintvars'"'	
	disp `"k = `k'"'
	
	disp _n `"-----------"'
	
	disp trim(itrim(`"`command'"'))
	disp trim(itrim(`"`cmdname' `depvar' `facvars' `adjvars'"'))
	local cmd2 `"`cmdname' `depvar' `facvars' `adjvars'"'
	assert `: list command === cmd2'
	disp `"`cmdname' `depvar' `usetrtvar' `useintvars' `covars' `adjfacvars' `adjvars'"'
	
	qui `cmdname' `depvar' `usetrtvar' `useintvars' `covars' `adjfacvars' `adjvars'
	local eN_2 = e(N)
	local eN_fail_2 = e(N_fail)
	local ell_0_2 = e(ll_0)
	local ell_2 = e(ll)

	_ms_omit_info e(b)
	mat omit_info_2 = r(omit)
	local yvars2 : colnames e(b)
	local c : word count `yvars2'
	forvalues i=1/`c' {
		local omit = omit_info_2[1, `i']
		if !`omit' {
			local newyvars2 `newyvars2' `: word `i' of `yvars2''
		}
	}

	qui xi, prefix(_J): `command'
	local eN_1 = e(N)
	local eN_fail_1 = e(N_fail)
	local ell_0_1 = e(ll_0)
	local ell_1 = e(ll)

	_ms_omit_info e(b)
	mat omit_info_1 = r(omit)
	local yvars1 : colnames e(b)
	local c : word count `yvars1'
	forvalues i=1/`c' {
		local omit = omit_info_1[1, `i']
		if !`omit' {
			local newyvars1 `newyvars1' `: word `i' of `yvars1''
		}
	}

	disp `"`newyvars1' == `newyvars2'"'
	assert `: list newyvars1 === newyvars2'
	disp `"`eN_1' == `eN_2'"'
	assert `eN_1'==`eN_2'
	disp `"`eN_fail_1' == `eN_fail_2'"'
	assert `eN_fail_1'==`eN_fail_2'
	disp `"`ell_0_1' == `ell_0_2'"'
	assert float(`ell_0_1')==float(`ell_0_2')
	disp `"`ell_1' == `ell_2'"'
	assert float(`ell_1')==float(`ell_2')
	
	exit 0
	// UP TO HERE
	
	// Note:
	// `command' should be the same as:
	// `cmdname' `depvar' `usetrtvar' `useintvars' `covars' `adjvars'
	if "`showmodels'"=="" local quietly quietly
	
	// Interaction effects (contrasts)
	preserve
		qui xi, prefix(_J) : mvmeta_make `command' if `touse', `cmdopts' ///
			usevars(`usetrtvar' `useintvars') /*names(_yInt _SInt)*/ by(`study') ppfix(check) `hard' `collapse_opt' nodetails clear
		local yvars : colnames e(b)
		
		disp `"`yvars'"'
		
		nois disp as text _n "Test of interaction(s):"
		local vlist
		foreach v of local useintvars {
			local vlist `vlist' y`v'
		}
		`quietly' mvmeta y S, vars(`vlist') fixed
		nois test `vlist'
		local chi2_int = r(chi2)
		local df_int = r(df)
		local p_int = r(p)

		tempname GammaHat VarGammaHat		
		matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix define `GammaHat'    = e(b)[1,       1..`k1']

		if `"`saving'"'!=`""' | `"`clear'"'!=`""' {
			
			tokenize useintvars
			if `"`sepgroups'"'==`""' {
				qui levelsof `intvars', local(levs)		// `intvars' is the (single) original factor variable
				rename y_`usetrtvar' _y_ES`ref'
				qui gen double _y_seES`ref' = sqrt(S_`usetrtvar'_`usetrtvar')
				
				foreach i of local levs {
					if `i'!=`ref' {
						qui gen double _y_ES`i' = _y_ES`ref' + y_``i''
						qui gen double _y_seES`i' = sqrt(S_`usetrtvar'_`usetrtvar' + S_``i''_``i'' + 2*S_`usetrtvar'_``i'')
					}
				}
				local extra _y_ES* _y_seES*
				local extra_reshape _y_ES _y_seES
			}
						
			foreach i of local levs {
				if `i'!=`ref' {
					rename y_``i'' _yInt_ES`i'
					qui gen double _yInt_seES`i' = sqrt(S_``i''_``i'')
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

	tokenize `useintvars'
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
		
		// Note:
		// `command' should be the same as:
		// `cmdname' `depvar' `usetrtvar' `useintvars' `covars' `other'

		// remove `useintvars' from `command':
		local command_new `cmdname' `depvar2' `usetrtvar' `covars' `other'
		local usevars `usetrtvar'
		
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
		// replace `usetrtvar' with new variable `refintvar'
		tempvar refintvar
		qui gen byte `refintvar' = `usetrtvar' * (`intvars'==`ref')		// `intvars' is the (single) original factor variable
		
		local command_new `cmdname' `depvar2' `refintvar' `useintvars' `covars' `other'
		local usevars `refintvar' `useintvars'
	}
	
	preserve
		qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`usevars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
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

		preserve
			qui mvmeta_make `command_new' if `touse', `offset_opts' `cmdopts' usevars(`usevars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
			qui mvmeta y S, fixed commonparm
			scalar `eff_pos' = e(b)[1,1]
		restore
		
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
		
		preserve
			qui mvmeta_make `command_new' if `touse', `offset_opts' `cmdopts' usevars(`usevars') /*names(_y _S)*/ by(`study') ppfix(check) `hard' clear
			qui mvmeta y S, fixed commonparm
			scalar `eff_neg' = e(b)[1,1]
		restore

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
		
		qui gen double _Int_xb = `ThetaHat'[1,1]*`usetrtvar'
		qui gen double _Int_stdp = .
		tokenize `useintvars'
		forvalues i = 1 / `k1' {
			qui replace _Int_xb = _Int_xb + `GammaHat'[1, `i']*``i''
		}
		
		mata: st_view(X=., ., "`usetrtvar'", "`touse'")
		mata: st_view(Z=., ., "`useintvars'", "`touse'")
		mata: VarThetaHat = st_matrix("`VarThetaHat'")
		mata: VarGammaHat = st_matrix("`VarGammaHat'")
		mata: GradZ = st_matrix("`Grad'")
		mata: GradZ = Z + J(rows(Z), 1, 1) * GradZ
		mata: Var_Xb = diag(VarThetaHat*X) + GradZ*VarGammaHat*GradZ'
		mata: st_store(st_viewobs(X), "_Int_stdp", sqrt(diagonal(Var_Xb)))

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




program define ExitErr 
	nois disp as err  "Failed to identify interaction variables"
	nois disp as err `"{it:command} should contain at least one term of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
	exit 198
end

