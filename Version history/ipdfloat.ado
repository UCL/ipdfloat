*** IPDFLOAT ***
* Routine for fitting IPD-based "floating subgroup" models
// by David Fisher, Peter Godolphin, Ian White
// (MRC Clinical Trials Unit at UCL)


* v0.1  David Fisher 17sep2021
// Very early version with a lot of functionalities hard-coded, mostly for test of principle

* v0.2  David Fisher 28sep2021
// Now generalisable enough to be applied (in testing) to real IPD meta-analysis data
// but still rather clunkily coded

* v0.3  David Fisher 17dec2021
// Now automatically identifies the treatment and interaction variables so long as they are supplied explicity as the first independent variable term
// Fully general estimation of predicted treatment effects for continuous interaction, using Mata

* v0.4  David Fisher 17jan2022
// Overhauled code that parses the command line, now much more general
//  - optional checks that compare direct use of "xi:command" with command built up from xi-derived covariates
// Checks version of -metan- (requires version 4.05  David Fisher  29nov2021;  N.B. no longer uses -ipdmetan- )
// Many functionalities from -metafloat- v0.9 copied across, e.g. random-effects, saved dataset

*! v0.5  David Fisher 03feb2022
// Various improvements, bug fixes, and tying up of loose ends
// Code for saved dataset moved to new subroutine "SavingClear"
//  - still to do: "design", "trend"


* Syntax:
*   ipdfloat, study(studyname) [options] : cmdname [depvar] i.trtvar*[i.]intvar [adjustment vars] [if] [in] [wgt] , [cmd_opts]

* where:
*  [depvar] = dependent variable (ONLY optional IF it's e.g. a survival model; otherwise NOT optional)
*  trtvar = treatment variable, taking on precisely two values (e.g. 0, 1)
*  intvar = covariate to be interacted with treatment
*    - if factor, write i.trt*i.intvar and only a single intvar is allowed
*    - if continuous, write i.trt*intvar1 [i.trt*intvar2 ...] and multiple intvars are allowed (e.g. to specify a FP or spline)

* Note: -ipdfloat- uses old-style "xi"-based factor notation, although "xi" itself should NOT be explicitly written within the command.
* (ultimately we may try to redesign this command so it handles modern factor-variable notation, but no immediate need)

* Reference level (for factor variable) or reference value (for first continuous variable) can be specified in advance using "char intvar[omit] #"; see -help char-
* Otherwise, -ipdfloat- will automatically use the lowest observed factor level, or the sample mean of (first) continuous covariate


program define ipdfloat, eclass
	
	version 11.0
	local version : di "version " string(_caller()) ":"	
	
	// Check that -metan- v4.05+ is installed
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
	// NOTE: if random-effects, also requires -mvmeta- version 3.5.1 # Ian White # 21dec2021
	// but currently no easy way to test for this
	
	_on_colon_parse `0'
	local 0       `"`s(before)'"'
	local command `"`s(after)'"'
	
	syntax, STUDY(varname) ///
		[ TRTVAR(passthru) INTVARs(passthru) ///
		  SHOWmodels HARD noSEPGroups AUGVARiance(real 1e5) EPS(real 1e-7) FORWards BACKWards Level(cilevel) SAVING(string) CLEAR ///
		  PRINT(string) * ]		// -mvmeta- options

	_get_eformopts, eformopts(`options') allowed(__all__) soptions
	local eformopt eform(`"`s(str)'"')
	local soptions `"`s(options)'"'
	
	`version' _prefix_command ipdfloat, level(`level') : `command'
	local cmdname `"`s(cmdname)'"'
	local cmdopts `"`s(options)'"'
	local rest `"`s(anything)'"'	
	if substr(`"`cmdname'"', 1, 2)!="st" {		// survival commands don't have a dependent variable
		gettoken depvar rest : rest
	}
	
	tempvar touse
	mark `touse' `s(if)' `s(in)'
	markout `touse' `study', strok
	
	
	*** Find treatment and interaction effect(s)  ( = contrasts if factor-variable)
	// First, identify `trtvar'
	// then search `command' for terms beginning "i.`trtvar'* ... "
	// and then finally run -xi- on the entire `command'
	// use the results to identify correct variables to use later, with -mvmeta-
	// Note: all variables must be numeric
	ProcessTerms `rest' if `touse', `trtvar' `intvars' `options' cmdname(`cmdname')
	if `"`checkexit'"'!=`""' exit	// for debugging
	
	local trtvar  `s(trtvar)'		// xi, prefix(_J) i.trtvar
	local intvars `s(intvars)'		// xi, prefix(_J) i.trtvar*[i.]intvar1 [i.trtvar*intvar2 ... ]  = _J[trt]X[int]_1_`i' if factor;  else = _J[trt]X[int]_1
	local covars  `s(covars)'		// main effect(s) of intvar(s)
	local adjvars `s(adjvars)'		// other adjustment variables (with xi applied if appropriate) not part of the main treatment-covariate interaction
	local facvar  `s(facvar)'		// if -factor- , this contains the original factor-variable
	local ref     `s(ref)'			// reference value and/or category
	
	local k1 : word count `intvars'
	local k = `k1' + 1
	
	// sreturn clear	// in advance of running s-class subroutines TransMatrix and CovStruct
	
	// Generate transformation matrix `T' based on choice of reference subgroup (= identity matrix if ref = 1st subgroup)
	// also define trend vector `d' if appropriate [NOT IMPLEMENTED YET IN -ipdfloat-]
	// [ subroutine taken directly from -metafloat- ]
	tempname T d
	if `"`facvar'"'!=`""' {
		TransMatrix `T' `d', subgroup(`facvar') `soptions'
		local soptions `"`s(options)'"'
		assert `k'==`s(k)'
		assert `ref'==`s(ref)'

		qui levelsof `facvar', local(sglist)
		local sgvarlab : variable label `facvar'
		local studylab : variable label `study'		// also save `study' variable label
		
		// save facvar value label
		// unlike in "metafloat", we need to use "label save" because we will manipulate the data before saving
		local sglab : value label `facvar'
		if `"`sglab'"'!=`""' {
			tempfile sglabfile
			qui label save `sglab' using `sglabfile'
		}
	}
	else matrix define `T' = I(`k')
	
	// Setup covariance structures
	// [ subroutine taken directly from -metafloat- ]
	tempname sgmat propU
	CovStruct `sgmat' `propU', k(`k') trans(`T') /*facvar(`facvar') `sepgroups'*/ `soptions'
	local structure  `s(structure)'		// must be one of: fixed, unstructured, exchangeable, randombeta, wscorrzero
	local sigmagamma `s(sigmagamma)'
	local sigmabeta  `s(sigmabeta)'
	local naive      `s(naive)'			// "naive" = undocumented option; don't constrain SigmaBeta and SigmaGamma to be linked
	local soptions `"`s(options)'"'

	if `"`facvar'"'==`""' | `"`sepgroups'"'!=`""' { 	// if non-factor (i.e. continuous); or factor without by-subgroup covariance structure
		if `"`naive'"'!=`""' {
			nois disp as err `"Option {bf:naive} is only valid with categorical covariates"'
			exit 198
		}
		if `"`structure'"'!=`"fixed"' {
			//local sigmabeta bscov(unstructured)
			//local naive naive
		}
	}
	
	/*	// TREND NOT IMPLEMENTED YET IN -ipdfloat-
	if `"`forcetrend'"'!=`""' {
		cap assert inlist(`"`structure'"', "fixed", "randombeta") | `"`s(default)'"'!=`""'		// either "fixed" or "randombeta" specified, or no structure was specified
		if _rc {
			nois disp as err `"cannot specify {bf:bscovariance(`structure')} with {bf:forcetrend}; random-effects structure is implied by use of trend"'
			exit 184
		}
	}
	*/
	
	local end = 0
	local r : list posof `"`ref'"' in sglist
	local usevars `trtvar' `intvars'
	forvalues i = 1/`k1' {
		local intvari : word `i' of `intvars'
	
		if `"`facvar'"'==`""' {		// continuous covariate
			local useJyvars `useJyvars' y`intvari'
		}
		else {		// factor variable
			local covari  : word `i' of `covars'

			// extract stem and form Ibase (only do this once)
			// from xi.ado:
			// local name = substr("`pre'`base'",1,11) + "_"
			// and if this name already exists, move on to [a-z][A-Z] in turn.
			// therefore, end of string is (_|[a-z]|[A-Z])[1-9]+				
			if !`end' {
				local len = length(`"`covari'"')
				local j = 1
				local el = substr(`"`covari'"', -`j', `j')
				cap confirm integer number `el'
				while !_rc {
					local el0 : copy local el
					local ++j
					local el = substr(`"`covari'"', -`j', `j')
					cap confirm integer number `el'
				}
				local n2 = `len' - `j' + 1
				local stem = substr(`"`covari'"', 1, `n2')				
				local Ibase = subinstr(`"`stem'`ref'"', `"_J"', `"_I"', 1)	// finally, switch prefix _J for _I
				local end = 1
			}
			
			local varlab : variable label `covari'
			tokenize `"`varlab'"', parse("=")
			args sgname eq sgval
			cap {
				assert `"`sgname'"'==`"`facvar'"'
				assert `"`eq'"'==`"=="'
				confirm number `sgval'
				assert `"`4'"'==`""'
			}
			if _rc {
				nois disp as err `"Something has gone wrong involving -xi-"'
				exit 198
			}
			
			// insert baseline expression `Ibase' into macro `Iyvars' at appropriate point
			// (for matrix labelling prior to displaying of results)
			if `i'==`r' local Iyvars `Iyvars' y`Ibase'
			local Iyvar = subinstr(`"y`covari'"', `"_J"', `"_I"', 1)
			local Iyvars `Iyvars' `Iyvar'
			local Jyvars `Jyvars' y`covari'

			if `sgval'!=`ref' {
				local useJyvars `useJyvars' y`intvari'
			}
		}
	}
	if `"`facvar'"'!=`""' & `r'==`k' {
		local Iyvars `Iyvars' y`Ibase'
	}	
		
	
	** Generate -mvmeta- dataset containing contrasts plus constant (reference subgroup)	
	if "`showmodels'"=="" local quietly quietly
	if "`print'"=="" local print bscov
	
	preserve
		qui mvmeta_make `cmdname' `depvar' `trtvar' `intvars' `covars' `adjvars' if `touse', `cmdopts' ///
			usevars(`trtvar' `intvars') by(`study') ppfix(check) `hard' `collapse_opt' nodetails clear

		local smiss = 0
		foreach s of varlist S_* {
			qui count if `s'>=`augvariance'
			if r(N) local smiss = 1
			qui replace `s' = `augvariance' if `s'>=`augvariance'		// REVISIT (or at least tell user what is happening and why)
		}
		if `smiss' {
			nois disp `"{error}missing covariances detected; replaced with augmented value `augvariance'"'
		}
		
		** STEP 1: ESTIMATE INTERACTION (MEAN AND BS VARIANCE)
		nois disp as text _n "Test of interaction(s):"		
		`quietly' mvmeta y S, vars(`useJyvars') print(`print') `sigmagamma' `soptions' `eformopt'
		nois test `useJyvars'
		local chi2_int = r(chi2)
		local df_int = r(df)
		local p_int = r(p)

		tempname GammaHat VarGammaHat SigmaGamma Chol
		matrix define `VarGammaHat' = e(V)[1..`k1', 1..`k1']
		matrix define `GammaHat'    = e(b)[1,       1..`k1']
		matrix define `SigmaGamma'  = e(Sigma)
		
		// extract cholesky matrix of Sigma
		if `"`structure'"'=="unstructured" {
			// special case: if `unstructured', -mvmeta- outputs the elements of `Chol' as part of e(b)
			// so don't bother packing them up into a matrix; we're only going to unpack them again later		
			matrix define `Chol' = e(b)[1, `k'...]
		}
		else if inlist(`"`structure'"', "exchangeable", "wscorrzero") {
			local tau = e(b)[1, `k']
			cholesky2 `SigmaGamma' `Chol'
			matrix define `Chol' = sign(`tau') * `Chol'
		}
		else if `"`structure'"'=="randombeta" {
			matrix define `Chol' = J(`k1', `k1', 0)
		}
		
		local rownames_eb : rownames e(b)
		matrix rownames `GammaHat' = `: word 1 of `rownames_eb''
		matrix colnames `GammaHat' = `Jyvars'
		matrix coleq `GammaHat' = :
		matrix coleq `VarGammaHat' = :
		matrix roweq `VarGammaHat' = :

		if `"`saving'"'!=`""' | `"`clear'"'!=`""' {

			// From the IPD, we have formed a dataset of "reference subgroup" + "contrasts"
			// However, for the saved data, we need *subgroups* rather than contrasts
			// Therefore, form non-reference subgroups "manually" using available y* and S* variables
			tokenize `intvars'

			rename y`trtvar' _y_ES`ref'
			qui gen double _y_seES`ref' = sqrt(S`trtvar'`trtvar')
			
			local i = 0
			foreach j of local sglist {
				if `j'!=`ref' {
					local ++i					
					
					qui gen double _y_ES`j' = _y_ES`ref' + y``i''
					qui gen double _y_seES`j' = sqrt(S`trtvar'`trtvar' + S``i''``i'' + 2*S`trtvar'``i'')
					
					rename y``i'' _yInt_ES`j'
					qui gen double _yInt_seES`j' = sqrt(S``i''``i'')					
				}
			}
			local extra _y_ES* _y_seES*
			local extra_reshape _y_ES _y_seES
			
			keep `study' _yInt_ES* _yInt_seES* `extra'
			qui reshape long _yInt_ES _yInt_seES `extra_reshape', i(`study') j(`facvar')
			
			tempfile main_effects
			qui save `main_effects'
		}

	restore
	
	
	** STEP 2: ESTIMATE TREATMENT EFFECT IN REF GROUP (MEAN, BS VARIANCE AND BS CORRELATION WITH INTERACTIONS)
	// subtract within-trial interactions from observed interaction data
	// set up equations that make shifted interactions have mean zero
	
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
	if inlist(`"`r(cmd)'"', "regress", "xtreg", "mixed") {		// commands with continuous outcomes, from which offset can be directly subtracted
		tempvar depvar_offset
		qui gen double `depvar_offset' = `depvar' - `offset'
		local depvar2 : copy local depvar_offset
	}
	else {
		local offset_opt offset(`offset')
		local depvar2 : copy local depvar
	}

	local eqlist
	if `"`facvar'"'==`""' | `"`sepgroups'"'!=`""' {			// non-factor (i.e. continuous); or factor without by-subgroup covariance structure
		// remove `intvars' from `command':
		local command_new `cmdname' `depvar2' `trtvar' `covars' `adjvars'
		local usevars `trtvar'
	}
		
	else {
		local command_new `cmdname' `depvar2' `trtvar' `intvars' `covars' `adjvars'
		local usevars `intvars' `trtvar'	// put `trtvar' last
				
		local i = 0
		tokenize `useJyvars'
		while `"`1'"'!=`""' {
			local ++i
			local eqlist `eqlist' `1':_Zero `designvars',
			macro shift
		}
		local eqlist equations(`eqlist' y`trtvar':_One `designvars') nocons
	}
	
	// set up constraints on the BS variance of the interactions
	// using the Cholesky matrix from STEP 1
	if `"`structure'"'!="fixed" & `"`naive'"'==`""' {
		SetupConstraints `Chol' `propU', k(`k') structure(`structure') `forcetrend'
		local constr_list `s(constr_list)'
		
		if `"`constr_list'"'!=`""' {
			if `"`listconstraints'"'!=`""' {
				nois disp as text _n "Constraints on elements of Cholesky factor for SigmaU:"
				nois constraint dir `constr_list'
			}
			local constr_opt constraints(`constr_list')
		}
	}
	
	preserve
		qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`usevars') by(`study') ppfix(check) `hard' clear
		foreach s of varlist S_* {
			qui replace `s' = `augvariance' if `s'>=`augvariance'		// see above
		}
		qui gen byte _Zero = 0
		qui gen byte _One  = 1
		
		// fit model: make sure we are using -mvmeta- version 3.5.1 # Ian White # 21dec2021 (or later)
		if "`showmodels'"!="" nois disp as text _n `"Estimation of floating subgroup for reference category:"'		
		`quietly' mvmeta y S, print(`print') commonparm `eqlist' `constr_opt' `sigmabeta' `soptions'

		tempname ThetaHat VarThetaHat SigmaU
		matrix define `ThetaHat' = e(b)[1,1]
		matrix define `SigmaU'   = e(Sigma)
		
		// correct VarThetaHat for uncertainty in the heterogeneity matrix if requested
		if `"`uncertainv'"'==`""' {
			matrix define `VarThetaHat' = e(V)[1,1]
		}
		else {
			matrix define `VarThetaHat' = e(V)
			matrix define `VarThetaHat' = invsym(`VarThetaHat')
			matrix define `VarThetaHat' = `VarThetaHat'[1,1]
			matrix define `VarThetaHat' = invsym(`VarThetaHat')
		}
	restore
	

	*** Obtain "delta method" gradients, in order to form the final variance-covariance matrix
	tempname Grad h eff_pos eff_neg
	matrix define `Grad' = J(1, `k1', .)

	tokenize `intvars'
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
		else local offset_opt offset(`offset_h')

		preserve
			qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`usevars') by(`study') ppfix(check) `hard' clear	
			foreach s of varlist S_* {
				qui replace `s' = `augvariance' if `s'>=`augvariance'		// see above
			}
			qui gen byte _Zero = 0
			qui gen byte _One  = 1			
			qui mvmeta y S, commonparm `eqlist' `constr_opt' `sigmabeta' `soptions'	// REVISIT -- do constraints still work with updated offset??
			scalar `eff_pos' = e(b)[1,1]
		restore
		
		if `"`forwards'`backwards'"'!=`""' {
			matrix `Grad'[1, `i'] = `sign'(`eff_pos' - `ThetaHat'[1,1])/`h'
			continue
		}

		// update offset
		qui replace `offset_h' = `offset' - `h'*``i''
		unabcmd `cmdname'
		if inlist(`"`r(cmd)'"', "regress", "xtreg", "mixed") {
			qui replace `depvar_offset' = `depvar' - `offset_h'
		}
		
		preserve
			qui mvmeta_make `command_new' if `touse', `offset_opt' `cmdopts' usevars(`usevars') by(`study') ppfix(check) `hard' clear
			foreach s of varlist S_* {
				qui replace `s' = `augvariance' if `s'>=`augvariance'		// see above
			}
			qui gen byte _Zero = 0
			qui gen byte _One  = 1
			qui mvmeta y S, commonparm `eqlist' `constr_opt' `sigmabeta' `soptions'	// REVISIT -- do constraints still work with updated offset??
			scalar `eff_neg' = e(b)[1,1]
		restore

		matrix `Grad'[1, `i'] = (`eff_pos' - `eff_neg')/(2*`h')
	}
		

	*** Assemble matrices
	if "`facvar'"!="" {
		tempname ones Z A BetaHat VarBetaHat
		matrix define `ones' = J(`k', 1, 1)
		matrix define `Z' = J(`k', `k1', 0)
		matrix define `Z'[2, 1] = I(`k1')

		matrix define `BetaHat' = `ones'*`ThetaHat' + `Z'*`GammaHat''
		matrix define `A' = `ones'*`Grad' + `Z'
		matrix define `VarBetaHat' = `VarThetaHat'*J(`k', `k', 1) + `A'*`VarGammaHat'*`A''
		
		// reverse transformation
		matrix define `BetaHat'    = `T''*`BetaHat'
		matrix define `VarBetaHat' = `T''*`VarBetaHat'*`T'
		
		// collect and display under "ereturn"
		matrix define `BetaHat' = `BetaHat''
		matrix coleq `BetaHat' = :
		matname `BetaHat'    `Iyvars', columns(.) explicit
		matname `VarBetaHat' `Iyvars', explicit
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
		mata: A = st_matrix("`Grad'")
		mata: A = J(rows(Z), 1, 1)*A + Z
		mata: Var_Xb = diag(VarThetaHat*X) + A*VarGammaHat*A'
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
		
		// generate matrices containing pooled estimates and variances
		tempname BetaTable GammaTable
		matrix define `GammaTable' = `GammaHat'', vecdiag(`VarGammaHat')'
		matrix define `BetaTable' = `BetaHat'', vecdiag(`VarBetaHat')'
		
		preserve
		use `main_effects', clear
		
		// re-instate original variable labels and value labels
		if `"`sglabfile'"'!=`""' {
			qui do `sglabfile'
			label values `facvar' `sglab'
		}
		// label variable `facvar' `"`sgvarlab'"'
		// label variable `study' `"`studylab'"'
		label variable `study' `"`sgvarlab'"'		// swap over, so they display correctly after -metan-
		label variable `facvar' `"`studylab'"'
		
		SavingClear `GammaTable' `BetaTable', saving(`saving') ///
			study(`study') subgroup(`facvar') augvariance(`augvariance')
		
		if `"`clear'"'!=`""' {
			restore, not
		}		
	}

	
	*** Finally, present the results
	if "`facvar'"!="" {

		ereturn post `BetaHat' `VarBetaHat', depname(Subgroup)

		ereturn matrix GammaHat    = `GammaHat'
		ereturn matrix VarGammaHat = `VarGammaHat'
		ereturn matrix SigmaGamma  = `SigmaGamma'
		
		tempname SigmaBeta
		if `"`sepgroups'"'==`""' {
			tempname L
			matrix define `L' = `Z', `ones'
			matrix define `SigmaBeta' = `L'*`SigmaU'*`L''
		}
		else matrix define `SigmaBeta' = `SigmaU'
		matrix define `SigmaBeta' = `T''*`SigmaBeta'*`T'
		ereturn matrix SigmaBeta  = `SigmaBeta'

		ereturn matrix Grad        = `Grad'
		ereturn scalar ThetaHat    = `ThetaHat'[1,1]
		ereturn scalar VarThetaHat = `VarThetaHat'[1,1]
		// ereturn scalar n = `n'		// number of studies
		ereturn scalar k = `k'		// number of subgroups

		ereturn scalar chi2_int = `chi2_int'
		ereturn scalar df_int = `df_int'
		ereturn scalar p_int = `p_int'
		
		ereturn local bscov_Gamma `"`sigmagamma'"' 
		ereturn local bscov_Beta `"`sigmabeta'"' 
		
		nois disp as text _n "Floating subgroups:"
		ereturn display, `eformopt'

		if `"`constr_list'"'!=`""' & `"`structure'"'!=`"user"' {	// only drop internally-defined constraints
			constraint drop `constr_list'
		}
	}
	
	ereturn local command     `"`command'"'
	ereturn local command_new `"`command_new'"'
	
end




program define ProcessTerms, sclass

	syntax anything [if] [in], [TRTVAR(varname numeric) INTVARs(varlist numeric) ///
		check cmdname(string) /* undocumented; for debugging */ * /* -mvmeta- options; ignore here*/ ]
	
	marksample touse
	
	// Assume first element of `command' is of the form: i.`trtvar'*`i'`intvar'
	if `"`trtvar'"'==`""' {
		tokenize `"`anything'"', parse(" .*")
		unab trtvar : `3', max(1)
		capture {
			assert `"`1'`2'"'==`"i."'
			assert `"`4'"'==`"*"'
		}
		if _rc {
			nois disp as err  "Could not identify treatment variable. If option {bf:trtvar()} is not specified,"
			nois disp as err `" then first element of {it:command} should be of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
			exit 198
		}
		confirm numeric variable `trtvar'
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
	
	// Now search `anything' for terms beginning "i.trt* ... "
	local specintvars : copy local intvars
	local remain : copy local intvars
	gettoken el rest : anything
	local k = 1						// for treatment
	while `"`el'"'!=`""' {
		tokenize `"`el'"', parse(".*|")
		
		if lower(`"`1'`2'"')==`"i."' {
			unab newvar : `3'
			local facvars `facvars' `el'
			if `"`factor'"'!=`""' {
				gettoken el rest : rest
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
					if inlist(`"`intvars'"', `"`newvar'"', `""') {
						local remain : list remain - newvar
						local intvars : copy local newvar
						local factor factor
						gettoken el rest : rest
						continue	// *not* "continue, break" because there may be `adjvars' to collect
					}
				}
				
				// continuous
				else {
					unab newvar : `5'
					if `"`intvars'"'==`""' | `: list newvar in intvars' {
						local remain : list remain - newvar
						local covars `covars' `newvar'
						local ++k
					}
				}
			}
		}
		else local adjvars `adjvars' `el'
		gettoken el rest : rest
	}
	
	if `"`remain'"'!=`""' {
		nois disp as err `"Failed to identify interaction(s) between treatment variable {bf:`trtvar'} and covariate(s) {bf:`remain'}"'
		nois disp as err `"Please check {it:command} and {bf:intvars()} option, and try again"'
		exit 198
	}
	confirm numeric variable `intvars' `covars'
	
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
		if `"`specintvars'"'!=`""' {
			local problem : list specintvars - covars
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
	
	// At this point (note, all variables are unchanged from data in memory; -xi- not yet applied):
	// `trtvar' identifies treatment variable
	// `intvars' identifies variables interacted *with treatment* AND in intvars() if supplied (*not* any other interactions)
	// `covars' == `intvars' if continuous; otherwise blank.
	
	
	** Now run -xi- on the supplied command
	// and use the results to check again for consistency and for correct variables to use
	qui xi, prefix(_J) `anything'
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
			cap assert `"`rest'"'==`""'
			if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
			
			// parse first part (treatment)
			local open
			local close
			tokenize `"`trtpart'"', parse("(=)")
			if `"`1'"'!=`"("' {
				args trtname eq trtval					// format:  trtname==trtval [& intname==intval]  ==> factor `intname' (if present)
				cap {
					assert `"`factor'"'!=`""' | `"`symbol'"'==`""' if `"`trtname'"'==`"`trtvar'"'
					assert `"`2'"'==`"=="'
					assert `"`4'"'==`""'
				}
			}
			else {
				cap {
					assert `"`1'`3'`5'"'==`"(==)"'
					assert `"`6'"'==`""'
				}
				args open trtname eq trtval close		// format:  (trtname==trtval)*intname  ==> continuous `intname'
			}
			if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
			
			// parse second part (interaction)
			if `"`symbol'"'!=`""' {
				tokenize `"`intpart'"', parse("=")
				args intname eq intval
				cap {
					assert (`"`symbol'"'==`"*"' & `"`open'"'!=`""') | (`"`symbol'"'==`"&"' & `"`open'"'==`""')
					assert `"`intname'"'!=`""' & `"`eq'"'==`""'   if `"`symbol'"'==`"*"'		// format:  (trtname==trtval)*intname  ==> continuous `intname'
					assert `"`intval'"'!=`""'  & `"`eq'"'==`"=="' if `"`symbol'"'==`"&"'		// format:  trtname==trtval [& intname==intval]  ==> factor `intname' (if present)
					assert `"`4'"'==`""'
				}
				if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
			}

			// if interaction with *treatment* (or treatment alone) ...
			if `"`trtname'"'==`"`trtvar'"' {
				cap assert `trtval' != `trtref'
				if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')

				// if treatment variable alone (i.e. *not* interaction), then continue to process next element
				if `"`symbol'"'==`""' {
					local usetrtvar : copy local eli
					continue
				}

				// if a treatment-covariate interaction, check whether `intname' was specified by intvars()
				// if so:			
				if (`"`specintvars'"'==`""' | `: list intname in specintvars') {
					if `"`factor'"'!=`""' {
						// if we have a factor-variable interaction defined, cannot have any futher interactions with treatment
						if `"`specintvars'"'==`""' {
							if `"`open'"'!=`""' {
								nois disp as err "Cannot interact treatment with both a factor and a continuous variable"
								qui ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
							}
							if `"`intname'"'!=`"`intvars'"' {
								nois disp as err "Cannot interact treatment with more than one factor variable"
								qui ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
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
					if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
					
					local useintvars `useintvars' `eli'
					local intnames `intnames' `intname'
					continue
				}
			}		// end if `"`trtname'"'==`"`trtvar'"'
				
			// otherwise, must be one of:
			//  - a main effect of a factor-variable covariate;
			//  - a treatment-covariate interaction not specified by intvars(), so not relevant to the procedure;
			//  - a non-treatment interaction, which is being adjusted-for but is not relevant to the procedure.
			cap assert `"`intname'"'!=`"`trtvar'"'
			if _rc {
				nois disp as err "Cannot specify an interaction term where the treatment variable is in second position"
				qui ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
			}
			
			// `intvars' in first position: quick check for consistency of ref value [IS THIS NEEDED?]
			if `"`trtname'"'==`"`intvars'"' {
				cap assert `trtval'!=`ref'
				if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
			}
			
			// if interaction with continuous covariate, need to save main-effect of this covariate as well as the interaction
			if `"`intpart'"'==`""' {
				if `: list trtname in intvars' & !`: list trtname in covars' local covars `covars' `eli'
				else local adjfacvars `adjfacvars' `eli'
			}
			else {
				if `"`symbol'"'==`"*"' local adjfacvars `adjfacvars' `intname'
				local adjfacvars `adjfacvars' `eli'
			}
		}		// end if `"`el'"'==`"_J"'
	}		// end foreach el of local OldVarsPrefix 

	local adjfacvars : list uniq adjfacvars

	cap {
		confirm numeric variable `trtvar' `usetrtvar'
		assert `: word count `trtvar'' == 1
		assert `: word count `usetrtvar'' == 1

		confirm numeric variable `intvars' `useintvars'	`covars'
		assert `: word count `useintvars'' == `k1'
		assert `: word count `covars'' == `k1'
		
		if `"`factor'"'==`""' {
			assert `: list intvars === covars'
			assert `: list intvars === intnames'		// this is trivial; remove??
			assert `: word count `intvars''==`: word count `useintvars''
		}
	}
	if _rc nois ExitErr `VarsToDrop', oldpre(`OldVarsPrefix') oldvars(`OldVarsToDrop')
	
	sreturn local trtvar `usetrtvar'		// xi, prefix(_J) i.trtvar
	sreturn local intvars `useintvars'		// xi, prefix(_J) i.trtvar*[i.]intvar1 [i.trtvar*intvar2 ... ]  = _J[trt]X[int]_1_`i' if factor;  else = _J[trt]X[int]_1
	if `"`factor'"'!=`""' {
		sreturn local facvar `intvars'		// if -factor- , this contains the original factor-variable
	}
	sreturn local covars `covars'					// main effect(s) of intvar(s)
	sreturn local adjvars `adjfacvars' `adjvars'	// other adjustment variables *not* part of the main treatment-covariate interaction(s)
													//  (including those generated by -xi- where appropriate)
	sreturn local ref `ref'							// reference value and/or category	


	
	************************
	* TEMP -- extra checks *
	************************
	
	* Fit original command and new command, and compare statistics to check they are the same model
	if `"`check'"'==`""' exit
	c_local checkexit checkexit

	disp `"trtvar: `trtvar'"'
	disp `"(use)intvars: `useintvars'"'
	disp `"int(fac)vars: `intvars'"'
	disp `"facvars: `facvars'"'
	disp `"covars: `covars'"'
	disp `"adjfacvars: `adjfacvars'"'
	disp `"adjvars: `adjvars'"'
	disp `"usevars: `usetrtvar' `useintvars'"'	
	disp `"k = `k'"'
	
	disp _n `"-----------"'
	
	disp trim(itrim(`"`cmdname' `anything'"'))
	disp trim(itrim(`"`cmdname' `depvar' `facvars' `adjvars'"'))
	local cmd2 `"`depvar' `facvars' `adjvars'"'
	assert `: list anything === cmd2'
	local cmd2 `"`cmdname' `cmd2'"'
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

	qui xi, prefix(_J): `cmdname' `anything'
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

	disp `"`: list sort newyvars1'"'
	disp `"`: list sort newyvars2'"'
	assert `: list newyvars1 === newyvars2'
	disp `"`eN_1' == `eN_2'"'
	assert `eN_1'==`eN_2'
	disp `"`eN_fail_1' == `eN_fail_2'"'
	assert `eN_fail_1'==`eN_fail_2'
	disp `"`ell_0_1' == `ell_0_2'"'
	assert float(`ell_0_1')==float(`ell_0_2')
	disp `"`ell_1' == `ell_2'"'
	assert float(`ell_1')==float(`ell_2')

end

	
	
// process covariance structures
// [ subroutine taken directly from -metafloat- ]
program define CovStruct, sclass
	syntax anything, k(integer) trans(name) [ /*facvar(passthru) sepgroups*/ BSCOVariance(string) naive ///
		fixed UNStructured EXCHangeable RANDOMBeta WSCORRZero WSC0 * ]
	
	local soptions : copy local options
	if `"`wsc0'"'!=`""' local wscorrzero wscorrzero
	opts_exclusive `"`fixed' `unstructured' `exchangeable' `randombeta' `wscorrzero'"' `""' 184
	local structure `fixed'`unstructured'`exchangeable'`randombeta'`wscorrzero'

	if `"`structure'"'!=`""' {
		if `"`bscovariance'"'!=`""' {
			nois disp as err `"Cannot specify both {bf:`structure'} and {bf:bscovariance()}"'
			exit 198
		}
	}
	else {
		local 0 `bscovariance'
		syntax [ , FIXED UNStructured EXCHangeable RANDOMBeta WSCORRZero WSC0 ]
		if `"`wsc0'"'!=`""' local wscorrzero wscorrzero
		opts_exclusive `"`fixed' `unstructured' `exchangeable' `randombeta' `wscorrzero'"' `""' 184
		local structure `fixed'`unstructured'`exchangeable'`randombeta'`wscorrzero'
	}
	
	// parse pre-defined tempname matrices
	tokenize `anything'
	args sgmat propU
	local T : copy local trans
	
	if `"`structure'"'=="fixed" {
		local sigmagamma fixed
		local sigmabeta fixed
	}
	else if `"`structure'"'=="randombeta" {
		local sigmagamma fixed
		matrix define `propU' = J(`k', `k', 0)
		matrix define `propU'[`k', `k'] = 1
		matrix define `propU' = `T''*`propU'*`T'
		local sigmabeta bscov(proportional `propU')
	}		
	else {
		local sigmagamma bscov(unstructured)
		local sigmabeta bscov(unstructured)

		if inlist(`"`structure'"', "exchangeable", "wscorrzero") {
			local sigmagamma bscov(exchangeable .5)
			
			/*
			if `"`sepgroups'"'!=`""' {	// revisit whether == or !=
				local sigmabeta bscov(exchangeable)			// using "NEW" mvmeta version 3.5.1 # Ian White # 21dec2021 -- REVISIT!!
			}
			*/
			
			if `"`structure'"'=="wscorrzero" {
				local sigmabeta bscov(equals `propU')	// matrix yet to be defined -- depends on SigmaGamma
				// local sigmabeta bscov(exchangeable 0)		// "naive" form: proportional to identity, zero covariances
			}
		}
		else if `"`structure'"'==`""' {
			sreturn local default default
			local structure unstructured
		}
	}

	if `"`naive'"'!=`""' & `"`structure'"'!=`"unstructured"' {
		nois disp as err `"Option {bf:naive} is only valid with {bf:unstructured} covariance structure"'
		exit 198
	}
	
	sreturn local structure    `structure'
	sreturn local sigmagamma `"`sigmagamma'"'
	sreturn local sigmabeta  `"`sigmabeta'"'
	sreturn local naive      `"`naive'"'	
	sreturn local options    `"`soptions'"'	

end


// Generate transformation matrix based on choice of reference subgroup
// ...and set up trend if appropriate
// [ subroutine taken directly from -metafloat- ]
program define TransMatrix, sclass

	syntax namelist, subgroup(varname) [ noTRend * ]		// "notrend" not yet implemented in -ipdfloat-
	tokenize `namelist'
	args T d
	
	qui tab `subgroup'
	local k = r(r)
	local k1 = `k' - 1
	qui levelsof `subgroup', local(sglist)		

	if `"`: char `subgroup'[omit]'"'==`""' {
		local r = 1
		gettoken ref sgnoref : sglist
		matrix define `T' = I(`k')
	}
	else {
		local ref : char `subgroup'[omit]
		local r : list posof `"`ref'"' in sglist
		cap assert `r' > 0
		if _rc {
			nois disp as err `"Reference value `ref' not taken on by variable {bf:`subgroup'} within the analysis sample"'
			exit 198
		}
		if `r'==1 matrix define `T' = I(`k')
		else {
			local kr = `k' - `r'
			local kr1 = `kr' + 1
			matrix define `T' = J(1, `k', 0)
			matrix define `T'[1, `r'] = 1
			if `r'==`k' {
				matrix define `T' = `T' \ ( I(`k1'), J(`k1', 1, 0))
			}
			else {
				tempname Id
				matrix define `Id' = I(`k1')
				matrix define `T' = `T' \ ( `Id'[1..`k1', 1..`kr'], J(`k1', 1, 0), `Id'[1..`k1', `kr1'..`k1'] )
			}
		}
	}

	// Note:
	// `sglist' is an **ordered** (according to -levelsof- ) list of subgroup identifiers
	// `r' is the index w.r.t. `sglist'
	// `ref' is the relevant element of `sglist' (e.g. for labelling the matrix stripe elements)
	// (and, unless `string', `ref' will also be a numeric value, which may have a value label)
	
	// "trend" not yet implemented in -ipdfloat-
	/*
	* Setup test for trend
	if `"`trend'"'==`""' {
		forvalues i = 1 / `k' {
			matrix define `d' = nullmat(`d') \ `i'	// linear trend (could in theory allow more general forms for `d')
		}

		// Convert `d' into a matrix of contrasts
		// ... and apply transformation if reference is not first subgroup
		tempname M
		matrix define `M' = J(`k1', 1, -1), I(`k1')
		matrix define `d' = `M' * `T' * `d'
	}
	*/
	
	sreturn local k `"`k'"'
	sreturn local ref  `"`ref'"'	
	sreturn local options `"`options'"'
	
end


program define cholesky2
* produce an approximate cholesky decomposition, even if matrix is not positive definite
* (Note: taken from -mvmeta- by Ian White)
args M C
cap matrix `C' = cholesky(`M')
if _rc {
    local eps = trace(`M')/1000
    if `eps'<=0 {
        di as error "cholesky2: matrix `M' has non-positive trace"
        matrix list `M'
        exit 498
    }
    local p = rowsof(`M')
    cap matrix `C' = cholesky(`M'+`eps'*I(`p'))
    if _rc {
        di as error "cholesky2: can't decompose matrix `M'"
        matrix list `M'
        exit 498
    }
}
end


program define SetupConstraints, sclass

	syntax namelist(min=2 max=2), k(integer) structure(string) [forcetrend]
	tokenize `namelist'
	args Chol propU
	
	local k1 = `k' - 1
	if `"`structure'"'!="wscorrzero" {
		local c = 0
		forvalues i = 1 / `k1' {
			forvalues j = 1 / `i' {
				local ++c
				
				if `"`structure'"'=="unstructured" & `"`forcetrend'"'==`""' {
					local el = `Chol'[1, `c']	// special case: see above
				}
				else local el = `Chol'[`i', `j']
				
				constraint free
				constraint `r(free)' [chol`j'`i']_b[_cons] = `el'
				local constr_list `constr_list' `r(free)'
			}
		}
	}
	
	// if exchangeable, need to do a little more work
	if `"`structure'"'=="exchangeable" {
		local tau = `Chol'[1, 1]
		tempname Chol_k
		matrix define `Chol_k' = inv(`Chol') * J(`k1', 1, -.5*(`tau')^2)
		forvalues j = 1 / `k1' {
			constraint free
			constraint `r(free)' [chol`j'`k']_b[_cons] = `=`Chol_k'[`j', 1]'
			local constr_list `constr_list' `r(free)'
		}
	}
	
	// similarly if trend
	else if `"`forcetrend'"'!=`""' {
		forvalues j = 1 / `k1' {
			constraint free
			constraint `r(free)' [chol`j'`k']_b[_cons] = 0
			local constr_list `constr_list' `r(free)'
		}
	}
	
	// finally, wscorrzero is a special case: SigmaBeta is completely defined by SigmaGamma
	else if `"`structure'"'=="wscorrzero" {
		local tau = `Chol'[1, 1]
		
		tempname LinvT		// = Linv * Linv-transpose, where `L' is as defined in the main routine
		matrix define `LinvT' = I(`k') + J(`k', `k', 1)
		matrix define `LinvT'[1, `k'] = J(`k1', 1, -1) 
		matrix define `LinvT'[`k', 1] = J(1, `k1', -1) 
		matrix define `LinvT'[`k', `k'] = 1

		matrix define `propU' = ((`tau')^2)*`LinvT'
	}
	
	// constraint 4 y_Isubgroup_2=0
	// constraint 5 y_Isubgroup_3=0

	sreturn local constr_list `"`constr_list'"'
	
end


program define SavingClear

	syntax anything, study(varname) subgroup(varname) ///
		[ SAVING(string asis) AUGVARiance(real 1e5) ] 
	
	tokenize `anything'
	args GammaTable BetaTable
	
	// identify "abbreviated" form of `subgroup' used by -xi-
	// from xi.ado:
	// local name = substr("`pre'`base'",1,11) + "_"
	// and if this name already exists, move on to [a-z][A-Z] in turn.
	// therefore, end of string is (_|[a-z]|[A-Z])[1-9]+	
	local Iyvars : rownames `BetaTable'
	local Jyvars : rownames `GammaTable'
	local k  : word count `Iyvars'
	local k1 : word count `Jyvars'
	assert `k1' == `k' - 1

	local pre = 3		// length of prefix; here hard-coded as y_I or y_J
	local i = 1
	foreach name in `Iyvars' `Jyvars' {
		local len = length(`"`name'"')
		local j = 1
		local el = substr(`"`name'"', -`j', `j')
		cap confirm integer number `el'
		while !_rc {
			local el0 : copy local el
			local ++j
			local el = substr(`"`name'"', -`j', `j')
			cap confirm integer number `el'
		}
		local n1 = `pre' + 1
		local n2 = `len' - `pre' - `j'
		local sgab = substr(`"`name'"', `n1', `n2')
		cap assert `"`sgab'"'==substr(`"`subgroup'"', 1, `n2')
		if _rc {
			nois disp as err `"Something has gone wrong involving -xi-"'
			exit 198
		}
		if `i' <= `k' local sglist `sglist' `el0'
		else local sgnoref `sgnoref' `el0'
		local ++i
	}
	local ref : list sglist - sgnoref
	local r : list posof `"`ref'"' in sglist	
	
	// insert missing reference row into `GammaTable' so that its elements match with those of `BetaTable'
	if `r'==1 {
		matrix define `GammaTable' = (., .) \ `GammaTable'[1..`k1', 1..2]
	}
	else if `r'==`k' {
		matrix define `GammaTable' = `GammaTable'[1..`k1', 1..2] \ (., .)
	}
	else {
		local r1 = `r' - 1
		matrix define `GammaTable' = `GammaTable'[1..`r1', 1..2] \ (., .) \ `GammaTable'[`r'..`k1', 1..2]
	}
	matrix rownames `GammaTable' = `sgnoref'
	matrix colnames `GammaTable' = _yInt _yInt_S
	matrix rownames `BetaTable' = `sglist'
	matrix colnames `BetaTable' = _y _y_S
		
	// add in pooled estimates
	tempvar pooled obs
	qui gen byte `pooled' = 1
	local oldN = _N
	local newN = `oldN' + `k'
	qui set obs `newN'
	qui replace `pooled' = 0 if _n > `oldN'
	
	foreach x of local sglist {
		local ++oldN
		qui replace `subgroup' = `x' in `oldN'
	}
	
	cap confirm numeric variable `study'
	if !_rc {
		summ `study', meanonly
		qui replace `study' = `r(max)' + 1 if !`pooled'
	}
	else {
		qui replace `study' = "_Pooled" if !`pooled'
	}
	
	// -svmat- places results starting at the first observation
	// so need to re-sort so that new observations are at the top, then re-sort again afterwards
	qui gen `obs' = _n
	sort `pooled' `obs'
	
	svmat `BetaTable', names(col)
	qui replace _y_ES = _y if !`pooled'
	qui replace _y_seES = sqrt(_y_S) if !`pooled'
	drop _y _y_S
	
	svmat `GammaTable', names(col)
	qui replace _yInt_ES = _yInt if !`pooled'
	qui replace _yInt_seES = sqrt(_yInt_S) if !`pooled'
	drop _yInt _yInt_S

	qui replace `pooled' = 1 - `pooled'
	sort `pooled' `obs'
	drop `obs'
	
	// use -metan- to restructure
	qui replace _yInt_ES = . if _yInt_seES^2 >= `augvariance'
	qui replace    _y_ES = . if    _y_seES^2 >= `augvariance'
	
	qui ds `study' `subgroup' _y_ES _y_seES, not
	local lcolvars `"`r(varlist)'"'

	qui metan _y_ES _y_seES, nograph nohet nosubgroup nooverall keeporder ///
		study(`subgroup') by(`study') rcols(`lcolvars') ///
		clear prefix(_y)

	qui drop *_EFFECT
	qui rename _y_BY _BY
	qui rename _y_STUDY _STUDY
	qui rename _y_LABELS _LABELS
	
	summ _BY if `pooled'==1, meanonly
	qui replace `pooled' = (_BY==`r(min)')
	
	qui replace _BY = . if `pooled'
	qui replace _LABELS = "Pooled estimates" if _y_USE==0 & `pooled'
	qui drop in L
	
	qui order _y_USE, before(_y_ES)
	qui replace _y_USE = 5 if inlist(_y_USE, 1, 2) & `pooled'

	summ _y_WT if !`pooled', meanonly
	qui replace _y_WT = 100*_y_WT / r(sum) if !`pooled'
	summ _y_WT if `pooled', meanonly
	qui replace _y_WT = 100*_y_WT / r(sum) if `pooled'
	
	qui gen double _yInt_LCI = _yInt_ES - invnormal(.975)*_yInt_seES
	qui gen double _yInt_UCI = _yInt_ES + invnormal(.975)*_yInt_seES
	qui gen double _yInt_WT = 1/_yInt_seES^2
	
	summ _yInt_WT if !`pooled', meanonly
	qui replace _yInt_WT = 100*_yInt_WT / r(sum) if !`pooled' & _STUDY!=`ref'
	summ _yInt_WT if `pooled', meanonly
	qui replace _yInt_WT = 100*_yInt_WT / r(sum) if `pooled' & _STUDY!=`ref'
	label variable _yInt_WT "% Weight"
	format _yInt_WT %6.2f

	qui gen byte _yInt_USE = _y_USE
	qui replace _yInt_USE = 9 if _STUDY==`ref'
	qui replace _yInt_USE = 2 if _yInt_USE==1 & missing(_yInt_seES)
	
	order _yInt_USE _yInt_ES _yInt_seES _yInt_LCI _yInt_UCI _yInt_WT, after(_y_WT)
	qui compress	
	
	if `"`saving'"'!=`""' {
		_prefix_saving `saving'
		qui save `s(filename)', `s(replace)'
	}
	
end


program define ExitErr
	syntax [varlist], [oldpre(namelist) oldvars(namelist)]

	disp as err  "Failed to identify interaction variables"
	disp as err `"{it:command} should contain at least one term of the form {bf:i.}{it:trtvar}{bf:*}[{bf:i.}]{it:intvar}"'
	
	if `"`varlist'"'!=`""" {'
		qui drop `varlist'
		char define _dta[__xi__Vars__Prefix__] `OldVarsPrefix'
		char define _dta[__xi__Vars__To__Drop__] `OldVarsToDrop'	
	}
	
	exit 198
end
