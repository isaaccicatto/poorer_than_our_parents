* ============================================================================
* 01_POSITIVE — raw CSV -> RA life-cycles -> Stata_Cache/$TAG_pos.dta
* All parameters are READ FROM THE CSV (no hardcoded values).
* Region suffix is normalized to _x throughout.
* Run via 00_master (needs $TAG $REGION $CSV $YEAR0 $HORIZON $COHORTS).
* ============================================================================
if "$TAG" == "" {
    di as error "globals not set — run via 00_master.do"
    exit 198
}
local out "Plots/$TAG"
cap mkdir "`out'"

import delimited "$CSV", clear case(lower) asdouble

* ---- normalize region suffix -> _x -----------------------------------------
local regsuf = lower("$REGION")
rename *_`regsuf' *_x

sort time
gen year = $YEAR0 + time
assert year[1] == $YEAR0

* ---- parameters from the data (constant columns; bare or suffixed names) ---
capture program drop _getpar
program define _getpar, rclass
    args name
    capture confirm variable `name'
    if !_rc {
        assert `name' == `name'[1]
        return scalar value = `name'[1]
        exit
    }
    capture confirm variable `name'_x
    if !_rc {
        assert `name'_x == `name'_x[1]
        return scalar value = `name'_x[1]
        exit
    }
    di as error "parameter `name' not found (tried `name' and `name'_x)"
    error 111
end

foreach p in beta sigma nu omega {
    _getpar `p'
    local `p' = r(value)
}
di as text "params (from CSV): beta=`beta' sigma=`sigma' nu=`nu' omega=`omega'"

* Delta (reference; used in the V-based appendix checks)
gen Delta_y_x = mu_y_x^(`sigma'/(1 - `sigma'))
gen Delta_m_x = mu_m_x^(`sigma'/(1 - `sigma'))
gen Delta_o_x = mu_o_x^(`sigma'/(1 - `sigma'))

* ---- macro scalars for the results table (BEFORE any panel surgery) --------
quietly summ year
local ylast = r(max)
foreach v in r_star f_x k_x y_x {
    capture confirm variable `v'
    if _rc gen `v' = .                       // RoW files may lack some
}
scalar rstar_ini = r_star[1]
quietly summ r_star if year == `ylast'
scalar rstar_fin = r(mean)
scalar nfa_ini   = f_x[1]
quietly summ f_x if year == `ylast'
scalar nfa_fin   = r(mean)
scalar ky_ini    = k_x[1] / y_x[1]
quietly summ k_x if year == `ylast'
local kk = r(mean)
quietly summ y_x if year == `ylast'
scalar ky_fin = `kk' / r(mean)

* ---- interior extreme of each macro series ("--" if monotone) --------------
capture program drop _ext
program define _ext, rclass
    args v
    tempvar d
    gen `d' = `v' - `v'[_n-1]
    quietly count if `d' > 1e-9  & !missing(`d')
    local nup = r(N)
    quietly count if `d' < -1e-9 & !missing(`d')
    local ndn = r(N)
    drop `d'
    if `nup' == 0 | `ndn' == 0 {
        return scalar ext = .
        return scalar t   = .
        exit
    }
    tempvar dev
    gen `dev' = abs(`v' - `v'[1])
    quietly summ `dev', meanonly
    quietly summ year if abs(`dev' - r(max)) < 1e-12
    local tx = r(min)
    drop `dev'
    quietly summ `v' if year == `tx', meanonly
    return scalar ext = r(mean)
    return scalar t   = `tx'
end

_ext r_star
scalar rstar_ext  = r(ext)
scalar rstar_extt = r(t)

_ext f_x
scalar nfa_ext  = r(ext)
scalar nfa_extt = r(t)

tempvar kyv
gen `kyv' = k_x / y_x
_ext `kyv'
scalar ky_ext  = r(ext)
scalar ky_extt = r(t)
drop `kyv'

* ============================================================================
* RA LIFE-CYCLES (Markov)
* ============================================================================
foreach c of global COHORTS {
    gen pi_y_`c' = .
    gen pi_m_`c' = .
    gen pi_o_`c' = .
    gen exp_c_uncond_`c'   = .
    gen exp_a_uncond_`c'   = .
    gen exp_c_cond_`c'     = .
    gen exp_a_cond_`c'     = .
    gen exp_ann_uncond_`c' = .
    gen exp_ann_cond_`c'   = .
    gen exp_a_own_cond_`c' = .
    gen pi_alive_`c'       = .

    replace pi_y_`c' = 1 if year == `c'
    replace pi_m_`c' = 0 if year == `c'
    replace pi_o_`c' = 0 if year == `c'
    replace exp_ann_uncond_`c' = 0 if year == `c'    // anchor

    forvalues age = 0/99 {
        local curr_yr = `c' + `age'
        local next_yr = `c' + `age' + 1

        replace pi_alive_`c' = pi_y_`c' + pi_m_`c' + pi_o_`c' if year == `curr_yr'

        replace exp_c_uncond_`c' = pi_y_`c'*c_y_indiv_x + pi_m_`c'*c_m_indiv_x ///
                                 + pi_o_`c'*c_o_indiv_x if year == `curr_yr'
        replace exp_a_uncond_`c' = pi_y_`c'*a_y_indiv_x + pi_m_`c'*a_m_indiv_x ///
                                 + pi_o_`c'*a_o_indiv_x if year == `curr_yr'
        replace exp_ann_uncond_`c' = (1 - gamma_x)*tilde_r_x*pi_o_`c'[_n-1]*a_o_indiv_x[_n-1] ///
            + tilde_r_x*gamma_x*exp_ann_uncond_`c'[_n-1] ///
            if year == `curr_yr' & year > `c'

        replace exp_c_cond_`c'     = exp_c_uncond_`c'   / pi_alive_`c' if year == `curr_yr'
        replace exp_a_cond_`c'     = exp_a_uncond_`c'   / pi_alive_`c' if year == `curr_yr'
        replace exp_ann_cond_`c'   = exp_ann_uncond_`c' / pi_alive_`c' if year == `curr_yr'
        replace exp_a_own_cond_`c' = exp_a_cond_`c' - exp_ann_cond_`c' if year == `curr_yr'

        summ gamma_x if year == `next_yr', meanonly
        local gamma_next = r(mean)
        replace pi_y_`c' = pi_y_`c'[_n-1]*`nu' if year == `next_yr'
        replace pi_m_`c' = pi_y_`c'[_n-1]*(1-`nu') + pi_m_`c'[_n-1]*`omega' if year == `next_yr'
        replace pi_o_`c' = pi_m_`c'[_n-1]*(1-`omega') + pi_o_`c'[_n-1]*`gamma_next' if year == `next_yr'
    }
}

* ---- expected lifetime consumption + per-cohort stats ----------------------
foreach c of global COHORTS {
    egen npv_`c' = total(exp_c_cond_`c') if inrange(year, `c', `c' + $HORIZON)

    quietly summ npv_`c'
    scalar s_npv_`c' = r(mean)
    quietly summ exp_c_cond_`c' if inrange(year, `c', `c' + $HORIZON)
    scalar s_pk_`c' = r(max)
    scalar s_mn_`c' = r(min)
    quietly summ year if inrange(year, `c', `c' + $HORIZON) ///
        & abs(exp_c_cond_`c' - s_pk_`c') < 1e-10
    scalar s_pka_`c' = r(min) - `c'
    quietly summ year if inrange(year, `c', `c' + $HORIZON) ///
        & abs(exp_c_cond_`c' - s_mn_`c') < 1e-10
    scalar s_mna_`c' = r(min) - `c'
}

* ---- bar chart -------------------------------------------------------------
graph bar npv_1973 npv_1990 npv_2006 npv_2022, ///
    title(" ", size(medium)) ///
    ytitle("Consumption (Efficiency Units)", size(medsmall)) ///
    legend(order(1 "1973 (Baby Boomer)" 2 "1990 (Gen X)" 3 "2006 (Gen Y)" 4 "2022 (Gen Z)") ///
           ring(1) position(6) cols(4) size(medsmall) symxsize(*0.5)) ///
    bar(1, color(navy)) bar(2, color(maroon)) bar(3, color(forest_green)) bar(4, color(dkorange)) ///
    blabel(bar, format(%9.2f)) ylabel(, format(%9.2f) labsize(medsmall)) ///
    graphregion(color(white)) name(npv_bar, replace)
graph export "`out'/lifetime_consumption_bar.pdf", as(pdf) name(npv_bar) replace
*"Expected Lifetime Consumption ($HORIZON Years)"

* ---- persist everything the next stage / table needs -----------------------
char _dta[tag]    "$TAG"
char _dta[region] "$REGION"
foreach p in beta sigma nu omega {
    char _dta[`p'] "``p''"
}
foreach s in rstar_ini rstar_fin nfa_ini nfa_fin ky_ini ky_fin ///
             rstar_ext rstar_extt nfa_ext nfa_extt ky_ext ky_extt {
    local tmp = string(`s', "%21.15g")
    char _dta[`s'] "`tmp'"
}
foreach c of global COHORTS {
    foreach st in npv pk pka mn mna {
        local tmp = string(s_`st'_`c', "%21.15g")
        char _dta[`st'_`c'] "`tmp'"
    }
}

save "Stata_Cache/${TAG}_pos.dta", replace
di as text "01_positive done -> Stata_Cache/${TAG}_pos.dta"


* ============================================================================
* RESULTS CSV — SA life-cycle, long by cohort x age -> Stata_Results/results_positive_$TAG.csv
* ============================================================================
cap mkdir "Stata_Results"
preserve
    local stubs "pi_y_ pi_m_ pi_o_ pi_alive_ exp_c_uncond_ exp_c_cond_ exp_a_uncond_ exp_a_cond_ exp_ann_cond_ exp_a_own_cond_ npv_"
    local keepv "year"
    foreach s of local stubs {
        foreach c of global COHORTS {
            local keepv "`keepv' `s'`c'"
        }
    }
    keep `keepv'
    reshape long `stubs', i(year) j(cohort)
    gen age = year - cohort
    keep if age >= 0 & age <= $HORIZON
    foreach v of varlist *_ {
        local nv = substr("`v'", 1, length("`v'") - 1)
        rename `v' `nv'
    }
    order cohort age year
    sort cohort age
    export delimited "Stata_Results/results_positive_${TAG}.csv", replace
    di as text "  -> Stata_Results/results_positive_${TAG}.csv"
restore
