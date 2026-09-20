* ============================================================================
* 02_NORMATIVE — Stata_Cache/$TAG_pos.dta -> C_life, Rawls (C*), KH (S), decomposition
*                -> Stata_Cache/$TAG_norm.dta  +  Stata_Cache/res_$TAG.dta (one-row table)
* Parameters come from the _dta chars written by 01_positive.
* ============================================================================
if "$TAG" == "" {
    di as error "globals not set — run via 00_master.do"
    exit 198
}
local out "Plots/$TAG"
cap mkdir "`out'"

use "Stata_Cache/${TAG}_pos.dta", clear
sort year

local beta  : char _dta[beta]
local sigma : char _dta[sigma]
local nu    : char _dta[nu]
local omega : char _dta[omega]

* ---- discounting, weights, entry wealth (BEFORE panel extension) -----------
gen ln_q_x     = ln(q_x)
gen sum_ln_q_x = sum(ln_q_x)
gen d          = exp(-sum_ln_q_x)            // D_j (discount to date 0)
gen n_entrants = phi_x * tilde_s_m_x         // N_j (mass of newborn young)
gen z_indiv    = c_y_indiv_x / mu_y_x        // Z_{j,0}^{RA}

* ---- extend the panel with the final BGP (flagged) -------------------------
local N0 = _N
local yr_last = year[`N0']
gen byte extended = 0
set obs `=_N + $HORIZON + 1'
replace extended = 1 if _n > `N0'
replace year = year[`N0'] + (_n - `N0') if extended
foreach v of varlist gamma_x c_y_indiv_x c_m_indiv_x c_o_indiv_x {
    replace `v' = `v'[`N0'] if extended
}

* ---- C_life for every entry year -------------------------------------------
gen C_life = .
forvalues j = $YEAR0/`yr_last' {
    local r0 = `j' - $YEAR0 + 1
    local py = 1
    local pm = 0
    local po = 0
    local tot = 0
    forvalues k = 0/$HORIZON {
        local r = `r0' + `k'
        local alive = `py' + `pm' + `po'
        local tot = `tot' + (`py'*c_y_indiv_x[`r'] + `pm'*c_m_indiv_x[`r'] ///
                             + `po'*c_o_indiv_x[`r'])/`alive'
        local g   = gamma_x[`=`r'+1']
        local pyn = `py'*`nu'
        local pmn = `py'*(1-`nu') + `pm'*`omega'
        local pon = `pm'*(1-`omega') + `po'*`g'
        local py = `pyn'
        local pm = `pmn'
        local po = `pon'
    }
    replace C_life = `tot' in `r0'
}

* ---- Rawlsian target and transfers -----------------------------------------
gen wgt   = d * n_entrants * z_indiv if !missing(C_life) & year <= $J_END
gen wgt_c = wgt / C_life
quietly summ wgt
scalar W_pool = r(sum)
quietly summ wgt_c
scalar sum_w  = r(sum)
scalar C_star = W_pool / sum_w
assert !missing(C_star)

gen v      = z_indiv * (C_star/C_life - 1) if year <= $J_END
gen v_rate = 100 * v / z_indiv

* self-financing check (IBC)
gen s_rawls = d * n_entrants * v
quietly summ s_rawls
di as text "Rawls IBC check (~0): " as result r(sum)
drop s_rawls

* break-even cohort
gen abs_v = abs(v) if inrange(year, $YEAR0, $J_END)
quietly summ abs_v, meanonly
gen aux = year if abs_v == r(min)
quietly summ aux, meanonly
scalar break_year = r(mean)
drop abs_v aux

* Rawls transfer plots (nominal + % of Z)
foreach spec in "v nominal Efficiency_Units" "v_rate percent %_of_Initial_Wealth" {
    tokenize `spec'
    local var  = "`1'"
    local suf  = "`2'"
    local ylab = subinstr("`3'","_"," ",.)

    gen v_tax_`suf' = min(`var', 0)
    gen v_sub_`suf' = max(`var', 0)
    gen zero_`suf'  = 0

    twoway ///
        (rarea zero_`suf' v_tax_`suf' year if inrange(year, $YEAR0, $J_END), ///
            color(navy%25) lcolor(navy%25)) ///
        (rarea zero_`suf' v_sub_`suf' year if inrange(year, $YEAR0, $J_END), ///
            color(maroon%25) lcolor(maroon%25)) ///
        (line `var' year if inrange(year, $YEAR0, $J_END), ///
            lcolor(black) lwidth(medthick)), ///
        yline(0, lcolor(gs8) lwidth(thin)) ///
        xline(`=break_year', lcolor(gs11) lpattern(dot)) ///
        text(0 `=break_year' "------> Break-even (`=break_year')", ///
             place(ne) size(medsmall) color(gs6)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year of Entry (Cohort)", size(medsmall)) ///
        ytitle("Rawlsian Transfer (`ylab')", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(1 "Taxed cohorts" 2 "Subsidized cohorts") ///
               ring(1) position(6) cols(2) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(rawls_transfer_`suf', replace)
    graph export "`out'/rawls_transfer_`suf'.pdf", as(pdf) ///
        name(rawls_transfer_`suf') replace 
		*width(2006)
		*"Optimal Rawlsian Transfer by Cohort"

    drop v_tax_`suf' v_sub_`suf' zero_`suf'
}

* ---- benchmarks -------------------------------------------------------------
* (A) Boomer cohort along the actual transition
quietly summ C_life if year == 1973
scalar C_boomer = r(mean)

* (B) pre-transition BGP (frozen 1970 consumptions and gamma)
scalar cy0 = c_y_indiv_x[1]
scalar cm0 = c_m_indiv_x[1]
scalar co0 = c_o_indiv_x[1]
scalar g0  = gamma_x[1]
scalar py = 1
scalar pm = 0
scalar po = 0
scalar C_BGP = 0
forvalues k = 0/$HORIZON {
    scalar alive = py + pm + po
    scalar C_BGP = C_BGP + (py*cy0 + pm*cm0 + po*co0)/alive
    scalar pyn = py*`nu'
    scalar pmn = py*(1-`nu') + pm*`omega'
    scalar pon = pm*(1-`omega') + po*g0
    scalar py = pyn
    scalar pm = pmn
    scalar po = pon
}

di as text _n "======= BENCHMARKS ======="
di "C_boomer (1973)        = " C_boomer
di "C_BGP    (pre-1970 SS) = " C_BGP
di "C_star   (Rawlsian)    = " C_star

* ---- KH test + welfare decomposition, per benchmark ------------------------
gen zz = d * n_entrants * z_indiv if !missing(C_life) & year <= $J_END
quietly summ zz
scalar Ztot = r(sum)

foreach bench in boomer BGP {
    scalar C_bench = C_`bench'
    local  suf     = "`bench'"

    gen v_kh_`suf' = z_indiv * (C_bench / C_life - 1) if !missing(C_life) & year <= $J_END
    gen s_kh_`suf' = d * n_entrants * v_kh_`suf'

    quietly summ s_kh_`suf'
    scalar Scost_`suf' = r(sum)
    scalar S_`suf'     = -Scost_`suf'
    scalar Spct_`suf'  = 100 * S_`suf' / Ztot

    di as text _n "----- KH test (SS0 = `suf') -----"
    di "S = -Delta_KH (>0 pass)   = " S_`suf'
    di "S as % of total resources = " Spct_`suf' "%"

    scalar eff_effect_`suf' = 100 * ((C_star / C_bench) - 1)
    gen total_effect_`suf'  = 100 * ((C_life / C_bench) - 1)
    gen equity_effect_`suf' = 100 * ((C_life / C_star) - 1)

    di "Aggregate efficiency effect (%) = " eff_effect_`suf'

    local benchlabel = cond("`suf'"=="boomer", "Boomer", "Pre-transition BGP")
    twoway ///
        (function y = `=eff_effect_`suf'', range($YEAR0 $J_END) ///
            lcolor(forest_green) lpattern(shortdash) lwidth(medthick)) ///
        (line total_effect_`suf'  year if inrange(year, $YEAR0, $J_END), ///
            lcolor(navy) lwidth(medthick)) ///
        (line equity_effect_`suf' year if inrange(year, $YEAR0, $J_END), ///
            lcolor(maroon) lpattern(dash) lwidth(medthick)), ///
        yline(0, lcolor(gs10) lwidth(thin)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year of Entry (Cohort)", size(medsmall)) ///
        ytitle("% Variation in Lifetime Consumption", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(2 "Total Effect (vs SS0)" ///
                     3 "Equity Effect (vs Rawlsian)" ///
                     1 "Efficiency Effect") ///
               ring(1) position(6) cols(3) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) ///
        name(welfare_decomp_`suf', replace)
    graph export "`out'/welfare_decomp_`suf'.pdf", as(pdf) ///
        name(welfare_decomp_`suf') replace 
		*width(2006)
		
}
*"Intergenerational Welfare Decomposition (SS0 = `benchlabel')"

di as text _n "  SUMMARY               Boomer         BGP"
di "  S = -Delta_KH   :   " %10.4f S_boomer "  " %10.4f S_BGP
di "  Efficiency eff. :   " %10.4f eff_effect_boomer "  " %10.4f eff_effect_BGP

save "Stata_Cache/${TAG}_norm.dta", replace

* ============================================================================
* ONE-ROW RESULTS DATASET (feeds 04_table)
* ============================================================================
local tag    : char _dta[tag]
local region : char _dta[region]
foreach s in rstar_ini rstar_fin nfa_ini nfa_fin ky_ini ky_fin ///
             rstar_ext rstar_extt nfa_ext nfa_extt ky_ext ky_extt {
    local `s' : char _dta[`s']
}
foreach c of global COHORTS {
    foreach st in npv pk pka mn mna {
        local `st'`c' : char _dta[`st'_`c']
    }
}

clear
set obs 1
gen str20 tag    = "`tag'"
gen str8  region = "`region'"
gen double beta  = `beta'
gen double sigma = `sigma'
foreach s in rstar_ini rstar_fin ky_ini ky_fin nfa_ini nfa_fin ///
             rstar_ext rstar_extt ky_ext ky_extt nfa_ext nfa_extt {
    gen double `s' = ``s''
}
foreach c of global COHORTS {
    gen double npv_`c' = `npv`c''
    gen double pk_`c'  = `pk`c''
    gen double pka_`c' = `pka`c''
    gen double mn_`c'  = `mn`c''
    gen double mna_`c' = `mna`c''
}
gen double C_star     = scalar(C_star)
gen double C_boomer   = scalar(C_boomer)
gen double C_BGP      = scalar(C_BGP)
gen double break_year = scalar(break_year)
gen double S_pct_boomer = scalar(Spct_boomer)
gen double S_pct_BGP    = scalar(Spct_BGP)
gen double eff_boomer   = scalar(eff_effect_boomer)
gen double eff_BGP      = scalar(eff_effect_BGP)

save "Stata_Cache/res_${TAG}.dta", replace
di as text "02_normative done -> Stata_Cache/${TAG}_norm.dta + Stata_Cache/res_${TAG}.dta"
