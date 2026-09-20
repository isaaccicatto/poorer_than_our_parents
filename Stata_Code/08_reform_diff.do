* ============================================================================
* 08_REFORM_DIFF — counterfactual minus baseline, RoW only.
*   Pairs (edit REFORMS / REGION below if needed):
*       base_RoW  vs  repl_RoW      (pension de-indexation)
*       base_RoW  vs  by_RoW        (public debt)
*   Inputs : Stata_Cache/base_RoW_norm.dta, Stata_Cache/<reform>_RoW_norm.dta
*   Outputs: Plots/counterfactual_diff/
*     <reform>_diff_panels.pdf          consumption | net assets | earnings
*                                       (index rel. to YEAR0, cf - base, p.p.;
*                                        x-axis starts at the reform year)
*     <reform>_diff_cons_ratio_age.pdf  left : ratio vs Boomer by age, cf - base
*                                             (dotted line = age at which the
*                                              Boomer, the denominator, is hit)
*                                       right: lifetime consumption, cf vs base
*                                              (% change) per cohort, bars
*   Reform year = first year the two runs diverge (unanticipated shock).
*   Run AFTER the scenario loop in 00_master. Missing cache -> pair skipped.
*   Never saves data.
* ============================================================================
if "$YEAR0" == "" | "$J_END" == "" | "$HORIZON" == "" | "$COHORTS" == "" {
    di as error "globals not set (YEAR0 J_END HORIZON COHORTS) — run via 00_master.do"
    exit 198
}
local REGION  "RoW"
local REFORMS "repl by"
local out "Plots/counterfactual_diff"
cap mkdir "Plots"
cap mkdir "`out'"
local win "inrange(year, $YEAR0, $J_END)"

* ---- loader: keep what the figures need, tag with a suffix ----------------
capture program drop _prep
program define _prep
    args tag sfx
    use "Stata_Cache/`tag'_norm.dta", clear
    keep if !extended
    sort year
    foreach c in y m o {
        capture confirm variable tilde_s_`c'_x
        if _rc {
            capture confirm variable share_`c'_x
            if !_rc rename share_`c'_x tilde_s_`c'_x
        }
    }
    capture confirm variable e_indiv_x
    if _rc {
        capture confirm variable e_x
        if !_rc gen e_indiv_x = e_x / tilde_s_o_x
        else    gen e_indiv_x = .
    }
    capture confirm variable w_net_indiv_x
    if _rc gen w_net_indiv_x = w_indiv_x - t_indiv_x
    keep year c_y_indiv_x c_m_indiv_x c_o_indiv_x ///
              a_y_indiv_x a_m_indiv_x a_o_indiv_x ///
              w_indiv_x t_indiv_x w_net_indiv_x e_indiv_x exp_c_cond_*
    rename (c_y_indiv_x c_m_indiv_x c_o_indiv_x a_y_indiv_x a_m_indiv_x a_o_indiv_x ///
            w_indiv_x t_indiv_x w_net_indiv_x e_indiv_x) ///
           (c_y_`sfx' c_m_`sfx' c_o_`sfx' a_y_`sfx' a_m_`sfx' a_o_`sfx' ///
            w_`sfx' t_`sfx' wn_`sfx' e_`sfx')
    rename exp_c_cond_* exp_c_`sfx'_*
end

local base "base_`REGION'"
capture confirm file "Stata_Cache/`base'_norm.dta"
if _rc {
    di as error "*** 07_reform_diff: Stata_Cache/`base'_norm.dta AUSENTE — nada a fazer ***"
    exit
}

* cohort labels (bar chart + legends)
capture program drop _cohlab
program define _cohlab, rclass
    args c
    local lab "`c'"
    if "`c'" == "1973" local lab "Boomer (1973)"
    if "`c'" == "1990" local lab "Gen X (1990)"
    if "`c'" == "2006" local lab "Gen Y (2006)"
    if "`c'" == "2022" local lab "Gen Z (2022)"
    return local lab "`lab'"
end

foreach rf of local REFORMS {
    local cf "`rf'_`REGION'"
    capture confirm file "Stata_Cache/`cf'_norm.dta"
    if _rc {
        di as error "*** 07_reform_diff: Stata_Cache/`cf'_norm.dta AUSENTE — par `cf' PULADO ***"
        continue
    }

    di as text _n "{hline 60}"
    di as text "  07_reform_diff   `cf'  minus  `base'   -> `out'/`rf'_diff_*.pdf"
    di as text "{hline 60}"

    _prep `base' b
    tempfile B
    save `B'
    _prep `cf' c
    merge 1:1 year using `B'
    capture assert _merge == 3
    if _rc di as error "07_reform_diff[`cf']: AVISO — anos nao coincidem; uso a intersecao"
    keep if _merge == 3
    drop _merge
    sort year
    capture assert year[1] == $YEAR0
    if _rc di as error "07_reform_diff[`cf']: AVISO — primeiro ano != $YEAR0"

    * ---- reform year: first year the two runs diverge ---------------------
    tempvar dev
    gen `dev' = max(abs(c_y_c - c_y_b), abs(c_m_c - c_m_b), abs(c_o_c - c_o_b), ///
                    abs(w_c - w_b), abs(t_c - t_b), abs(e_c - e_b))
    quietly summ year if `dev' > 1e-9, meanonly
    if r(N) == 0 {
        di as error "*** 07_reform_diff[`cf']: base e cf identicos — par PULADO ***"
        continue
    }
    local t_ref = r(min)
    local xstep = max(5, round(($J_END - `t_ref') / 3))
    di as text "07_reform_diff[`cf']: reform detected at year `t_ref'"

    * ---- differences of the plotted ratios (p.p. of the YEAR0 level) ------
    foreach v in c_y c_m c_o a_y a_m a_o w t wn e {
        gen d_`v' = 100 * (`v'_c / `v'_c[1] - `v'_b / `v'_b[1])
    }
    di as text "  max |diff| over the window (p.p. of $YEAR0 level):"
    foreach v in c_y c_m c_o a_y a_m a_o w t wn e {
        quietly summ d_`v' if `win'
        if r(N) > 0 di as text "    d_`v'" _col(12) as result %8.3f max(abs(r(min)), abs(r(max)))
    }

    * ====================================================================
    * PANELS — consumption | net assets | earnings, from the reform year on
    * ====================================================================
    local post "year >= `t_ref' & year <= $J_END"

    twoway (line d_c_y year if `post', lcolor(navy)         lwidth(medthick)) ///
           (line d_c_m year if `post', lcolor(dkorange)     lwidth(medthick)) ///
           (line d_c_o year if `post', lcolor(forest_green) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(black)) ///
        title("Consumption", size(medsmall)) xtitle("") ///
        xlabel(`t_ref'(`xstep')$J_END, labsize(small)) ///
        ylabel(, format(%9.2f) labsize(small) angle(0)) ///
        ytitle("Counterfactual - Baseline (p.p. of $YEAR0 level)", size(small)) ///
        legend(order(1 "Young" 2 "Middle-Aged" 3 "Retirees") ///
               ring(1) position(6) rows(1) size(small) symxsize(*0.5) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(p_cons, replace)

    twoway (line d_a_y year if `post', lcolor(navy)         lwidth(medthick)) ///
           (line d_a_m year if `post', lcolor(dkorange)     lwidth(medthick)) ///
           (line d_a_o year if `post', lcolor(forest_green) lwidth(medthick)), ///
        yline(0, lpattern(dash) lcolor(black)) ///
        title("Net assets", size(medsmall)) xtitle("") ///
        xlabel(`t_ref'(`xstep')$J_END, labsize(small)) ///
        ylabel(, format(%9.2f) labsize(small) angle(0)) ///
        ytitle("", size(small)) ///
        legend(off) ///
        graphregion(color(white)) plotregion(color(white)) name(p_asst, replace)

    * earnings: palette deliberately distinct from the cohort colours
    quietly count if !missing(d_e) & `post'
    if r(N) > 0 {
        local p_e `"(line d_e year if `post', lcolor(purple) lwidth(medthick) lpattern(dash))"'
        local leg `"1 "Gross Wage" 2 "Taxes" 3 "Net Wage" 4 "Pensions""'
    }
    else {
        local p_e ""
        local leg `"1 "Gross Wage" 2 "Taxes" 3 "Net Wage""'
    }
    twoway (line d_w  year if `post', lcolor(gs4)       lwidth(medthick)) ///
           (line d_t  year if `post', lcolor(cranberry) lwidth(medthick)) ///
           (line d_wn year if `post', lcolor(teal)      lwidth(medthick)) ///
           `p_e', ///
        yline(0, lpattern(dash) lcolor(black)) ///
        title("Earnings", size(medsmall)) xtitle("") ///
        xlabel(`t_ref'(`xstep')$J_END, labsize(small)) ///
        ylabel(, format(%9.2f) labsize(small) angle(0)) ///
        ytitle("", size(small)) ///
        legend(order(`leg') ///
               ring(1) position(6) rows(2) size(small) symxsize(*0.5) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(p_earn, replace)

    graph combine p_cons p_asst, cols(2) iscale(0.85) ///
        xsize(10) ysize(5) graphregion(color(white)) name(diff_panels, replace)
    graph export "`out'/`rf'_diff_panels.pdf", as(pdf) name(diff_panels) replace

    * Exporta Earnings (p_earn) separadamente
    graph export "`out'/`rf'_diff_earnings.pdf", as(pdf) name(p_earn) replace

    * ====================================================================
    * RATIO vs Boomer by age (cf - base)  |  lifetime consumption, % change
    * ====================================================================
    local ref : word 1 of $COHORTS
    local others ""
    foreach c of global COHORTS {
        if "`c'" != "`ref'" local others "`others' `c'"
    }
    local a_ref = `t_ref' - `ref'          // age at which the denominator is hit

    preserve
        keep year exp_c_b_* exp_c_c_*
        reshape long exp_c_b_ exp_c_c_, i(year) j(cohort)
        gen age = year - cohort
        keep if age >= 0 & age <= $HORIZON
        drop year
        reshape wide exp_c_b_ exp_c_c_, i(age) j(cohort)

        * lifetime consumption (sum over ages 0..HORIZON), % change cf vs base
        local k = 1
        foreach c of global COHORTS {
            quietly summ exp_c_b_`c'
            local sb = r(sum)
            quietly summ exp_c_c_`c'
            local sc = r(sum)
            local pct`k' = 100 * (`sc' / `sb' - 1)
            _cohlab `c'
            local lab`k' "`r(lab)'"
            di as text "  lifetime consumption `c': " as result %7.3f `pct`k'' as text " %"
            local ++k
        }
        local ncoh = `k' - 1

        * left panel: ratio differences by age
        local plots ""
        local leg   ""
        local k = 1
        local cols "maroon forest_green dkorange navy"
        foreach c of local others {
            local col : word `k' of `cols'
            gen d_r`c' = 100 * (exp_c_c_`c' / exp_c_c_`ref' - exp_c_b_`c' / exp_c_b_`ref')
            local plots `"`plots' (line d_r`c' age, lcolor(`col') lwidth(medthick))"'
            _cohlab `c'
            local leg `"`leg' `k' "`r(lab)'""'
            local ++k
        }
        local nleg = `k' - 1
        quietly summ d_r`: word 1 of `others''
        local ymax = r(max)
        local ymin = r(min)
        foreach c of local others {
            quietly summ d_r`c'
            if r(max) > `ymax' local ymax = r(max)
            if r(min) < `ymin' local ymin = r(min)
        }
        local ytxt = `ymin' + 0.06 * (`ymax' - `ymin')

        twoway `plots', ///
            yline(0, lpattern(dash) lcolor(black)) ///
            xline(`a_ref', lcolor(gs8) lpattern(dot)) ///
            text(`ytxt' `a_ref' " ---> Impact on Boomer (denom.)", place(e) size(small) color(gs6)) ///
            title("Consumption ratio vs Boomer", size(medsmall)) ///
            xtitle("Periods Since Entry", size(small)) ///
            ytitle("Counterfactual - Baseline (p.p. of ratio to `ref')", size(small)) ///
            ylabel(, format(%9.2f) labsize(small)) xlabel(, labsize(small)) ///
            legend(order(`leg') ///
                   ring(1) position(6) rows(1) size(small) symxsize(*0.5) region(lcolor(white))) ///
            graphregion(color(white)) plotregion(color(white)) name(p_ratio, replace)

        * right panel: bars, % change in lifetime consumption per cohort
        clear
        set obs `ncoh'
        gen byte idx = _n
        gen str24 lab = ""
        gen double pct = .
        forvalues k = 1/`ncoh' {
            replace lab = "`lab`k''" in `k'
            replace pct = `pct`k''   in `k'
        }
        local bcols "navy maroon forest_green dkorange"
        local baropts ""
        forvalues k = 1/`ncoh' {
            local bc : word `k' of `bcols'
            local baropts "`baropts' bar(`k', color(`bc'))"
        }
        graph bar (asis) pct, over(lab, sort(idx) label(labsize(small))) ///
            asyvars showyvars legend(off) `baropts' ///
            yline(0, lcolor(gs8) lwidth(thin)) ///
            blabel(bar, format(%9.2f) size(small)) ///
            title("Lifetime consumption", size(medsmall)) ///
            ytitle("Counterfactual vs Baseline (%)", size(small)) ///
            ylabel(, format(%9.2f) labsize(small) angle(0)) ///
            graphregion(color(white)) plotregion(color(white)) name(p_bar, replace)
    restore

    graph combine p_ratio p_bar, cols(2) iscale(0.9) ///
        xsize(12) ysize(5) graphregion(color(white)) name(diff_ratio, replace)
    graph export "`out'/`rf'_diff_cons_ratio_age.pdf", as(pdf) name(diff_ratio) replace

    di as text "07_reform_diff[`cf'] done -> `out'/`rf'_diff_panels.pdf, `rf'_diff_cons_ratio_age.pdf"
}
