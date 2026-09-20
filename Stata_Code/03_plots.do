* ============================================================================
* 03_PLOTS — Stata_Cache/$TAG_norm.dta -> earnings, indiv index, cross-age and
*            ratio plots. The destructive reshapes live HERE and only here.
*            Never saves data.
* Chronological plots: x = year, window $YEAR0..$J_END (t <= 100), no xtitle.
* Pensions are omitted from the earnings plots: E_indiv = repl * W_indiv, so
* they sit exactly on top of the gross-wage line.
* ============================================================================
if "$TAG" == "" {
    di as error "globals not set — run via 00_master.do"
    exit 198
}
local out "Plots/$TAG"
cap mkdir "`out'"
local win "inrange(year, $YEAR0, $J_END)"

use "Stata_Cache/${TAG}_norm.dta", clear
sort year

* ---- earnings block (skipped automatically if the file lacks wage vars) ----
* Pensions: if e_indiv / w_indiv is constant over the window (E = repl x W),
* they grow exactly like gross wages and are omitted (they would sit on top of
* the wage line). Otherwise (decoupled pensions, e.g. repl/by scenarios) a
* fourth line is added to both plots.
capture confirm variable w_indiv_x
if !_rc {
    capture confirm variable w_net_indiv_x
    if _rc gen w_net_indiv_x = w_indiv_x - t_indiv_x

    * pension per retiree: use e_indiv_x if present, else e_x / share_o
    capture confirm variable e_indiv_x
    if _rc {
        capture confirm variable e_x
        if !_rc {
            capture confirm variable tilde_s_o_x
            if !_rc gen e_indiv_x = e_x / tilde_s_o_x
            else {
                capture confirm variable share_o_x
                if !_rc gen e_indiv_x = e_x / share_o_x
            }
        }
    }
    capture confirm variable e_indiv_x
    local havepens = (_rc == 0)

    * do pensions track gross wages? (ratio constant up to rounding)
    local samegrowth = 1
    if `havepens' {
        tempvar rat
        gen `rat' = e_indiv_x / w_indiv_x if !extended & `win'
        quietly summ `rat'
        if r(N) > 0 & (r(max) - r(min)) > 1e-6 * abs(r(mean)) local samegrowth = 0
        drop `rat'
    }
    di as text "03_plots[$TAG]: pensions " cond(`samegrowth', "track gross wages (E = repl x W) -> omitted", "decoupled from wages -> plotted")

    if `samegrowth' {
        local p_e_g   ""
        local p_e_l   ""
        local leg     `"1 "Gross Wage (= Pensions)" 2 "Taxes" 3 "Net Wage""'
        local ncol    3
        local noteopt `"note("Pensions grow at the same rate as gross wages, by construction (E = repl x W).", size(medsmall))"'
    }
    else {
        local p_e_g   `"(line e_growth year if `win', lcolor(purple) lwidth(medthick) lpattern(dash))"'
        local p_e_l   `"(line e_lvl    year if `win', llcolor(purple) lwidth(medthick) lpattern(dash))"'
        local leg     `"1 "Gross Wage" 2 "Taxes" 3 "Net Wage" 4 "Pensions""'
        local ncol    4
        local noteopt ""
    }

    * ---- yearly growth ------------------------------------------------------
    gen w_growth     = 100 * (w_indiv_x     / w_indiv_x[_n-1]     - 1) if !extended
    gen t_growth     = 100 * (t_indiv_x     / t_indiv_x[_n-1]     - 1) if !extended
    gen w_net_growth = 100 * (w_net_indiv_x / w_net_indiv_x[_n-1] - 1) if !extended
    if !`samegrowth' gen e_growth = 100 * (e_indiv_x / e_indiv_x[_n-1] - 1) if !extended

    twoway (line w_growth     year if `win', lcolor(gs4)    lwidth(medthick)) ///
           (line t_growth     year if `win', lcolor(cranberry)   lwidth(medthick)) ///
           (line w_net_growth year if `win', lcolor(teal) lwidth(medthick)) ///
           `p_e_g', ///
        title(" ", size(medium)) ///
        xtitle("") ///
        ytitle("Yearly Growth (%)", size(medsmall)) ///
        legend(order(`leg') ///
               ring(1) position(6) cols(`ncol') size(medsmall) symxsize(*0.5)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall)) ///
        `noteopt' ///
        graphregion(color(white)) name(earn_time, replace)
    graph export "`out'/wage_chronological_cond.pdf", as(pdf) name(earn_time) replace
    *"Earnings Growth Over Time"

    * ---- levels relative to $YEAR0 -------------------------------------------
    local lvlvars "w_indiv_x t_indiv_x w_net_indiv_x"
    if !`samegrowth' local lvlvars "`lvlvars' e_indiv_x"
    foreach v of local lvlvars {
        gen `v'_0 = `v' if year == $YEAR0
        replace `v'_0 = `v'_0[_n-1] if year > $YEAR0
    }
    gen w_lvl     = w_indiv_x     / w_indiv_x_0
    gen t_lvl     = t_indiv_x     / t_indiv_x_0
    gen w_net_lvl = w_net_indiv_x / w_net_indiv_x_0
    if !`samegrowth' gen e_lvl = e_indiv_x / e_indiv_x_0

    twoway (line w_lvl     year if `win', lcolor(gs4)     lwidth(medthick)) ///
           (line t_lvl     year if `win', lcolor(cranberry)  lwidth(medthick)) ///
           (line w_net_lvl year if `win', lcolor(teal) lwidth(medthick)) ///
           `p_e_l', ///
        yline(1, lpattern(dash) lcolor(black)) ///
        title(" ", size(medium)) ///
        xtitle("") ///
        ytitle("Ratio (Current Year / $YEAR0)", size(medsmall)) ///
        legend(order(`leg') ///
               ring(1) position(6) cols(`ncol') size(medsmall) symxsize(*0.5)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall)) ///
        graphregion(color(white)) name(earn_2_time, replace)
    graph export "`out'/wage_2_chronological_cond.pdf", as(pdf) name(earn_2_time) replace
    *"Earnings Over Time"
}
* ============================================================================
* INDIVIDUAL CONSUMPTION | NET ASSETS by cohort, index ($YEAR0 = 1), by(panel)
*   left  = c_y/c_m/c_o (indiv) relative to $YEAR0
*   right = a_y/a_m/a_o (indiv) relative to $YEAR0
* ============================================================================
local haveind = 1
foreach v in c_y_indiv_x c_m_indiv_x c_o_indiv_x a_y_indiv_x a_m_indiv_x a_o_indiv_x {
    capture confirm variable `v'
    if _rc local haveind = 0
}
if !`haveind' di as error "03_plots[$TAG]: sem c_/a_*_indiv_x -> indice individual PULADO"

if `haveind' {
preserve
    keep year c_y_indiv_x c_m_indiv_x c_o_indiv_x a_y_indiv_x a_m_indiv_x a_o_indiv_x
    keep if `win'
    sort year
    foreach v in c_y c_m c_o a_y a_m a_o {
        gen `v'_lvl = `v'_indiv_x / `v'_indiv_x[1]
    }

    * stack into two panels sharing the same three series names
    expand 2, gen(dup)
    gen byte panel = 1 + dup
    drop dup
    label define pnl2 1 "Consumption" 2 "Net assets", replace
    label values panel pnl2
    gen v_y = cond(panel == 1, c_y_lvl, a_y_lvl)
    gen v_m = cond(panel == 1, c_m_lvl, a_m_lvl)
    gen v_o = cond(panel == 1, c_o_lvl, a_o_lvl)

    twoway ///
        (line v_y year, lcolor(navy)         lwidth(medthick)) ///
        (line v_m year, lcolor(dkorange)     lwidth(medthick)) ///
        (line v_o year, lcolor(forest_green) lwidth(medthick)), ///
        by(panel, cols(2) note("") legend(position(6) ring(1)) iscale(1) graphregion(color(white)) yrescale) ///
        xsize(11) ysize(5) ///
        yline(1, lpattern(dash) lcolor(black)) ///
        subtitle(, size(medsmall) nobox) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        xtitle("") ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        ytitle("Ratio (Current Year / $YEAR0)", size(medsmall)) ///
        legend(order(1 "Young" 2 "Middle-Aged" 3 "Retirees") ///
               cols(3) size(medsmall) symxsize(*0.5) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(indiv_idx, replace)
    graph export "`out'/indiv_cons_assets_index.pdf", as(pdf) name(indiv_idx) replace
    *"Individual Consumption and Net Assets by Cohort (index, $YEAR0 = 1)"
restore
}


* ============================================================================
* CROSS-AGE PLOTS (destructive from here on — data is never saved)
* ============================================================================
keep year exp_c_cond_* exp_a_cond_* exp_a_own_cond_* pi_alive_*
reshape long exp_c_cond_ exp_a_cond_ exp_a_own_cond_ pi_alive_, i(year) j(cohort)
gen age = year - cohort
keep if age >= 0 & age <= $HORIZON

replace pi_alive_ = 100 * pi_alive_

twoway (line exp_c_cond_ age if cohort == 1973, lcolor(navy)         lwidth(medthick)) ///
       (line exp_c_cond_ age if cohort == 1990, lcolor(maroon)       lwidth(medthick)) ///
       (line exp_c_cond_ age if cohort == 2006, lcolor(forest_green) lwidth(medthick)) ///
       (line exp_c_cond_ age if cohort == 2022, lcolor(dkorange)     lwidth(medthick)), ///
    title(" ", size(medium)) ///
    xtitle("Periods Since Entry", size(medsmall)) ///
    ytitle("Consumption (Efficiency Units)", size(medsmall)) ///
    legend(order(1 "1973 (Baby Boomer)" 2 "1990 (Gen X)" 3 "2006 (Gen Y)" 4 "2022 (Gen Z)") ///
           ring(1) position(6) cols(4) size(medsmall) symxsize(*0.5)) ///
    xlabel(, labsize(medsmall)) ylabel(, format(%9.2f) labsize(medsmall)) ///
    graphregion(color(white)) name(cons_age, replace)
graph export "`out'/cons_cross_age_cond.pdf", as(pdf) name(cons_age) replace
*"Expected Life-Cycle Consumption (Conditional on Survival)"

twoway (line exp_a_cond_ age if cohort == 1973, lcolor(navy)         lwidth(medthick)) ///
       (line exp_a_cond_ age if cohort == 1990, lcolor(maroon)       lwidth(medthick)) ///
       (line exp_a_cond_ age if cohort == 2006, lcolor(forest_green) lwidth(medthick)) ///
       (line exp_a_cond_ age if cohort == 2022, lcolor(dkorange)     lwidth(medthick)) ///
       (line exp_a_own_cond_ age if cohort == 1973, lcolor(navy)         lwidth(medthick) lpattern(dash)) ///
       (line exp_a_own_cond_ age if cohort == 1990, lcolor(maroon)       lwidth(medthick) lpattern(dash)) ///
       (line exp_a_own_cond_ age if cohort == 2006, lcolor(forest_green) lwidth(medthick) lpattern(dash)) ///
       (line exp_a_own_cond_ age if cohort == 2022, lcolor(dkorange)     lwidth(medthick) lpattern(dash)), ///
    title(" ", size(medium)) ///
    xtitle("Periods Since Entry", size(medsmall)) ///
    ytitle("Net Assets (Efficiency Units)", size(medsmall)) ///
    legend(order(1 "1973 (Baby Boomer)" 2 "1990 (Gen X)" 3 "2006 (Gen Y)" ///
                 4 "2022 (Gen Z)" 5 "Ex-Annuity (Own)" 6 "" 7 "" 8 "") ///
           ring(1) position(6) rows(2) size(medsmall) symxsize(*0.5)) ///
    xlabel(, labsize(medsmall)) ylabel(, format(%9.2f) labsize(medsmall)) ///
    graphregion(color(white)) name(asst_age, replace)
graph export "`out'/asst_cross_age_cond.pdf", as(pdf) name(asst_age) replace
*"Expected Life-Cycle Net Asset Position (Conditional on Survival)"

twoway (line pi_alive_ age if cohort == 1973, lcolor(navy)         lwidth(medthick)) ///
       (line pi_alive_ age if cohort == 1990, lcolor(maroon)       lwidth(medthick)) ///
       (line pi_alive_ age if cohort == 2006, lcolor(forest_green) lwidth(medthick)) ///
       (line pi_alive_ age if cohort == 2022, lcolor(dkorange)     lwidth(medthick)), ///
    title(" ", size(medium)) ///
    xtitle("Periods Since Entry", size(medsmall)) ///
    ytitle("Probability of Being Alive (%)", size(medsmall)) ///
    ylabel(0(20)100, format(%9.2f) labsize(medsmall)) xlabel(, labsize(medsmall)) ///
    legend(order(1 "1973 (Baby Boomer)" 2 "1990 (Gen X)" 3 "2006 (Gen Y)" 4 "2022 (Gen Z)") ///
           ring(1) position(6) cols(4) size(medsmall) symxsize(*0.5)) ///
    graphregion(color(white)) name(survival, replace)
graph export "`out'/survival.pdf", as(pdf) name(survival) replace
*"Survival Probability by Generation"

* ---- ratios (descendant / Baby Boomer) -------------------------------------
drop year exp_a_own_cond_ exp_a_cond_ pi_alive_
reshape wide exp_c_cond_, i(age) j(cohort)

gen ratio_c_1990 = 100 * (exp_c_cond_1990 / exp_c_cond_1973)
gen ratio_c_2006 = 100 * (exp_c_cond_2006 / exp_c_cond_1973)
gen ratio_c_2022 = 100 * (exp_c_cond_2022 / exp_c_cond_1973)

twoway (line ratio_c_1990 age, lcolor(maroon)       lwidth(medthick)) ///
       (line ratio_c_2006 age, lcolor(forest_green) lwidth(medthick)) ///
       (line ratio_c_2022 age, lcolor(dkorange)     lwidth(medthick)), ///
    yline(100, lpattern(dash) lcolor(black)) ///
    title(" ", size(medium)) ///
    xtitle("Periods Since Entry", size(medsmall)) ///
    ytitle("Ratio (Descendant / Baby Boomer, %)", size(medsmall)) ///
    ylabel(, format(%9.2f) labsize(medsmall)) xlabel(, labsize(medsmall)) ///
    legend(order(1 "Gen X (1990)" 2 "Gen Y (2006)" 3 "Gen Z (2022)") ///
           ring(1) position(6) cols(4) size(medsmall) symxsize(*0.5)) ///
    graphregion(color(white)) name(cons_ratio_age, replace)
graph export "`out'/cons_ratio_age_cond.pdf", as(pdf) name(cons_ratio_age) replace
*"Consumption Ratio vs 1973 Cohort (Conditional on Survival)"

di as text "03_plots done -> `out'/"


