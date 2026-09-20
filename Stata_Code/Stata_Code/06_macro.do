* ============================================================================
* 06_MACRO — Stata_Cache/${PFX}_P_pos.dta + Stata_Cache/${PFX}_A_pos.dta
*            -> cross-country macro plots in Plots/macro_analysis/
*   M1) interest-rate differential vs RoW:  R_tilde - R*        (P, A)
*   M2) NFA trajectories                                         (P, A)
*   M3) INDIVIDUAL net wage & pension                            (P, A) x 2
*   M4) AGGREGATE pensions/GDP & tax burden/GDP                  (P, A) x 2
* Runs AFTER the scenario loop in 00_master (needs both caches). No $TAG.
* $PFX picks the scenario family ("base", "beta_high", ...); default "base".
* Country labels follow the TAG suffix (paper labels), exactly as 04_table:
*   P = ${PFX}_P  (<- CSV suffix _A)     A = ${PFX}_A  (<- CSV suffix _B)
* Colour rule: one hue family per country (P orange, A purple);
*   DARK  (solid)  = net wage / tax burden  (W-T, T)
*   LIGHT (dashed) = pensions               (E)
* All graphs in Times New Roman (graph set window fontface).
* ============================================================================
if "$YEAR0" == "" | "$J_END" == "" {
    di as error "globals not set — run via 00_master.do"
    exit 198
}
if "$PFX" == "" global PFX "base"
graph set window fontface "Times New Roman"
local out "Plots/macro_analysis"
cap mkdir "Plots"
cap mkdir "`out'"
local win "inrange(year, $YEAR0, $J_END)"
local lw  "medthick"

* ---- palette (RGB) ---------------------------------------------------------
* preset 1 (default): P orange, A purple
local P_d "203 106 43"
local P_l "245 190 150"
local A_d "94 60 153"
local A_l "180 160 215"
* preset 2: P green, A red  -> comment preset 1 and uncomment these
*local P_d "63 124 72"
*local P_l "150 200 160"
*local A_d "144 53 59"
*local A_l "220 140 145"

* ---- helpers ---------------------------------------------------------------
capture program drop _hasall
program define _hasall, rclass
    local ok = 1
    foreach v of local 0 {
        capture confirm variable `v'
        if _rc local ok = 0
    }
    return scalar ok = `ok'
end

capture program drop _chkfile
program define _chkfile
    args f
    capture confirm file "`f'"
    if _rc di as error "*** export FALHOU: `f' ***"
    else   di as text  "  -> `f' OK"
end

* ---- load both SOEs, keep what we need, suffix every column by country -----
local vars "tilde_r_x r_star f_x y_x e_x t_x w_indiv_x t_indiv_x e_indiv_x w_net_indiv_x"

foreach e in P A {
    capture confirm file "Stata_Cache/${PFX}_`e'_pos.dta"
    if _rc {
        di as error "*** 06_macro: Stata_Cache/${PFX}_`e'_pos.dta AUSENTE — rode o cenario ${PFX}_`e' antes ***"
        exit 601
    }
    use "Stata_Cache/${PFX}_`e'_pos.dta", clear
    sort year

    * net wage: use the CSV column if present, otherwise build it (03_plots style)
    capture confirm variable w_net_indiv_x
    if _rc {
        capture confirm variable w_indiv_x
        local rc1 = _rc
        capture confirm variable t_indiv_x
        if `rc1' == 0 & _rc == 0 gen w_net_indiv_x = w_indiv_x - t_indiv_x
    }

    local keep "year"
    foreach v of local vars {
        capture confirm variable `v'
        if !_rc local keep "`keep' `v'"
        else    di as error "06_macro[`e']: coluna `v' ausente — grafico que a usa sera pulado"
    }
    keep `keep'
    rename *_x *_`e'
    capture rename r_star r_star_`e'
    tempfile d_`e'
    save `d_`e''
}

use `d_P', clear
merge 1:1 year using `d_A'
capture assert _merge == 3
if _rc di as error "06_macro: AVISO — anos nao coincidem entre P e A; uso a intersecao"
keep if _merge == 3
drop _merge
sort year
capture assert year[1] == $YEAR0
if _rc di as error "06_macro: AVISO — primeiro ano != $YEAR0"

* R* is a world price: the two caches must agree
_hasall r_star_P r_star_A
if r(ok) {
    capture assert abs(r_star_P - r_star_A) < 1e-8
    if _rc di as error "06_macro: AVISO — r_star difere entre P e A (caches de mundos diferentes?)"
}

* ============================================================================
* M1 — interest-rate differential vs RoW (percentage points)
* ============================================================================
_hasall tilde_r_P tilde_r_A r_star_P r_star_A
if r(ok) {
    gen rdiff_P = 100 * (tilde_r_P - r_star_P)
    gen rdiff_A = 100 * (tilde_r_A - r_star_A)
    twoway ///
        (line rdiff_P year if `win', lcolor("`P_d'") lwidth(`lw')) ///
        (line rdiff_A year if `win', lcolor("`A_d'") lwidth(`lw')), ///
        yline(0, lcolor(gs8) lwidth(thin)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("Adjusted return minus R* (p.p.)", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(1 "R{sub:P} - R{sup:*}" 2 "R{sub:A} - R{sup:*}") ring(1) position(6) cols(2) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(mac_rdiff, replace)
    graph export "`out'/${PFX}_rate_differential.pdf", as(pdf) name(mac_rdiff) replace 
	*width(2006)
    _chkfile "`out'/${PFX}_rate_differential.pdf"
}
else di as error "*** 06_macro: M1 PULADO (faltam tilde_r / r_star) ***"
*Interest Rate Differential vs. Rest of the World

* ============================================================================
* M2 — NFA trajectories (levels; set NFA_AS_SHARE = 1 for % of GDP)
* ============================================================================
local NFA_AS_SHARE 0
_hasall f_P f_A y_P y_A
if r(ok) {
    if `NFA_AS_SHARE' {
        gen nfa_P = 100 * f_P / y_P
        gen nfa_A = 100 * f_A / y_A
        local nfa_lab "Net Foreign Assets (% of GDP)"
    }
    else {
        gen nfa_P = f_P
        gen nfa_A = f_A
        local nfa_lab "Net Foreign Assets (Efficiency Units)"
    }
    twoway ///
        (line nfa_P year if `win', lcolor("`P_d'") lwidth(`lw')) ///
        (line nfa_A year if `win', lcolor("`A_d'") lwidth(`lw')), ///
        yline(0, lcolor(gs8) lwidth(thin)) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("`nfa_lab'", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(1 "F{sub:P}" 2 "F{sub:A}") ring(1) position(6) cols(2) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(mac_nfa, replace)
    graph export "`out'/${PFX}_nfa.pdf", as(pdf) name(mac_nfa) replace 
	*width(2006)
    _chkfile "`out'/${PFX}_nfa.pdf"
}
else di as error "*** 06_macro: M2 PULADO (faltam f / y) ***"
*Net Foreign Asset Position

* ============================================================================
* M3 — INDIVIDUAL net wage (dark, solid) and pension (light, dashed)
* ============================================================================
_hasall w_net_indiv_P w_net_indiv_A e_indiv_P e_indiv_A
if r(ok) {
    twoway ///
        (line w_net_indiv_P year if `win', lcolor("`P_d'") lwidth(`lw')) ///
        (line e_indiv_P     year if `win', lcolor("`P_l'") lwidth(`lw') lpattern(dash)) ///
        (line w_net_indiv_A year if `win', lcolor("`A_d'") lwidth(`lw')) ///
        (line e_indiv_A     year if `win', lcolor("`A_l'") lwidth(`lw') lpattern(dash)), ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("Efficiency Units", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(1 "W{sub:P} - T{sub:P} (net wage)" 2 "E{sub:P} (pensions)" ///
                     3 "W{sub:A} - T{sub:A} (net wage)" 4 "E{sub:A} (pensions)") ///
               ring(1) position(6) cols(4) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(mac_wage_pens, replace)
    graph export "`out'/${PFX}_wage_pension_indiv.pdf", as(pdf) name(mac_wage_pens) replace 
	*width(2006)
    _chkfile "`out'/${PFX}_wage_pension_indiv.pdf"
}
else di as error "*** 06_macro: M3 PULADO (faltam w_indiv / t_indiv / e_indiv) ***"
*Individual Net Wage and Pension

* ============================================================================
* M4 — AGGREGATE tax burden (dark, solid) and pension spending (light, dashed), % of GDP
* ============================================================================
_hasall e_P e_A t_P t_A y_P y_A
if r(ok) {
    gen tax_y_P  = 100 * t_P / y_P
    gen tax_y_A  = 100 * t_A / y_A
    gen pens_y_P = 100 * e_P / y_P
    gen pens_y_A = 100 * e_A / y_A
    twoway ///
        (line tax_y_P  year if `win', lcolor("`P_d'") lwidth(`lw')) ///
        (line pens_y_P year if `win', lcolor("`P_l'") lwidth(`lw') lpattern(dash)) ///
        (line tax_y_A  year if `win', lcolor("`A_d'") lwidth(`lw')) ///
        (line pens_y_A year if `win', lcolor("`A_l'") lwidth(`lw') lpattern(dash)), ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("Year", size(medsmall)) ///
        ytitle("% of GDP", size(medsmall)) ///
        title(" ", size(medium)) ///
        legend(order(1 "T{sub:P} (taxes)" 2 "E{sub:P} (pensions)" ///
                     3 "T{sub:A} (taxes)" 4 "E{sub:A} (pensions)") ///
               ring(1) position(6) cols(4) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(mac_fiscal, replace)
    graph export "`out'/${PFX}_pensions_taxes_gdp.pdf", as(pdf) name(mac_fiscal) replace 
	*width(2006)
    _chkfile "`out'/${PFX}_pensions_taxes_gdp.pdf"
}
else di as error "*** 06_macro: M4 PULADO (faltam e / t / y agregados) ***"

*"Tax Burden and Pension Spending"

* ---- endpoint readout ($YEAR0 -> $J_END) ------------------------------------
di as text _n "===== MACRO ANALYSIS (${PFX}) — $YEAR0 -> $J_END ====="
foreach v in rdiff_P rdiff_A nfa_P nfa_A w_net_indiv_P w_net_indiv_A ///
             e_indiv_P e_indiv_A tax_y_P tax_y_A pens_y_P pens_y_A {
    capture confirm variable `v'
    if _rc continue
    quietly summ `v' if year == $YEAR0, meanonly
    local v0 = r(mean)
    quietly summ `v' if year == $J_END, meanonly
    di as text "  " %-14s "`v'" " : " as result %8.3f `v0' as text "  ->  " as result %8.3f r(mean)
}

* ---- persist (compact) ------------------------------------------------------
char _dta[pfx] "$PFX"
save "Stata_Cache/${PFX}_macro.dta", replace
di as text "06_macro done -> `out'/ + Stata_Cache/${PFX}_macro.dta"
