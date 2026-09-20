* ============================================================================
* 07.1_DECOMP — Stata_Cache/${PFX}_{P,A}_nfa.dta (+ y_x, f_x from _pos caches)
*             -> Plots/macro_analysis/
*   D1) Aggregate effect flows (sum over cohorts), P | A            [main text]
*       bands: dilution | net returns | savings   (turnover nets to zero -> omitted)
*   D2) Per-cohort effect flows, P | A, one figure per cohort       [appendix]
*       bands: turnover | dilution | net returns | savings
*       D1/D2 are drawn as ONE by(country) graph each (shared y-axis, one legend).
*   D3.1) P - A of aggregate effect flows (annual), single panel
*   D3.2) P - A of annual flows incl. -Delta k, -Delta b (sum = Delta f)
*   D4)   P - A of NFA contributions (cumulative since $YEAR0), single panel
* Sign convention: P - A  (>0 = channel favours P's accumulation relative to A)
* Runs AFTER 05_nfa (needs its caches). No $TAG. $PFX default "base".
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
local bw  0.8
local YCOMMON 1            // P|A panels share the y-axis
local DIFF_AS_SHARE 0      // D3.1 only: 1 = each country's flows in % of own GDP before differencing
local yr = cond(`YCOMMON', "", "yrescale")

* by() layout shared by every P|A figure
local byopts `"cols(2) note("") legend(position(6) ring(1)) iscale(1) graphregion(color(white)) `yr'"'

* ---- effect palette (kept distinct from the cohort blue/orange/green) -------
local c_chn "160 160 160"
local c_dil "131 197 190"
local c_ret "0 109 119"
local c_act "226 149 120"
* G1 contributions (RGB of the named colours used in 05_nfa)
local c_mb  "193 5 52"
local c_ay  "26 71 111"
local c_am  "227 126 35"
local c_ao  "85 117 47"
local c_mk  "128 128 128"

* ---- helpers ---------------------------------------------------------------
capture program drop _addband
program define _addband
    args name
    gen `name'_bot = cond(`name' >= 0, _pos, _neg)
    gen `name'_top = `name'_bot + `name'
    replace _pos = _pos + max(`name', 0)
    replace _neg = _neg + min(`name', 0)
end

capture program drop _chkfile
program define _chkfile
    args f
    capture confirm file "`f'"
    if _rc di as error "*** export FALHOU: `f' ***"
    else   di as text  "  -> `f' OK"
end

* ============================================================================
* LOAD: both caches, aggregate over cohorts, suffix by country, merge on year
* ============================================================================
foreach e in P A {
    foreach f in nfa pos {
        capture confirm file "Stata_Cache/${PFX}_`e'_`f'.dta"
        if _rc {
            di as error "*** 07_decomp: Stata_Cache/${PFX}_`e'_`f'.dta AUSENTE — rode 05_nfa antes ***"
            exit 601
        }
    }
    use "Stata_Cache/${PFX}_`e'_nfa.dta", clear
    capture confirm variable chn_y
    if _rc {
        di as error "*** 07_decomp: ${PFX}_`e'_nfa.dta sem fluxos (05_nfa pulou G2) ***"
        exit 601
    }
    merge 1:1 year using "Stata_Cache/${PFX}_`e'_pos.dta", keepusing(y_x f_x) keep(1 3) nogen
    sort year

    * aggregate effects (sum over cohorts); turnover must net to zero
    foreach ef in chn dil ret act {
        gen `ef'_agg = `ef'_y + `ef'_m + `ef'_o
    }
    capture assert abs(chn_agg) < 1e-8
    if _rc {
        quietly summ chn_agg
        di as error "07_decomp[`e']: AVISO — turnover agregado != 0 (max |.| = " %9.2e max(abs(r(min)), abs(r(max))) ")"
    }

    * per-cohort annual total = sum of its four effects (identity)
    foreach c in y m o {
        gen d1_`c' = chn_`c' + dil_`c' + ret_`c' + act_`c'
    }
    * annual Delta f from the cumulative series
    gen d1_f = d_f - d_f[_n-1]
    replace d1_f = 0 in 1

    foreach v of varlist _all {
        if "`v'" != "year" rename `v' `v'_`e'
    }
    tempfile d_`e'
    save `d_`e''
}
use `d_P', clear
merge 1:1 year using `d_A'
capture assert _merge == 3
if _rc di as error "07_decomp: AVISO — anos nao coincidem entre P e A; uso a intersecao"
keep if _merge == 3
drop _merge
sort year

* ============================================================================
* D1 / D2 — by(country) figures on a stacked (long) copy of the data
* ============================================================================
tempfile longP
preserve
    keep year *_P
    foreach v of varlist *_P {
        local nv = substr("`v'", 1, length("`v'") - 2)
        rename `v' `nv'
    }
    gen byte country = 1
    save `longP'
restore
preserve
    keep year *_A
    foreach v of varlist *_A {
        local nv = substr("`v'", 1, length("`v'") - 2)
        rename `v' `nv'
    }
    gen byte country = 2
    append using `longP'
    label define cty 1 "P" 2 "A", replace
    label values country cty
    sort country year
    gen _pos = 0
    gen _neg = 0

    * ---- D1: aggregate effect flows ----------------------------------------
    foreach ef in dil ret act {
        _addband `ef'_agg
    }
    twoway ///
        (rbar dil_agg_bot dil_agg_top year if `win', barwidth(`bw') color("`c_dil'") lwidth(none)) ///
        (rbar ret_agg_bot ret_agg_top year if `win', barwidth(`bw') color("`c_ret'") lwidth(none)) ///
        (rbar act_agg_bot act_agg_top year if `win', barwidth(`bw') color("`c_act'") lwidth(none)) ///
        (line d1_a year if `win', lcolor(black) lwidth(medthick)), ///
        by(country, `byopts') ///
        xsize(11) ysize(5) ///
        yline(0, lcolor(gs6) lwidth(thin)) ///
        subtitle(, size(medium) nobox) ///
        xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
        ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
        xtitle("", size(medsmall)) ///
        ytitle("Annual flow (efficiency units)", size(medsmall)) ///
        legend(order(1 "Dilution" 2 "Net Returns" 3 "Savings" 4 "{&Delta}a (total)") ///
               cols(4) size(medsmall) region(lcolor(white))) ///
        graphregion(color(white)) plotregion(color(white)) name(dc_agg, replace)
    graph export "`out'/${PFX}_hh_assets_agg_flows.pdf", as(pdf) name(dc_agg) replace
    _chkfile "`out'/${PFX}_hh_assets_agg_flows.pdf"
    drop *_bot *_top

    * ---- D2: per-cohort effect flows, one figure per cohort ----------------
    local lab_y "Young"
    local lab_m "Middle-Aged"
    local lab_o "Retirees"
    foreach c in y m o {
        replace _pos = 0
        replace _neg = 0
        foreach ef in chn dil ret act {
            _addband `ef'_`c'
        }
        twoway ///
            (rbar chn_`c'_bot chn_`c'_top year if `win', barwidth(`bw') color("`c_chn'") lwidth(none)) ///
            (rbar dil_`c'_bot dil_`c'_top year if `win', barwidth(`bw') color("`c_dil'") lwidth(none)) ///
            (rbar ret_`c'_bot ret_`c'_top year if `win', barwidth(`bw') color("`c_ret'") lwidth(none)) ///
            (rbar act_`c'_bot act_`c'_top year if `win', barwidth(`bw') color("`c_act'") lwidth(none)) ///
            (line d1_`c' year if `win', lcolor(black) lwidth(medthick)), ///
            by(country, `byopts') ///
        xsize(11) ysize(5) ///
            yline(0, lcolor(gs6) lwidth(thin)) ///
            subtitle(, size(medium) nobox) ///
            xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
            ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
            xtitle("", size(medsmall)) ///
            ytitle("Annual flow (efficiency units)", size(medsmall)) ///
            legend(order(1 "Turnover" 2 "Dilution" 3 "Net Returns" 4 "Savings" 5 "{&Delta}({&eta}{sup:`c'}a)") ///
                   cols(5) size(medsmall) region(lcolor(white))) ///
            graphregion(color(white)) plotregion(color(white)) name(dc_`c', replace)
        graph export "`out'/${PFX}_hh_assets_flows_`c'.pdf", as(pdf) name(dc_`c') replace
        _chkfile "`out'/${PFX}_hh_assets_flows_`c'.pdf"
        drop *_bot *_top
    }
restore

* wide data is back from here on
gen _pos = 0
gen _neg = 0

* ============================================================================
* D3.1 — P - A of aggregate effect flows (annual), single panel
* ============================================================================
if `DIFF_AS_SHARE' {
    gen sc_P = 100 / y_x_P
    gen sc_A = 100 / y_x_A
    local dunit "% of own GDP"
}
else {
    gen sc_P = 1
    gen sc_A = 1
    local dunit "efficiency units"
}
foreach ef in dil ret act {
    gen dif_`ef' = `ef'_agg_P*sc_P - `ef'_agg_A*sc_A
}
gen dif_a = d1_a_P*sc_P - d1_a_A*sc_A
gen dif_f = d1_f_P*sc_P - d1_f_A*sc_A
capture assert abs(dif_dil + dif_ret + dif_act - dif_a) < 1e-4 if `win'
if _rc di as error "07_decomp: AVISO — barras de D3 nao somam a dif_a"

replace _pos = 0
replace _neg = 0
foreach ef in dil ret act {
    _addband dif_`ef'
}
twoway ///
    (rbar dif_dil_bot dif_dil_top year if `win', barwidth(`bw') color("`c_dil'") lwidth(none)) ///
    (rbar dif_ret_bot dif_ret_top year if `win', barwidth(`bw') color("`c_ret'") lwidth(none)) ///
    (rbar dif_act_bot dif_act_top year if `win', barwidth(`bw') color("`c_act'") lwidth(none)) ///
    (function y = 0, range($YEAR0 $J_END) lcolor(gs6) lwidth(thin)) ///
    (line dif_a year if `win', lcolor(black) lwidth(medthick)) ///
    (line dif_f year if `win', lcolor(black) lwidth(medthick) lpattern(dash)), ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
    xtitle("", size(medsmall)) ///
    ytitle("Annual flow, P - A (`dunit')", size(medsmall)) ///
    title(" ", size(medium)) ///
    legend(order(1 "Dilution" 2 "Net Returns" 3 "Savings" ///
                 5 "{&Delta}a{sub:P} - {&Delta}a{sub:A}" 6 "{&Delta}f{sub:P} - {&Delta}f{sub:A}") ///
           ring(1) position(6) cols(5) size(medsmall) region(lcolor(white))) ///
    xsize(9) ysize(4.5) ///
    graphregion(color(white)) plotregion(color(white)) name(dc_diff_flows, replace)
graph export "`out'/${PFX}_diff_effects_flows_1.pdf", as(pdf) name(dc_diff_flows) replace
_chkfile "`out'/${PFX}_diff_effects_flows_1.pdf"
drop *_bot *_top

* ============================================================================
* D3.2 — P - A of annual flows, single panel
*      bars: dilution | net returns | savings | -Delta k | -Delta b   (sum = Delta f)
*      line: Delta f_P - Delta f_A   (dashed)
* ============================================================================
foreach e in P A {
    gen fl_mk_`e' = cA_mk_`e' - cA_mk_`e'[_n-1]     // -Delta k (already signed in cA_mk)
    gen fl_mb_`e' = cA_mb_`e' - cA_mb_`e'[_n-1]     // -Delta b
    replace fl_mk_`e' = 0 in 1
    replace fl_mb_`e' = 0 in 1
}
foreach ef in dil ret act {
    capture drop dif_`ef'
    gen dif_`ef' = `ef'_agg_P - `ef'_agg_A
}
capture drop dif_mk
gen dif_mk = fl_mk_P - fl_mk_A
capture drop dif_mb
gen dif_mb = fl_mb_P - fl_mb_A
capture drop dif_f
gen dif_f  = d1_f_P - d1_f_A
capture assert abs(dif_dil + dif_ret + dif_act + dif_mk + dif_mb - dif_f) < 1e-4 if `win'
if _rc di as error "07_decomp: AVISO — barras de D3.2 nao somam a dif_f"

* boring greys for the two accounting counterparts
local c_mk_boring "150 150 150"
local c_mb_boring "90 90 90"

replace _pos = 0
replace _neg = 0
foreach ef in dil ret act mk mb {
    _addband dif_`ef'
}
twoway ///
    (rbar dif_dil_bot dif_dil_top year if `win', barwidth(`bw') color("`c_dil'") lwidth(none)) ///
    (rbar dif_ret_bot dif_ret_top year if `win', barwidth(`bw') color("`c_ret'") lwidth(none)) ///
    (rbar dif_act_bot dif_act_top year if `win', barwidth(`bw') color("`c_act'") lwidth(none)) ///
    (rbar dif_mk_bot  dif_mk_top  year if `win', barwidth(`bw') color("`c_mk_boring'") lwidth(none)) ///
    (rbar dif_mb_bot  dif_mb_top  year if `win', barwidth(`bw') color("`c_mb_boring'") lwidth(none)) ///
    (function y = 0, range($YEAR0 $J_END) lcolor(gs6) lwidth(thin)) ///
    (line dif_f year if `win', lcolor(black) lwidth(medthick) lpattern(dash)), ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
    xtitle("", size(medsmall)) ///
    ytitle("Annual flow, P - A (efficiency units)", size(medsmall)) ///
    title(" ", size(medium)) ///
    legend(order(1 "Dilution" 2 "Net Returns" 3 "Savings" ///
                 4 "-{&Delta}k (capital)" 5 "-{&Delta}b (debt)" ///
                 7 "{&Delta}f{sub:P} - {&Delta}f{sub:A}") ///
           ring(1) position(6) cols(6) size(medsmall) region(lcolor(white))) ///
    xsize(9) ysize(4.5) ///
    graphregion(color(white)) plotregion(color(white)) name(dc_diff_flows2, replace)
graph export "`out'/${PFX}_diff_effects_flows_2.pdf", as(pdf) name(dc_diff_flows2) replace
_chkfile "`out'/${PFX}_diff_effects_flows_2.pdf"
drop *_bot *_top

* ============================================================================
* D4 — P - A of NFA contributions (cumulative since $YEAR0), single panel
*      always in efficiency units (cumulative changes over current GDP are odd)
* ============================================================================
foreach s in mb ay am ao mk {
    gen dif_cA_`s' = cA_`s'_P - cA_`s'_A
}
gen dif_df = d_f_P - d_f_A
capture assert abs(dif_cA_mb+dif_cA_ay+dif_cA_am+dif_cA_ao+dif_cA_mk - dif_df) < 1e-5 if `win'
if _rc di as error "07_decomp: AVISO — barras de D4 nao somam a dif_df"

replace _pos = 0
replace _neg = 0
foreach s in mb ay am ao mk {
    _addband dif_cA_`s'
}
twoway ///
    (rbar dif_cA_mb_bot dif_cA_mb_top year if `win', barwidth(`bw') color("`c_mb'") lwidth(none)) ///
    (rbar dif_cA_ay_bot dif_cA_ay_top year if `win', barwidth(`bw') color("`c_ay'") lwidth(none)) ///
    (rbar dif_cA_am_bot dif_cA_am_top year if `win', barwidth(`bw') color("`c_am'") lwidth(none)) ///
    (rbar dif_cA_ao_bot dif_cA_ao_top year if `win', barwidth(`bw') color("`c_ao'") lwidth(none)) ///
    (rbar dif_cA_mk_bot dif_cA_mk_top year if `win', barwidth(`bw') color("`c_mk'") lwidth(none)) ///
    (function y = 0, range($YEAR0 $J_END) lcolor(gs6) lwidth(thin)) ///
    (line dif_df year if `win', lcolor(black) lwidth(medthick)), ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
    xtitle("", size(medsmall)) ///
    ytitle("Change since $YEAR0, P - A (efficiency units)", size(medsmall)) ///
    title(" ", size(medium)) ///
    legend(order(1 "-{&Delta}b (debt)" 2 "{&Delta}a{sub:y} (young)" 3 "{&Delta}a{sub:m} (middle)" ///
                 4 "{&Delta}a{sub:o} (old)" 5 "-{&Delta}k (capital)" 7 "{&Delta}f{sub:P} - {&Delta}f{sub:A}") ///
           ring(1) position(6) cols(3) size(medsmall) region(lcolor(white))) ///
    xsize(9) ysize(4.5) ///
    graphregion(color(white)) plotregion(color(white)) name(dc_diff_nfa, replace)
graph export "`out'/${PFX}_diff_nfa_contributions.pdf", as(pdf) name(dc_diff_nfa) replace
_chkfile "`out'/${PFX}_diff_nfa_contributions.pdf"
drop *_bot *_top _pos _neg

* ---- channel weights at key years ----
foreach v in dif_dil dif_ret dif_act {
    quietly summ `v' if inrange(year,1970,2030), meanonly
    di "`v': cum = " r(sum)
}
quietly summ dif_f if inrange(year,1970,2030), meanonly
di "dif_f cum = " r(sum)

foreach v in dif_dil dif_ret dif_act {
    quietly summ `v' if inrange(year,2030,2050), meanonly
    di "`v': cum = " r(sum)
}
quietly summ dif_f if inrange(year,2030,2050), meanonly
di "dif_f cum = " r(sum)


foreach v in dif_dil dif_ret dif_act {
    quietly summ `v' if inrange(year,2050,2070), meanonly
    di "`v': cum = " r(sum)
}
quietly summ dif_f if inrange(year,2050,2070), meanonly
di "dif_f cum = " r(sum)


* ---- persist (compact) --------------------------------------------------------
keep year *_agg_P *_agg_A d1_a_P d1_a_A d1_f_P d1_f_A dif_* d_f_P d_f_A
char _dta[pfx] "$PFX"
save "Stata_Cache/${PFX}_decomp.dta", replace
di as text "07_decomp[${PFX}] done -> `out'/"


* ============================================================================
* CHANNEL WEIGHTS BY PHASE — cumulated share of (Delta f_P - Delta f_A)
*   share_j = (sum_t channel_j) / (sum_t dif_f)   over each window
*   >100% / negative shares are normal (offsetting channels); flagged as such.
* ============================================================================
di as text _n "{hline 64}"
di as text "  CHANNEL DECOMPOSITION OF (Delta f_P - Delta f_A), BY PHASE"
di as text "{hline 64}"

local windows "1970 2030 | 2030 2050 | 2050 2070"
local labs    "dif_act Savings | dif_ret Net_returns | dif_dil Dilution | dif_mk Capital | dif_mb Debt"

foreach w in "1970 2030" "2030 2050" "2050 2070" {
    tokenize `w'
    local y0 `1'
    local y1 `2'
    quietly summ dif_f if inrange(year,`y0',`y1'), meanonly
    local tot = r(sum)

    di as text _n "  Phase `y0'-`y1'   (cumulated {&Delta}f_P-{&Delta}f_A = " %6.3f `tot' ")"
    di as text "  " %-14s "channel" "   " %10s "cum. level" "   " %8s "share"
    di as text "  " "{hline 40}"

    * order: savings, net returns, dilution, capital, debt
    foreach pair in "dif_act Savings" "dif_ret Net-returns" "dif_dil Dilution" "dif_mk Capital" "dif_mb Debt" {
        tokenize "`pair'"
        local v    `1'
        local name `2'
        quietly summ `v' if inrange(year,`y0',`y1'), meanonly
        local lvl = r(sum)
        local pct = 100*`lvl'/`tot'
        di as text "  " %-14s "`name'" "   " as result %10.4f `lvl' ///
           as text "   " as result %7.1f `pct' as text "%"
    }
}
di as text _n "  Note: shares are cumulated channel / cumulated dif_f per window."
di


* ============================================================================
* 07.2_SHARES — Stata_Cache/${PFX}_{P,A}_pos.dta -> Plots/macro_analysis/
*   S1) Population shares:   tilde_s^y, tilde_s^m, tilde_s^o  (P & A -> 6 lines)
*   S2) Dependency ratios:   psi^y, psi^o                     (P & A -> 4 lines)
* Two DIFFERENT panels side by side -> by() does not apply; graph combine kept.
* ============================================================================
local lw medthick

* P orange tones (y light -> o dark), A purple tones
local P_y "245 180 120"
local P_m "214 118 45"
local P_o "140 70 20"
local A_y "190 170 220"
local A_m "128 90 175"
local A_o "70 40 110"

foreach e in P A {
    capture confirm file "Stata_Cache/${PFX}_`e'_pos.dta"
    if _rc {
        di as error "*** 08_shares: Stata_Cache/${PFX}_`e'_pos.dta AUSENTE ***"
        exit 601
    }
    use "Stata_Cache/${PFX}_`e'_pos.dta", clear
    sort year
    foreach c in y m o {
        capture confirm variable tilde_s_`c'_x
        if _rc {
            capture confirm variable share_`c'_x
            if !_rc rename share_`c'_x tilde_s_`c'_x
            else di as error "08_shares[`e']: sem tilde_s_`c'_x nem share_`c'_x"
        }
    }
    keep year tilde_s_y_x tilde_s_m_x tilde_s_o_x psi_y_x psi_o_x
    rename *_x *_`e'
    tempfile s_`e'
    save `s_`e''
}
use `s_P', clear
merge 1:1 year using `s_A', nogen
sort year

twoway ///
    (line tilde_s_y_P year if `win', lcolor("`P_y'") lwidth(`lw')) ///
    (line tilde_s_m_P year if `win', lcolor("`P_m'") lwidth(`lw')) ///
    (line tilde_s_o_P year if `win', lcolor("`P_o'") lwidth(`lw')) ///
    (line tilde_s_y_A year if `win', lcolor("`A_y'") lwidth(`lw')) ///
    (line tilde_s_m_A year if `win', lcolor("`A_m'") lwidth(`lw')) ///
    (line tilde_s_o_A year if `win', lcolor("`A_o'") lwidth(`lw')), ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%3.1f) labsize(medsmall)) ///
    xtitle("", size(medsmall)) ytitle("Share of population", size(medsmall)) ///
    title("Cohort Shares", size(medium)) ///
    legend(order(1 "s̃{sup:y}{sub:P}" 2 "s̃{sup:m}{sub:P}" 3 "s̃{sup:o}{sub:P}" ///
                 4 "s̃{sup:y}{sub:A}" 5 "s̃{sup:m}{sub:A}" 6 "s̃{sup:o}{sub:A}") ///
           ring(1) position(6) cols(6) size(medsmall) region(lcolor(white))) ///
    graphregion(color(white)) plotregion(color(white)) name(sh_shares, replace)

twoway ///
    (line psi_y_P year if `win', lcolor("`P_y'") lwidth(`lw')) ///
    (line psi_o_P year if `win', lcolor("`P_o'") lwidth(`lw')) ///
    (line psi_y_A year if `win', lcolor("`A_y'") lwidth(`lw')) ///
    (line psi_o_A year if `win', lcolor("`A_o'") lwidth(`lw')), ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%3.1f) labsize(medsmall)) ///
    xtitle("", size(medsmall)) ytitle("{&psi}", size(medsmall)) ///
    title("Dependency Ratios", size(medium)) ///
    legend(order(1 "{&psi}{sup:y}{sub:P}" 2 "{&psi}{sup:o}{sub:P}" ///
                 3 "{&psi}{sup:y}{sub:A}" 4 "{&psi}{sup:o}{sub:A}") ///
           ring(1) position(6) cols(4) size(medsmall) region(lcolor(white))) ///
    graphregion(color(white)) plotregion(color(white)) name(sh_psi, replace)

graph combine sh_shares sh_psi, cols(2) ///
    iscale(1) xsize(12) ysize(5) graphregion(color(white)) name(sh_comb, replace)
graph export "`out'/${PFX}_shares_psi.pdf", as(pdf) name(sh_comb) replace
_chkfile "`out'/${PFX}_shares_psi.pdf"

di as text "08_shares[${PFX}] done -> `out'/"

