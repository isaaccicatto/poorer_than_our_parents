* ============================================================================
* 05_NFA — Stata_Cache/${PFX}_P_pos.dta + Stata_Cache/${PFX}_A_pos.dta
*          -> asset decompositions, P | A in ONE by(country) graph
*             (shared y-axis, one legend, one font), into Plots/macro_analysis/
*   G1)  Decomposition of Variations in the Net Asset Position
*        Delta f = -Delta b + Delta a_y + Delta a_m + Delta a_o - Delta k
*   G2)  Decomposition of Household Assets (cumulated since $YEAR0)
*        per cohort, 4 effects: turnover | dilution | net returns | savings
*          turnover    = A_in - a(t-1)               (cohort flows y->m->o)
*          dilution    = -((g-1)/g) * A_in           (population/productivity growth)
*          net returns =  ((R_tilde-1)/g) * A_in     (gross return)
*          savings     =  budget residual            (-c_y | w-t-c_m | e-c_o)
*        dilution + net returns = (q-1) * A_in, q = R_tilde/g
*   G2b) same, raw annual flows
* Loop: per country, build series + save Stata_Cache/${PFX}_`e'_nfa.dta.
* Then: stack P (country=1) and A (country=2) and plot with by(country).
* Runs AFTER the scenario loop in 00_master. No $TAG. $PFX default "base".
* Identities are capture'd asserts: warn, don't abort.
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
foreach f in nfa_decomp_variations hh_assets_decomp hh_assets_decomp_flows {
    cap erase "`out'/${PFX}_`f'.pdf"
}
local win "inrange(year, $YEAR0, $J_END)"
local bw  0.8
local YCOMMON 1        // 1 = same y-axis for P and A; 0 = free axes (yrescale)
local yr = cond(`YCOMMON', "", "yrescale")

* by() layout shared by every P|A figure
local byopts `"cols(2) note("") legend(position(6) ring(1)) iscale(1) graphregion(color(white)) `yr'"'

* ---- palette: 4 tones per cohort, light -> dark = turnover, dilution, net returns, savings
local c_y_1 "200 215 235"
local c_y_2 "150 180 215"
local c_y_3 "90 130 185"
local c_y_4 "44 87 148"
local c_m_1 "250 215 185"
local c_m_2 "240 175 130"
local c_m_3 "220 135 80"
local c_m_4 "203 106 43"
local c_o_1 "200 225 205"
local c_o_2 "150 195 160"
local c_o_3 "105 165 120"
local c_o_4 "63 124 72"

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
* PER-COUNTRY SERIES -> CACHES (no plotting here)
* ============================================================================
local ok_P = 0
local ok_A = 0
local hf_P = 0
local hf_A = 0

foreach e in P A {
    di as text _n "{hline 60}"
    di as text "  05_nfa   ${PFX}_`e'"
    di as text "{hline 60}"

    capture confirm file "Stata_Cache/${PFX}_`e'_pos.dta"
    if _rc {
        di as error "*** 05_nfa: Stata_Cache/${PFX}_`e'_pos.dta AUSENTE — painel `e' PULADO ***"
        continue
    }
    use "Stata_Cache/${PFX}_`e'_pos.dta", clear
    sort year
    capture assert year[1] == $YEAR0
    if _rc di as error "05_nfa[`e']: AVISO — primeiro ano != $YEAR0"

    * ---- nu/omega: char -> column -> default --------------------------------
    local nu    : char _dta[nu]
    local omega : char _dta[omega]
    if "`nu'" == "" {
        capture confirm variable nu
        if !_rc {
            quietly summ nu in 1, meanonly
            local nu = r(mean)
        }
    }
    if "`omega'" == "" {
        capture confirm variable omega
        if !_rc {
            quietly summ omega in 1, meanonly
            local omega = r(mean)
        }
    }
    if "`nu'" == ""    local nu    = 0.8
    if "`omega'" == "" local omega = 0.975
    di as text "05_nfa[`e']: nu = `nu'   omega = `omega'"

    * ---- resolve columns (aggregates preferred; reconstruct if missing) ------
    foreach c in y m o {
        capture confirm variable tilde_s_`c'_x
        if _rc {
            capture confirm variable share_`c'_x
            if !_rc rename share_`c'_x tilde_s_`c'_x
        }
    }
    local faltam ""
    foreach v in f_x k_x y_x a_y_indiv_x a_m_indiv_x a_o_indiv_x ///
                 tilde_s_y_x tilde_s_m_x tilde_s_o_x {
        capture confirm variable `v'
        if _rc local faltam "`faltam' `v'"
    }
    if "`faltam'" != "" {
        di as error "*** 05_nfa[`e']: colunas AUSENTES:`faltam' — painel `e' PULADO ***"
        continue
    }
    foreach c in y m o {
        capture confirm variable a_`c'_x
        if _rc gen a_`c'_x = a_`c'_indiv_x * tilde_s_`c'_x
    }
    capture confirm variable b_x
    if _rc {
        capture confirm variable b_y_x
        if !_rc gen b_x = b_y_x * y_x
        else {
            capture confirm variable b_y
            if !_rc gen b_x = b_y * y_x
            else {
                di as error "*** 05_nfa[`e']: nem b_x, b_y_x nem b_y — painel `e' PULADO ***"
                continue
            }
        }
    }

    * ---- q / g / tilde_r: any two give the third  (q = tilde_r / g) ---------
    capture confirm variable q_x
    local hq = (_rc == 0)
    capture confirm variable g_x
    local hg = (_rc == 0)
    capture confirm variable tilde_r_x
    local hr = (_rc == 0)
    if !`hq' & `hg' & `hr' {
        gen q_x = tilde_r_x / g_x
        di as text "05_nfa[`e']: q_x reconstruido de tilde_r_x/g_x"
    }
    if !`hg' & `hq' & `hr' {
        gen g_x = tilde_r_x / q_x
        di as text "05_nfa[`e']: g_x reconstruido de tilde_r_x/q_x"
    }
    if `hq' & `hg' & `hr' {
        capture assert abs(q_x - tilde_r_x/g_x) < 1e-8
        if _rc {
            di as error "05_nfa[`e']: AVISO — q_x != tilde_r_x/g_x no CSV (timing?). Mantenho q_x e refaco g_x = tilde_r_x/q_x"
            replace g_x = tilde_r_x / q_x
        }
    }

    * ---- G1 series: contributions to Delta f (relative to $YEAR0) ------------
    foreach v in a_y_x a_m_x a_o_x k_x b_x f_x {
        gen d_`v' = `v' - `v'[1]
    }
    gen cA_mb = -d_b_x
    gen cA_ay =  d_a_y_x
    gen cA_am =  d_a_m_x
    gen cA_ao =  d_a_o_x
    gen cA_mk = -d_k_x
    gen d_f   =  d_f_x
    capture assert abs(cA_mb+cA_ay+cA_am+cA_ao+cA_mk - d_f) < 1e-6 if !missing(d_f)
    if _rc {
        quietly gen double __dev = abs(cA_mb+cA_ay+cA_am+cA_ao+cA_mk - d_f)
        quietly summ __dev
        di as error "05_nfa[`e']: AVISO identidade G1 viola tol (max dev = " %9.2e r(max) ") — prossigo"
        drop __dev
    }
    gen d_a = d_a_y_x + d_a_m_x + d_a_o_x

    * Shapley pair (cache/readout only, not plotted)
    foreach c in y m o {
        gen sav_`c'  = 0.5*(tilde_s_`c'_x[1] + tilde_s_`c'_x) ///
                         *(a_`c'_indiv_x - a_`c'_indiv_x[1])
        gen size_`c' = 0.5*(a_`c'_indiv_x[1] + a_`c'_indiv_x) ///
                         *(tilde_s_`c'_x - tilde_s_`c'_x[1])
        capture assert abs(sav_`c' + size_`c' - d_a_`c'_x) < 1e-6 if !missing(d_a_`c'_x)
        if _rc di as error "05_nfa[`e']: AVISO Shapley `c' viola tol — prossigo"
    }

    * ---- G2 series: 4-effect flow accounting per cohort ----------------------
    local haveflows = 1
    local faltam ""
    foreach v in q_x g_x tilde_r_x e_x w_x t_x c_y_x c_m_x c_o_x {
        capture confirm variable `v'
        if _rc {
            local haveflows = 0
            local faltam "`faltam' `v'"
        }
    }
    if !`haveflows' di as error "*** 05_nfa[`e']: SEM colunas de fluxo:`faltam' -> G2/G2b PULADOS ***"
    else            di as text  "05_nfa[`e']: colunas de fluxo OK -> G2/G2b"

    if `haveflows' {
        gen Ain_o = a_o_x[_n-1] + (1-`omega')*a_m_x[_n-1]
        gen Ain_m = `omega'*a_m_x[_n-1] + (1-`nu')*a_y_x[_n-1]
        gen Ain_y = `nu'*a_y_x[_n-1]

        * turnover: net inflow from cohort transitions
        gen chn_o = Ain_o - a_o_x[_n-1]
        gen chn_m = Ain_m - a_m_x[_n-1]
        gen chn_y = Ain_y - a_y_x[_n-1]

        * (q-1)*Ain = net returns + dilution, split exactly
        foreach c in y m o {
            gen ret_`c' =  ((tilde_r_x - 1) / g_x) * Ain_`c'
            gen dil_`c' = -((g_x - 1) / g_x) * Ain_`c'
        }

        * active savings: cohort budget residual
        gen act_o = e_x - c_o_x
        gen act_m = w_x - t_x - c_m_x
        gen act_y = - c_y_x

        * step identity (warn, don't abort)
        foreach c in y m o {
            capture assert abs(chn_`c'+ret_`c'+dil_`c'+act_`c' - (a_`c'_x - a_`c'_x[_n-1])) < 1e-5 if _n > 1
            if _rc {
                quietly gen double __dev = abs(chn_`c'+ret_`c'+dil_`c'+act_`c' - (a_`c'_x - a_`c'_x[_n-1])) if _n > 1
                quietly summ __dev
                di as error "05_nfa[`e']: AVISO identidade step `c' (max dev = " %9.2e r(max) ") — prossigo"
                drop __dev
            }
            foreach ef in chn dil ret act {
                replace `ef'_`c' = 0 in 1
            }
        }
        gen d1_a = d_a - d_a[_n-1]
        replace d1_a = 0 in 1

        * cumulate from $YEAR0 (warn, don't abort)
        foreach c in y m o {
            foreach ef in chn dil ret act {
                gen c`ef'_`c' = sum(`ef'_`c')
            }
            capture assert abs(cchn_`c'+cdil_`c'+cret_`c'+cact_`c' - d_a_`c'_x) < 1e-4 if !missing(d_a_`c'_x)
            if _rc {
                quietly gen double __dev = abs(cchn_`c'+cdil_`c'+cret_`c'+cact_`c' - d_a_`c'_x)
                quietly summ __dev
                di as error "05_nfa[`e']: AVISO identidade cumulada `c' (max dev = " %9.2e r(max) ") — prossigo"
                drop __dev
            }
        }
    }

    * ---- endpoint readout at $J_END ------------------------------------------
    di as text _n "===== DECOMPOSITION (${PFX}_`e') — Delta $YEAR0 -> $J_END ====="
    local list1 "cA_mb cA_ay cA_am cA_ao cA_mk d_f d_a"
    if `haveflows' {
        foreach c in y m o {
            local list1 "`list1' cchn_`c' cdil_`c' cret_`c' cact_`c'"
        }
    }
    foreach v of local list1 {
        quietly summ `v' if year == $J_END, meanonly
        di as text "  " %-8s "`v'" " = " as result %7.3f r(mean)
    }

    * ---- persist the decomposition series (compact) --------------------------
    local keeplist "year cA_mb cA_ay cA_am cA_ao cA_mk d_f d_a size_y sav_y size_m sav_m size_o sav_o"
    if `haveflows' {
        foreach c in y m o {
            local keeplist "`keeplist' chn_`c' dil_`c' ret_`c' act_`c' cchn_`c' cdil_`c' cret_`c' cact_`c'"
        }
        local keeplist "`keeplist' d1_a"
    }
    keep `keeplist'
    char _dta[pfx]    "$PFX"
    char _dta[region] "`e'"
    save "Stata_Cache/${PFX}_`e'_nfa.dta", replace
    di as text "05_nfa[`e'] series -> Stata_Cache/${PFX}_`e'_nfa.dta"

    local ok_`e' = 1
    local hf_`e' = `haveflows'
}

* ============================================================================
* STACK P (country=1) and A (country=2) -> by(country) figures
* ============================================================================
if !(`ok_P' & `ok_A') {
    di as error "*** 05_nfa: painel P ou A ausente — nada plotado ***"
    exit
}
use "Stata_Cache/${PFX}_P_nfa.dta", clear
gen byte country = 1
append using "Stata_Cache/${PFX}_A_nfa.dta"
replace country = 2 if missing(country)
label define cty 1 "P" 2 "A", replace
label values country cty
sort country year
gen _pos = 0
gen _neg = 0

* ---- G1: 5 bands + line, by(country) ----------------------------------------
foreach s in cA_mb cA_ay cA_am cA_ao cA_mk {
    _addband `s'
}
twoway ///
 (rbar cA_mb_bot cA_mb_top year if `win', barwidth(`bw') color(cranberry)    lwidth(none)) ///
 (rbar cA_ay_bot cA_ay_top year if `win', barwidth(`bw') color(navy)         lwidth(none)) ///
 (rbar cA_am_bot cA_am_top year if `win', barwidth(`bw') color(dkorange)     lwidth(none)) ///
 (rbar cA_ao_bot cA_ao_top year if `win', barwidth(`bw') color(forest_green) lwidth(none)) ///
 (rbar cA_mk_bot cA_mk_top year if `win', barwidth(`bw') color(gs8)          lwidth(none)) ///
 (line d_f year if `win', lcolor(black) lwidth(medthick)), ///
    by(country, `byopts') ///
        xsize(11) ysize(5) ///
    yline(0, lcolor(gs6) lwidth(thin)) ///
    subtitle(, size(medium) nobox) ///
    xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
    ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
    xtitle("", size(medsmall)) ///
    ytitle("Change since $YEAR0 (efficiency units)", size(medsmall)) ///
    legend(order(1 "-{&Delta}b (debt)" 2 "{&Delta}a{sub:y} (young)" ///
                 3 "{&Delta}a{sub:m} (middle)" 4 "{&Delta}a{sub:o} (old)" ///
                 5 "-{&Delta}k (capital)" 6 "{&Delta}f (NFA)") ///
           cols(3) size(medsmall) region(lcolor(white))) ///
    graphregion(color(white)) plotregion(color(white)) name(nfa_G1, replace)
graph export "`out'/${PFX}_nfa_decomp_variations.pdf", as(pdf) name(nfa_G1) replace
_chkfile "`out'/${PFX}_nfa_decomp_variations.pdf"
drop *_bot *_top

* ---- G2 / G2b: 12 bands + line, by(country) ---------------------------------
if `hf_P' & `hf_A' {
    foreach mode in cum flow {
        if "`mode'" == "cum" {
            local vp   "c"
            local lv   "d_a"
            local gnm  "nfa_G2"
            local fnm  "hh_assets_decomp"
            local ytit "Change since $YEAR0 (efficiency units)"
            local llab "{&Delta}a (total)"
        }
        else {
            local vp   ""
            local lv   "d1_a"
            local gnm  "nfa_G2b"
            local fnm  "hh_assets_decomp_flows"
            local ytit "Annual flow (efficiency units)"
            local llab "Annual {&Delta}a"
        }
        replace _pos = 0
        replace _neg = 0
        local plots ""
        foreach c in y m o {
            local k = 1
            foreach ef in chn dil ret act {
                _addband `vp'`ef'_`c'
                local plots `"`plots' (rbar `vp'`ef'_`c'_bot `vp'`ef'_`c'_top year if `win', barwidth(`bw') color("`c_`c'_`k''") lwidth(none))"'
                local ++k
            }
        }
        twoway `plots' ///
            (line `lv' year if `win', lcolor(black) lwidth(medthick)), ///
            by(country, `byopts') ///
        xsize(11) ysize(5) ///
            yline(0, lcolor(gs6) lwidth(thin)) ///
            subtitle(, size(medium) nobox) ///
            xlabel($YEAR0(20)$J_END, labsize(medsmall)) ///
            ylabel(, format(%9.2f) labsize(medsmall) angle(0)) ///
            xtitle("", size(medsmall)) ///
            ytitle("`ytit'", size(medsmall)) ///
            legend(order(1 "Y: Turnover" 2 "Y: Dilution" 3 "Y: Net Returns" 4 "Y: Savings" ///
                         5 "M: Turnover" 6 "M: Dilution" 7 "M: Net Returns" 8 "M: Savings" ///
                         9 "O: Turnover" 10 "O: Dilution" 11 "O: Net Returns" 12 "O: Savings" ///
                         13 "`llab'") ///
                   cols(4) size(small) region(lcolor(white))) ///
            graphregion(color(white)) plotregion(color(white)) name(`gnm', replace)
        graph export "`out'/${PFX}_`fnm'.pdf", as(pdf) name(`gnm') replace
        _chkfile "`out'/${PFX}_`fnm'.pdf"
        drop *_bot *_top
    }
}
else di as error "*** 05_nfa: G2/G2b PULADOS (colunas de fluxo ausentes em P ou A) ***"

drop _pos _neg
di as text "05_nfa[${PFX}] done -> `out'/"
