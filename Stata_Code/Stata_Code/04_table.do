* ============================================================================
* 04.1_TABLE — appends every Stata_Cache/res_*.dta into the sensitivity table.
* Re-running one scenario just overwrites its own res_ file; this stays true.
* ============================================================================
clear
local files : dir "Stata_Cache" files "res_*.dta"
local n : word count `files'
if `n' == 0 {
    di as error "no res_*.dta found in Stata_Cache/ — run scenarios first"
    exit 601
}
local first = 1
foreach f of local files {
    if `first' use "Stata_Cache/`f'", clear
    else append using "Stata_Cache/`f'"
    local first = 0
}

order tag region beta sigma ///
      rstar_ini rstar_fin ky_ini ky_fin nfa_ini nfa_fin
order C_star C_boomer C_BGP break_year ///
      S_pct_boomer S_pct_BGP eff_boomer eff_BGP, last
sort region beta sigma

replace rstar_ext = (rstar_ext - 1) * 100
replace rstar_ini = (rstar_ini - 1) * 100
replace rstar_fin = (rstar_fin - 1) * 100

save "Stata_Results/results_all.dta", replace
export delimited "Stata_Results/results_all.csv", replace

di as text _n "======= SENSITIVITY SUMMARY ======="
format rstar_ini rstar_fin %9.4f
format S_pct_boomer S_pct_BGP eff_boomer eff_BGP break_year %9.2f
list tag region beta sigma rstar_ini rstar_fin S_pct_BGP eff_BGP C_star break_year, ///
    sepby(region) abbreviate(12) noobs

di as text _n "full table -> Stata_Results/results_all.csv"



* ============================================================================
* 04.2_TEXTABLES — Stata_Results/results_all.dta -> Stata_Results/table_beta.tex
*                                           Stata_Results/table_sigma.tex
* UMA tabela por parametro: grade fixa de 3 subcolunas por cenario;
* sub-cabecalho muda por secao (Macro: t=0/t=100 | Cohorts: Life/Peak/Entry |
* Normative: celula mesclada).
* Overleaf: \usepackage{booktabs}  +  \input{table_beta}
* Math via \( \) — NUNCA "$" em string no Stata (expande global).
* ============================================================================
use "Stata_Results/results_all.dta", clear

capture program drop _vf
program define _vf, rclass
    args var tag fmt
    quietly summ `var' if tag == "`tag'", meanonly
    if r(N) == 0 | missing(r(mean)) return local s "--"
    else return local s = strtrim(string(r(mean), "`fmt'"))
end

local ECON "A RoW P"

foreach dim in beta sigma {

    if "`dim'" == "beta" {
        local scens "beta_low base beta_high"
        local sym "\beta"
    }
    else {
        local scens "sigma_low base sigma_high"
        local sym "\sigma"
    }

    * ---- header: parameter values, full precision, leading zero ------------
    local heads ""
    foreach pf of local scens {
        quietly summ `dim' if tag == "`pf'_A", meanonly
        local pv = strtrim(string(round(r(mean), 0.0001)))
        if substr("`pv'",1,1) == "." local pv = "0`pv'"
        local heads "`heads' & \multicolumn{3}{c}{\(`sym' = `pv'\)}"
    }

    tempname fh
    file open `fh' using "Stata_Results/table_`dim'.tex", write replace

    file write `fh' "\begin{table}[htbp]" _n
    file write `fh' "\centering" _n
    file write `fh' "\caption{Sensitivity to \(`sym'\)}" _n
    file write `fh' "\label{tab:sens_`dim'}" _n
    file write `fh' "\footnotesize" _n
    file write `fh' "\begin{tabular}{l ccc ccc ccc}" _n
    file write `fh' "\toprule" _n
    file write `fh' "`heads' \\" _n
    file write `fh' "\cmidrule(lr){2-4}\cmidrule(lr){5-7}\cmidrule(lr){8-10}" _n

    * ============== MACRO (t=0 | ext. (t) | t=infty) ==================
    file write `fh' "\textit{Macro} & \(t=0\) & ext. (t) & \(t=\infty\) & \(t=0\) & ext. (t) & \(t=\infty\) & \(t=0\) & ext. (t) & \(t=\infty\) \\" _n
    file write `fh' "\midrule" _n

    local line "\(r^{\star} (\%)\) (RoW)"
    foreach pf of local scens {
        _vf rstar_ini `pf'_A "%9.2f"
        local line "`line' & `r(s)'"
        _vf rstar_ext `pf'_A "%9.2f"
        local v1 "`r(s)'"
        _vf rstar_extt `pf'_A "%9.0f"
        if "`v1'" == "--" local line "`line' & --"
        else              local line "`line' & `v1' (`r(s)')"
        _vf rstar_fin `pf'_A "%9.2f"
        local line "`line' & `r(s)'"
    }
    file write `fh' "`line' \\" _n

    foreach e of local ECON {
        local line "\((K/Y)_{`e'}\)"
        foreach pf of local scens {
            _vf ky_ini `pf'_`e' "%9.2f"
            local line "`line' & `r(s)'"
            _vf ky_ext `pf'_`e' "%9.2f"
            local v1 "`r(s)'"
            _vf ky_extt `pf'_`e' "%9.0f"
            if "`v1'" == "--" local line "`line' & --"
            else              local line "`line' & `v1' (`r(s)')"
            _vf ky_fin `pf'_`e' "%9.2f"
            local line "`line' & `r(s)'"
        }
        file write `fh' "`line' \\" _n
    }
    foreach e in  A P {
        local line "\(F_{`e'}\)"
        foreach pf of local scens {
            _vf nfa_ini `pf'_`e' "%9.2f"
            local line "`line' & `r(s)'"
            _vf nfa_ext `pf'_`e' "%9.2f"
            local v1 "`r(s)'"
            _vf nfa_extt `pf'_`e' "%9.0f"
            if "`v1'" == "--" local line "`line' & --"
            else              local line "`line' & `v1' (`r(s)')"
            _vf nfa_fin `pf'_`e' "%9.2f"
            local line "`line' & `r(s)'"
        }
        file write `fh' "`line' \\" _n
    }

    * ============== POSITIVE (Lifetime | Peak | Entry) =======================
    file write `fh' "\addlinespace" _n
    file write `fh' "\textit{Positive} & Life & Peak & Entry & Life & Peak & Entry & Life & Peak & Entry \\" _n
    file write `fh' "\midrule" _n

    foreach g in 1973 1990 2006 2022 {
        foreach e of local ECON {
            local line "\(C_{`e',`g'}\)"
            foreach pf of local scens {
                _vf npv_`g' `pf'_`e' "%9.2f"
                local line "`line' & `r(s)'"
                _vf pk_`g' `pf'_`e' "%9.2f"
                local v1 "`r(s)'"
                _vf pka_`g' `pf'_`e' "%9.0f"
                local line "`line' & `v1' (`r(s)')"
                _vf mn_`g' `pf'_`e' "%9.2f"
                local line "`line' & `r(s)'"
            }
            file write `fh' "`line' \\" _n
        }
    }

    * ================== NORMATIVE (merged 3 cols) ===========================
    file write `fh' "\addlinespace" _n
    file write `fh' "\textit{Normative} & \multicolumn{3}{c}{} & \multicolumn{3}{c}{} & \multicolumn{3}{c}{} \\" _n
    file write `fh' "\midrule" _n

    foreach e of local ECON {
        local line "\(\hat{C}^{\mathrm{Rawls}}_{`e'}\)"
        foreach pf of local scens {
            _vf C_star `pf'_`e' "%9.2f"
            local line "`line' & \multicolumn{3}{c}{`r(s)'}"
        }
        file write `fh' "`line' \\" _n
    }
    foreach e of local ECON {
        local line "\(\text{break-even}_{j,`e'}\)"
        foreach pf of local scens {
            _vf break_year `pf'_`e' "%9.0f"
            local line "`line' & \multicolumn{3}{c}{`r(s)'}"
        }
        file write `fh' "`line' \\" _n
    }
    foreach e of local ECON {
        local line "\(\Delta_{KH,`e'}^{\mathrm{SS0=BGP}}\) (\%)"
        foreach pf of local scens {
            _vf S_pct_BGP `pf'_`e' "%9.1f"
            local line "`line' & \multicolumn{3}{c}{`r(s)'}"
        }
        file write `fh' "`line' \\" _n
    }
    foreach e of local ECON {
        local line "\(\Delta_{KH,`e'}^{\mathrm{SS0=1973}}\) (\%)"
        foreach pf of local scens {
            _vf S_pct_boomer `pf'_`e' "%9.1f"
            local line "`line' & \multicolumn{3}{c}{`r(s)'}"
        }
        file write `fh' "`line' \\" _n
    }
    file write `fh' "\bottomrule" _n
    file write `fh' "\multicolumn{10}{l}{\footnotesize ext.: maximum or minimum; -- if monotonic. Model period: 1970-2070.} \\" _n
	    file write `fh' "\multicolumn{10}{l}{\footnotesize For Normative, consider 1970-2070, the demographic transition duration.} \\" _n
		    file write `fh' "\multicolumn{10}{l}{\footnotesize For Positive, consider expected lifetime of 80 periods after entry in the model.} \\" _n
    file write `fh' "\end{tabular}" _n
    file write `fh' "\end{table}" _n

    di as text "wrote Stata_Results/table_`dim'.tex"
}
