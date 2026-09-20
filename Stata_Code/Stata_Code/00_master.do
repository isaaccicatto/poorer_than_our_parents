* ============================================================================
* 00_MASTER — orchestrates: positive -> normative -> plots, per scenario,
*             then builds the sensitivity table from all accumulated results.
* Edit ONLY: ROOT, the globals block, and the SCEN list.
* ============================================================================
clear all
set more off
set graphics on
graph set window fontface "Times New Roman" 

global ROOT "/Users/isaaccicatto/Library/Mobile Documents/com~apple~CloudDocs/Faculdade/Bocconi - ESS/Thesis/Code/Code"
cd "$ROOT"

cap mkdir "Stata_Cache"
cap mkdir "Stata_Results"
cap mkdir "Plots"

global YEAR0   = 1970                       // t = 0
global J_END   = 2070                       // last cohort in the transfer pool
global HORIZON = 80                         // life-cycle window (ages 0..HORIZON)
global COHORTS "1973 1990 2006 2022"        // plotted generations

* ---- scenario list:  tag | region | csv path (relative to ROOT) ------------
* tag = short id used in every cached file / plot folder (no spaces).
* Sensitivity = add one line per (region, beta, sigma) run.
* ---- scenario list: one line per run ---------------------------------------
* format: "tag region csvpath"  (tag without spaces; path relative to ROOT)
* sensitivity = add one line per scenario, same format
local SCEN
local SCEN `"`SCEN' "base_RoW MM Extended_Path_final/output_MM/transition_results_MM.csv""'
*local SCEN `"`SCEN' "base_P A  Extended_Path_final/output_A/transition_results_A.csv""'
*local SCEN `"`SCEN' "base_A B  Extended_Path_final/output_B/transition_results_B.csv""'

* ---- Run reform scenarios ---- *
*local SCEN `"`SCEN' "repl_RoW MM Extended_Path_final/output_MM_repl/transition_results_MM_repl.csv""'
*local SCEN `"`SCEN' "repl_P A Extended_Path_final/output_A_repl/transition_results_A_repl.csv""'
*local SCEN `"`SCEN' "repl_A B Extended_Path_final/output_B_repl/transition_results_B_repl.csv""'

*local SCEN `"`SCEN' "by_RoW MM Extended_Path_final/output_MM_by/transition_results_MM_by.csv""'
*local SCEN `"`SCEN' "by_P A Extended_Path_final/output_A_by/transition_results_A_by.csv""'
*local SCEN `"`SCEN' "by_A B Extended_Path_final/output_B_by/transition_results_B_by.csv""'


/*
* ---- Run sensitivity analysis ---- *
local SCEN `"`SCEN' "beta_high_RoW MM  Extended_Path_final/output_MM_beta_high/transition_results_MM_beta_high.csv""'
local SCEN `"`SCEN' "beta_high_P A  Extended_Path_final/output_A_beta_high/transition_results_A_beta_high.csv""'
local SCEN `"`SCEN' "beta_high_A B  Extended_Path_final/output_B_beta_high/transition_results_B_beta_high.csv""'

local SCEN `"`SCEN' "beta_low_RoW MM  Extended_Path_final/output_MM_beta_low/transition_results_MM_beta_low.csv""'
local SCEN `"`SCEN' "beta_low_P A  Extended_Path_final/output_A_beta_low/transition_results_A_beta_low.csv""'
local SCEN `"`SCEN' "beta_low_A B  Extended_Path_final/output_B_beta_low/transition_results_B_beta_low.csv""'

local SCEN `"`SCEN' "sigma_high_RoW MM  Extended_Path_final/output_MM_sigma_high/transition_results_MM_sigma_high.csv""'
local SCEN `"`SCEN' "sigma_high_P A  Extended_Path_final/output_A_sigma_high/transition_results_A_sigma_high.csv""'
local SCEN `"`SCEN' "sigma_high_A B  Extended_Path_final/output_B_sigma_high/transition_results_B_sigma_high.csv""'

local SCEN `"`SCEN' "sigma_low_RoW MM  Extended_Path_final/output_MM_sigma_low/transition_results_MM_sigma_low.csv""'
local SCEN `"`SCEN' "sigma_low_P A  Extended_Path_final/output_A_sigma_low/transition_results_A_sigma_low.csv""'
local SCEN `"`SCEN' "sigma_low_A B  Extended_Path_final/output_B_sigma_low/transition_results_B_sigma_low.csv""'
*/

foreach s of local SCEN {
    tokenize `s'
    global TAG    "`1'"
    global REGION "`2'"
    global CSV    "`3'"

    di as text _n "{hline 60}"
    di as text "  SCENARIO $TAG   region=$REGION"
    di as text "  csv: $CSV"
    di as text "{hline 60}"

    do "Stata_Code/01_positive.do"
    do "Stata_Code/02_normative.do"
	do "Stata_Code/03_plots.do"
}

global PFX "base"
*do "Stata_Code/05_nfa.do" 
*do "Stata_Code/06_macro.do"
*do "Stata_Code/07_decomp.do"
*do "Stata_Code/04_table.do"

*do "Stata_Code/08_reform_diff.do"

