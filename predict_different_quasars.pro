pro predict_different_quasars

; Calculate the signal-to-noise ratio (S/N) of potentially predictable, yet acausal, colour-difference of two near-antipodal quasars from separate ground-based telescopes.  This provides the potential for a truly loophole-free Bell-test: for such quasars, a quantum-mechanical (QM) experiment could exploit the independent fluxes of these sources as a switching mechanism. Observations of their simultaneously-sampled colours should be unspoiled by interference between the observers at the two sites, as those can be separated by a light-travel-time larger than the necessary photon-sampling rate of the quasars to reach the necessary S/N. So the test is to experimentally verify no correlation between the simultaneously-measured colours of the two sources.  The code simply runs the individual programs, with their settings as set in their respective parameter input tables.  There is a test dataset provided from Gemini `Alopeke/Zorro fast-framerate imaging of two stars, which is analyzed, and although these objects are not acuasal they still provide real-world calibration of the estimates.  Output test files are saved, and figures stored as JPEGs; these are as reported in the associated journal paper. The output is a list of potential quasar pairs (which would be acausal) observable with Gemini, and the best nights in the coming year to achieve the required S/N to verify they are uncorrelated.
; Eric Steinbring, 7 November 2024

; first, calculate the baseline, required signal-to-noise ratio
calculate_signal
; now, perform photometry on the test dataset from Gemini Alopeke/Zorro
do_fast_photometry
; and compare with a Bell-inequality experimental result performed at those two sites, and limited by noise; hence calibrating what is possible
simulate_quasar_correlation
;; format a list of potential quasar sources within which to search, which must be run prior to the search; these tables are provided pre-generated in the public code
;cull_quasars
; and with that, search for quasar pairs that fall outside the horizon that would allow these sites to spoil the test
find_quasar_pairs

end
