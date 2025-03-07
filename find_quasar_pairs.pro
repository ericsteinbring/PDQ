; load some useful functions first
forward_function quantum_correlation
forward_function bell_theorem
forward_function airmass_simple
forward_function sf_str
forward_function sf_trans_dec
forward_function sigfig
forward_function gcirc

pro find_quasar_pairs

; Finds available quasar pairs testable in an Earth-wide Bell-theorem test; those are taken from the output catalog of cull_quasars. It has the potential to restrict to an angular offset from a global, common timing/photometric calibrator, which can be shut off (by setting this to 180 degrees).  It also queries Vizier to report if a star is within the field, in order to provide simultaneous photometric calibration individually for each quasar. Note that colour and redshift-limit settings only select a pre-determined catalog by file labels output by cull_quasars, so those catalogues must be pre-generated to run for each setting, and must match. Note also that colour here is the absolute difference between the two quasars from the catalogue, not the error on each - although, in principle, that is what a test measures, as a predicted colour-difference (within Earth-crossing time) constitutes a signal. Two known bugs are that there can be a parsing fault if the input quasar catalog has any unexpected symbol (such as a comma) and so failing to find pairs (even giving a negative number of those!) and crash if the list is too long. The second bug is avoided by restricting the input magnitude limits.
; Eric Steinbring, 29 January 2019.

; Updated: 21 January 2024.  The code now allows testing against different potential probabilities of correlation with target separation, either as in the original paper, or an exponentially-decaying one. In either case, these make a prediction for which sources to pick which - although presumably providing uncorrelated photons - might be able to spoil a Bell-Theorem test undertaken with the Gemini telescopes as the quasar-monitoring telescopes.  These would be found providing photons to each telescope (sampling at a rate faster than their signaling speed) at least 50% correlated by colour (or more); which is just testable with Alopeke/Zorro, at an exposure time under 0.03 seconds. The pair-correlation is also slightly refined, to instead give the hour within the best night of observations.  It is assumed that an observation block of something under an hour is still useful, even if both targets are near the horizon, and the target as viewed from Hawai'i (Gemini-North) is setting - in effect, some twilight is usable.  Note that, at the beginning of what is here called a "night" it would be just sunset at Gemini-North, although at Gemini-South (Chile) the sun would have already set about 7 hours before; so a night (with both sites in civil darkness) lasts under 5 hours at any time of year, and that of astronomical darkness is less.

; Version: 28 September 2024.  Compares with correlation of real photometry as well, in order to check with the model. So far, this is just test data, which are overlapping Alopeke/Zorro observations available from the archive.  Those have exposure times of 0.06 seconds, so cannot be used to look for a signal (of known-acausal, unspoiled switching photons) as defined here.  They are, however, sufficient to report what S/N is achievable with bright sources, and so can be plotted along with the simulation as a "calibration" set.

; parameters and controls

; sites
latitude_north    = 19.8238 ; 28.76 ; deg ; either Maunakea or La Palma
longitude_north   = 155.469 ; 17.88, deg, W ; either Maunakea or La Palma
altitude_north    = 4213. ; 2396. ; meters, a.s.l. ; either Maunakea or La Palma
latitude_south    = -30.2407 ; deg, N ; Cerro Pachon
longitude_south   = 70.7367 ; -70.7367 ; 70.7367 ; deg, W ; Cerro Pachon
altitude_south    = 2722. ; meters, a.s.l. ; Cerro Pachon
hours             = 5 ; 4 ; 5 ; 7 ; 12 ; 24 ; hours of a "night", the longest possible continuous duration for each pair, i.e. visible at both sites, and whether astronmical darkness is required
night_offset      = (4./24.) ; (6 PM local sunset on Maunkea, 4 AM UTC) ; -(3./24.) ; (6 PM local sunset on Pachon, 9 PM UTC); offset for combined-telescope nights (set to Hawai`i if both sites to be dark)

; model telescope and sky conditions, settings sampling and signal-to-noise requirement
alpha_sky         = 0.250 ; mag, extinction, at zenith
beta_sky          = 0.5 ; 0.2 ; 0.4 ; 0.5 ; 1. ; mag per square-arcsec, sky brightness beyond the darkest possible
gamma_sky         = 0.8 ; 0.375 ; arcsec, best seeing, at zenith
epsilon           = 0.1 ; 0.01 ; 0.02 ; 0.05 ; 0.1 ; 0.250 ; mag, photometric error
zeta              = 0.01 ; 0.005 ; 0.01 ; 0.020 ; mag, photometric zeropoint
omega             = 0.02 ; 1. ; 0.5 ; 0.2 ; 0.1 ; 0.05 ; 0.04 ; 0.025 ; 0.02 ; 0.01 ; 0.005 ; mag, difference sensed, either width of quasar brightness distribution or the colour; ideally, limited by zeropoint
correlation       = 'exp' ; 'bell' ; 'exp' ; whether to test that based on distribution derived from Bell's-theorem as in earlier paper, or on exponential decay
fold_over         = 'no' ; 'no' ; 'yes' ; whether to "fold" the observed distribution, looking for absolute differences, as was reported in the original paper
e_fold            = 2. ; 0. ; 0.5 ; 1. ; 1.5 ; 2. ; 2.71828173. ; 3. ; 5. ; 10. ; e-folding angular rate at which correlation must drop off from perfect
samples           = 3000 ; 1.e6 ; 1.e5 ; 50000 ; 10000 ; 5000 ; 3000 ; 1000 ; 100 ; 10 ; number of simultaneous samples taken for each quasar in the observation, default is 1000
ratio_wanted      = 10. ; 3. ; 5. ; 10. ; 100. ; 1000. ; the level of discrimination demanded for a successful detection with accausal pairs 

; parameters; for "generic" case all angles are allowed, but default values are reset by choosing "peak", "truly acausal" or "same-field" setting
calibration_star  = 16. ; 5. ; 7.5 ; 10. ; 12. ; 15. ; 15.5 ; 16. ; 16.7 ; 17. ; 18. ; 19. ; 20. ; mag; the brightness of the global calibrator, which sets the S/N limit of calibration
calibration_limit = 180. ; 45. ; 60. ; 90. ; 120. ; 180. ; degrees, max. sep. on sky from the target to a global calibrator, where default is 180 deg and so shut off
separation_limit  = 0. ; 60. ; 120. ; 170. ; 175. ; degrees, minimum angular sep. between targets, default is zero for the generic whole-sky case
field_limit       = 180. ; 0.1 ; 5. ; 180. ; degrees, in same-field case this is the FoV, otherwise maximum angular sep. range from the minimum, default 180 deg
airmass_limit     = 2. ; 1.5 ; 2. ; 2.5 ; 3. ; 5. ; airmass within which to allow either target, default is 2
magnitude_limit   = 15.5 ; 15.5 ; 20. ; 19.5 ; 19. ; 18.5 ; 18. ; 17.5 ; 17. ; 16.5 ; 16. ; 15.5 ; 15. ; mag., max. Red mag, default is 15.5 mag
colour_limit      = 0.1 ; 0. ; 0.1 ; 0.2 ; 0.25 ; 0.5 ; 1., 2., magnitudes, catalog-colour upper absolute difference, default 0.2
redshift_limit    = 1.48 ; 0. ; 0.5, 1., 1.48 ; 1.57 ; 2., 3., 3.65 ; 4., 5., 6. ; min. redshift limit: none, at Hubble sphere, or truly acausal
separation_star   = 60. ; 30. ; 35. ; 40. ; 60. ; 300. ; arcseconds, default of 35 centrally fits Alopeke/Zorro unvignetted widefield setting, 300 is for GMOS-N/GMOS-S
minimum_nights    = 10. ; 1. ; 7. ; 10. ; 14., 15. ; 21. ; 28. ; 30. ; min. number of nights targets are both visible per year
calendar_year     = 2025 ; 2019 ; 2020 ; 2021 ; 2022 ; 2023 ; 2024 ; 2025
semester          = 'AB' ; 'A' ; 'B' ; 'AB' ; A is northern spring/summer semester, B is fall/winter, AB is both

; experiment and plotting toggles
experiment        = 'all' ; 'generic' ; 'peak' ; 'truly_acausal' ; 'same_field' ; 'test' ; 'all' ; default is all, or limited to the generic and other cases, including the test observations
calibrate         = 'no' ; 'yes' ; 'no' ; whether to check for a similar-brightnesss star within each field; 'no' returns the response 'self-calibrated' instead.
catalog           = 'gaia' ; 'gaia' ; 'ucac' ; choice of star catalog to use in finding calibration stars
just_model        = 'no' ; 'yes' ; 'no' ; whether to just plot the model, and stop
print_long_info   = 'yes' ; 'no' ; 'yes' ; whether to print all the available information, or not
plotting          = 'all' ; 'ratio' ; 'sky' ; 'polar' ; 'all' ; whether to plot the S/N of the targets, positions on the sky, or their orientations as viewed from the two sites; default is to just plot the ratio

; settings
; set the generic case to be the initial values
separation_limit_generic = separation_limit
field_limit_generic = field_limit
magnitude_limit_generic = magnitude_limit
colour_limit_generic = colour_limit
redshift_limit_generic = redshift_limit
separation_star_generic = separation_star
; effect at peak, that is, either that of both targets being near zenith for the two telescopes, or at a potential peak in the causal photons visible
separation_limit_peak = 125. ; 125. ; 120. ; 110. ; 100. ; 95. ; 90. ; 80. ; 65. ; 60.
field_limit_peak = 5. ; 15. ; 10. ; 5. ; 1.
magnitude_limit_peak = 15.5 ; 17.5 ; 17. ; 16.5 ; 16. ; 15.5 ; 15. 
colour_limit_peak = 0.1
redshift_limit_peak = 1.48 ; 1.57
separation_star_peak = 60. ; 35.
; truly acausal
separation_limit_truly_acausal = 175. ; 165. ; 170. ; 175.
field_limit_truly_acausal = 5. ; 15. ; 10. ; 5.
magnitude_limit_truly_acausal = 18.5 ; 20. ; 19.5 ; 19. ; 18.5 ; 18. ; 17.5; default is 18.5 mag
colour_limit_truly_acausal = 0.5 ; 0.1 ; 0.5 ; 2. 
redshift_limit_truly_acausal = 4. ; 3.65 ; 4. ; 4.5 ; 5. ; 5.5 ; 6.
separation_star_truly_acausal = 60. ; 35. 
; same field
separation_limit_same_field = 0.
field_limit_same_field = 0.083 ; 0.017 ; 0.083 ; degrees, the FoV of Alopeke/Zorro (0.017 degrees in wide-field mode) or GMOS-N/GMOS-S (0.083 X 0.083 degrees)
magnitude_limit_same_field = 15.5 ; 18.5 ; 18. ; 17.5 ; 17. ; 16.5 ; 16. ; 15.5 ; 15.
colour_limit_same_field = 1. ; 1. ; 0.5 ; 0.1
redshift_limit_same_field = 0. ; 1.48 ; 0.
separation_star_same_field  = 60.*60. ; 30. ; 300./2. ; 60.*60. ; arcsec, the FoV of Alopeke/Zorro (about a 1' circle in wide-field mode) or GMOS-N/GMOS-S (square, 5.5 arcmin on a side) or a degree

; position of global calibration source
; Crab pulsar, R.A.=83.633212 deg, Dec.=22.014460 deg (J2000.0), R=16.17 mag
ra_c = 83.633212
dec_c = 22.014460

; find the maximum zenith angle at the airmass limit
zenith_limits = 90.*findgen(100)/float(100)
airmass_limits = zenith_limits
for i = 0, n_elements(zenith_limits) - 1 do begin
 airmass_limits(i) = airmass_simple(zenith_limits(i))
; print, zenith_limits(i), airmass_limits(i)
endfor
zenith_limit = max(zenith_limits(where(airmass_limits le airmass_limit)))
print, 'Zenith angle limit at airmass limit (degrees): ', zenith_limit

; plotting and experiment settings
if plotting eq 'polar' then plottings = ['polar']
if plotting eq 'sky' then plottings = ['sky']
if plotting eq 'ratio' then plottings = ['ratio']
if plotting eq 'all' then plottings = ['ratio', 'sky', 'polar']
if experiment eq 'test' then experiments = ['test']
if experiment eq 'generic' then experiments = ['generic']
if experiment eq 'peak' then experiments = ['peak']
if experiment eq 'truly_acausal' then experiments = ['truly_acausal']
if experiment eq 'same_field' then experiments = ['same_field']
if experiment eq 'all' then experiments = ['generic', 'peak', 'truly_acausal', 'same_field', 'test']
separation_star_default = 300. ; 300./2. ; arcseconds, a default setting unless reset above is a box the size of the GMOS field

; setup display
window, 10, xsize = 500, ysize = 200, title='Colours'
wset, 10
!p.multi=[0,1,1]
charsize=2.5
loadct, 0, /silent ; greyscale, for models
;loadct, 1, /silent ; bluescale, for models
;loadct, 39, /silent ; colour, for data
!p.background=255
device, set_font='helvetica bold', /tt_font
!p.font=1
!except=0

; plot a colour bar
loadct, 39, /silent ; colour
color_line = findgen(256)
plot, color_line, 0.*color_line, xstyle=1, xrange=[0., 256.], ystyle=1, yrange=[0., 1.], charsize=0.5*charsize, color=0 ;, color=255
for i = 0, 255 do begin
 oplot, [i, i], [0., 1.], thick=2, color=i
endfor

; calculations
calibration = exp(magnitude_limit_truly_acausal - calibration_star) ; how much brighter the calibration source is; i.e. calibrator compared to the accausal antipodal quasars of the test
;calibration = exp(magnitude_limit - calibration_star) ; how much brighter the calibration source is; i.e. calibrator compared to the generic-case quasars of the test
resolution = 1000 ; floor(sqrt(2e6)) ; 1000 ; 500 ; 100 ; resolution of calculations
phi = 180.*findgen(resolution)/float(resolution) ; degrees
print, 'Calibration factor: ', calibration

; probability distribution; as considered in original paper, defined via the Bell-inequality, and demanding a whole-sky average 50% probability
M = sqrt(3./4.) ; 1. ; 0.87 ; 0.866 ; sqrt(3./4.) ; 0.75
;R = 0.67 ; 1. ; 0.87 ; 0.68 ; 0.67 ; 0.667 ; 2./3. ; 0.66 ; 0.5 ; 0.
p = double(bell_theorem(phi))
q = double(quantum_correlation(phi))
R = 2.*(1./!pi)*(1./float(resolution))*total(M*sqrt(1. + (1./M^2.)*p^2. - (1./M^2.)*q^2.))
C = M^2.*(1./float(resolution))*total(sqrt(1. + (1./M^2.)*p^2. - (1./M^2.)*q^2.))
M_0 = C - M ; sqrt(3./4.)
;print, 'M, R, C, M_0;', M, R, C, M_0
D = sqrt(M^2. + p^2. - q^2.)
D_R = (3./2.)*(D - R)
D_0 = 0.5*(abs(D_R) - R)
D_p = 0.5*(D_R - R)
D_0_p = abs(median(D_0))/4.
P_bell = 0.5 + 0.5*D_0
if fold_over eq 'yes' then P_bell(where(abs(D_R) ge 2.*R)) = 0.5 + 2.*D_0(where(abs(D_R) ge 2.*R)) - (3./4.)*R
P_0 = 0.5 + 0.5*D_0
P_0_p = 0.5 + 0.5*D_p
if correlation eq 'bell' then begin
print, 'Based on Bell-theorem'
P = P_bell
;f_linear = phi/180. ; linear with angle, that is just unity at antipodes
;P = f_linear*P_bell
;P = (1.+max(P)-min(P)/mean(P))*P ; possibly renormalizing
endif

; or instead an exponentially decaying probability distribution, which can accommodate perfect correlation only at zero target offset; here labelled as the "prime" case
e = exp(1.)
e_fold_prime = 5. ; 2. ; 3. ; 5. ; 10. ; 20. ; 100.
;print, e, alog(e), exp(0.)
; first, assume some fraction of truly acausal target photons by separation
f_linear = phi/180. ; linear with angle, that is, zero when coincident (or, in fact, the same object) and exactly unity at the antipodes
f_flat = replicate(0.5, n_elements(phi)) ; replicate(1., n_elements(phi)) ; either half or unity at all angles, the latter impossible as it implies infinite redshift
f_exp = 0. + (1./e)*exp(e^e_fold*phi/180.) ; exponential, that is, strongly restricted by redshift
;f_exp = 0. + (1./e)*exp(e^e_fold_prime*phi/180.) ; exponential at the limit, i.e. prime case
f_exp = f_exp/max(f_exp)
; tests of swapping the conditions
;f_acausal = f_exp
f_acausal = f_linear
;f_acausal = f_flat
; the probability of target photons being correlated with angle, despite being acausal
;P_linear = 1.-phi/180. ; instead a simple-minded linear case, where probability is zero at antipodes
P_linear = 1.-0.5*phi/180. ; instead a simple-minded linear case, which can also be exactly 50% probability at antipodes
P_flat = replicate(1., n_elements(phi)) ; and the opposite extreme, that is, where flat and so 100% correlation across the whole sky
P_fair = replicate(0.5, n_elements(phi)) ; instead, a limit of exactly 50:50 at all angles
P_exp = 0.5+0.5*(1./e)*exp(1.-e^e_fold*phi/180.) ; exponentially decaying from zero separation
P_prime = 0.5+0.5*(1./e)*exp(1.-e^e_fold_prime*phi/180.) ; the strong case, such as, the two objects are in fact coincident
; combining together the fractions
;P_comb = f_linear*P_fair*P_exp ; exponential case
;P_comb = f_linear*P_fair*P_prime ; instead, using the primed case
P_comb = f_acausal*P_exp ; exponential case
if correlation eq 'exp' then begin
 print, 'Based on exponential decay'
; P = P_exp
 P = P_comb
; P = (1.+max(P)-min(P)/mean(P))*P ; possibly renormalizing
endif

; report probability statistics
print, 'Max, min, mean, median probability across the sky: ', max(P), min(P), mean(P), median(P) 

; set up the display
;window, 11, xsize = 1000, ysize = 1000, title='Fractions and probability of correlated photons'
window, 11, xsize = 1000, ysize = 750, title='Probability of correlated photons'
wset, 11
;!p.multi = [0, 1, 2]
!p.multi = [0, 1, 1]
charsize=2.5
loadct, 0, /silent ; greyscale

;; plot fraction of acausal photons
;plot, phi, f_flat, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)', ytitle='Fraction of acausal-pair photons', ystyle=1, yrange=[-0.05, 1.05], ymargin=[0,2], charsize=charsize, color=0
;; shading
;loadct, 0, /silent ; greyscale
;for i = 0, n_elements(phi)-1 do begin
; oplot, [phi(i), phi(i)], [-0.05, f_flat(i)], thick=2, color=225 ; light grey
;endfor
;; replot
;oplot, phi, f_flat, color=0
;oplot, phi, f_acausal, color=0
;; angle limits, with coloured delinating lines for those regions
;loadct, 39, /silent ; colour, angle limit
;oplot, [separation_limit_truly_acausal, separation_limit_truly_acausal], [-0.05, 1.05], thick=5, color=220 ; red
;arrow, separation_limit_truly_acausal, 0.5, separation_limit_truly_acausal+5., 0.5, thick=2, color=220, /solid, /data
;oplot, [separation_limit_truly_acausal-1., separation_limit_truly_acausal-1.], [-0.05, 1.05], thick=5, color=150 ; green
;arrow, separation_limit_truly_acausal-1., 0.5, separation_limit_truly_acausal-1.-5., 0.5, thick=2, color=150, /solid, /data
;oplot, [90., 90.], [-0.05, 1.05], thick=5, color=50 ; blue
;arrow, 90., 0.5, 85., 0.5, thick=2, color=50, /solid, /data
;; limits
;loadct, 0, /silent ; greyscale
;oplot, [182., 182.], [-0.05, 1.05], thick=25, color=75 ; dark grey, beyond horizon
;; labels
;;xyouts, 95., 0.85, 'Hubble radius', charsize=charsize, charthick=charthick, color=0
;;xyouts, 91., 0.85, '|', orientation=-90., charsize=charsize, charthick=charthick, color=0
;xyouts, 135., 0.55, 'Acausal limit', charsize=charsize, charthick=charthick, color=0
;xyouts, 175., 0.55, '|', orientation=90., charsize=charsize, charthick=charthick, color=0
;; clean up
;;axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
;;axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
;axis, yaxis=0, color=0, ystyle=1, yrange=[-0.05, 1.05], charsize=charsize
;axis, yaxis=1, color=0, ystyle=1, yrange=[-0.05, 1.05], ytickformat='(A1)'

; plot probability distribution
;plot, phi, P, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Probability of photon correlation', ystyle=1, yrange=[-0.05, 1.05], ymargin=[3.5, 0], charsize=charsize, thick=3, color=0
plot, phi, P, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Probability of photon correlation', ystyle=1, yrange=[-0.05, 1.05], charsize=charsize, thick=3, color=0
; shading
loadct, 0, /silent ; greyscale
for i = 0, n_elements(phi)-1 do begin
 oplot, [phi(i), phi(i)], [0., P_prime(i)], thick=2, color=225 ; light grey
endfor
; replot
oplot, phi, P_prime, color=0
oplot, phi, P_flat, color=0
oplot, phi, P_fair, linestyle=3, color=0
;oplot, phi, P_linear, color=0
oplot, phi, P_exp, thick=3, color=0
oplot, phi, P, thick=3, color=100
; angle limits, with coloured delinating lines for those regions
loadct, 39, /silent ; colour, angle limit
oplot, [separation_limit_truly_acausal, separation_limit_truly_acausal], [zeta/2., 1.05], thick=5, color=220 ; red
;arrow, separation_limit_truly_acausal, 0.75, separation_limit_truly_acausal+5., 0.75, thick=2, color=220, /solid, /data
oplot, [separation_limit_truly_acausal-1., separation_limit_truly_acausal-1.], [zeta/2., 1.05], thick=5, color=150 ; green
;arrow, separation_limit_truly_acausal-1., 0.75, separation_limit_truly_acausal-1.-5., 0.75, thick=2, color=150, /solid, /data
oplot, [90., 90.], [zeta/2., 1.05], thick=5, color=50 ; blue
;arrow, 90., 0.75, 85., 0.75, thick=2, color=50, /solid, /data
; limits
loadct, 0, /silent ; greyscale
oplot, [182., 182.], [0., 1.05], thick=25, color=75 ; dark grey, beyond horizon
oplot, [-5., 185.], [mean(P), mean(P)], linestyle=2, color=0
;oplot, [-5., 185.], [median(P), median(P)], linestyle=3, color=0
oplot, [-5., 185.], [0., 0.], linestyle=1, color=0
; labels
xyouts, 7.5, 0.4, 'Two telescopes and one source', charsize=charsize, charthick=charthick, color=0
xyouts, 7.5, 0.46, '|', charsize=charsize, charthick=charthick, color=0
xyouts, 120., 0.65, 'Due to:', charsize=charsize, charthick=charthick, color=0
xyouts, 120., 0.57, 'All photons', charsize=charsize, charthick=charthick, color=0
xyouts, 120., 0.51, '|', charsize=charsize, charthick=charthick, color=0
xyouts, 120., 0.27, 'Acausal photons', charsize=charsize, charthick=charthick, color=0
;xyouts, 120., 0.30, '|', charsize=charsize, charthick=charthick, color=0
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-0.05, 1.05], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[-0.05, 1.05], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
;write_jpeg, 'figure_fraction_and_probability.jpg', screen_output, quality=100, true=1
write_jpeg, 'figure_probability.jpg', screen_output, quality=100, true=1

; signal
;signal = omega*(P-mean(P))
;signal_prime = omega*(P_0_p-mean(P))
signal = omega*P
signal_flat = omega*P_flat
signal_fair = omega*P_fair
signal_prime = omega*P_prime
signal_causal = omega*(1.-f_acausal)*P_exp
signal_0 = 0.5 ; meaning that, otherwise samples are perfectly fair
max_signal = max([signal, signal_flat, signal_fair, signal_prime, signal_causal])
print, 'Omega (mag): ', omega
print, 'No signal: ', signal_0
print, 'Max. signal (mag): ', max_signal

; set up the display
window, 12, xsize = 1000, ysize = 1000, title='Ideal Signal and Observational Noise'
wset, 12
!p.multi = [0, 1, 2]
charsize=2.5
loadct, 0, /silent ; greyscale

; plot signal
plot, phi, signal, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)', ytitle='Min. correlated signal per exp. (mag)', charsize=charsize, ystyle=1, yrange=[0., 1.25*max_signal], ymargin=[0,2], thick=3, color=0
; shading
loadct, 0, /silent ; greyscale
for i = 0, n_elements(phi)-1 do begin
 oplot, [phi(i), phi(i)], [0., signal_prime(i)], thick=2, color=225 ; light grey
endfor
; replot
oplot, phi, signal_prime, color=0
oplot, phi, signal_flat, color=0
oplot, phi, signal_fair, linestyle=3, color=0
oplot, phi, signal_causal, thick=3, color=0
oplot, phi, signal, thick=3, color=100
; angle limits, with coloured delinating lines for those regions
loadct, 39, /silent ; colour, angle limit
oplot, [separation_limit_truly_acausal, separation_limit_truly_acausal], [0., 1.25*max_signal], thick=5, color=220 ; red
oplot, [separation_limit_truly_acausal-1., separation_limit_truly_acausal-1.], [0., 1.25*max_signal], thick=5, color=150 ; green
oplot, [90., 90.], [0., 1.25*max_signal], thick=5, color=50 ; blue
;arrow, 90., median(signal), 85., median(signal), thick=2, color=50, /solid, /data
; limits
loadct, 0, /silent ; greyscale
oplot, [182., 182.], [0., 1.25*max_signal], thick=25, color=75 ; dark grey, beyond horizon
oplot, [-5., 185.], [mean(signal), mean(signal)], linestyle=2, color=0
;oplot, [-5., 185.], [median(signal), median(signal)], linestyle=3, color=0
oplot, [-5., 185.], [0., 0.], linestyle=1, color=0
; labels
xyouts, 95., 0.89*min(signal_flat(where(phi gt 90.))), 'Perfect correlation of all sources', charsize=charsize, charthick=charthick, color=0
xyouts, 95., 0.96*min(signal_flat(where(phi gt 90.))), '|', charsize=charsize, charthick=charthick, color=0
xyouts, 95., 1.15*min(signal_prime(where(phi gt 90.))), 'Zeropoint limit for a single exp.', charsize=charsize, charthick=charthick, color=0
xyouts, 95., 1.02*min(signal_prime(where(phi gt 90.))), '|', charsize=charsize, charthick=charthick, color=0
; clean up
;axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
;axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[0., 1.25*max_signal], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[0., 1.25*max_signal], ytickformat='(A1)'

; noise
max_noise = 2. ; 1000. ; 10. ; mag
print, 'Samples: ', samples
print, 'Maximum noise (mag): ',  max_noise
error = airmass_simple(phi) ; one telescope at one site, error limit
;error_single = abs(2. - sqrt(2.)*airmass_simple(phi)) ; two telescopes at one site
error_single = (airmass_limit/2.)*abs(2. - sqrt(2.)*(2./airmass_limit)*airmass_simple(phi)) ; two telescopes at one site
;error_twin = abs(2. - sqrt(2.)*airmass_simple(phi/2.)) ; two telescopes at sites separated by 90 degrees, very close to Gemini's configuration
error_twin = (airmass_limit/2.)*abs(2. - sqrt(2.)*(2./airmass_limit)*airmass_simple(phi/2.)) ; two telescopes at sites separated by 90 degrees, very close to Gemini's configuration
error_antipodal = airmass_simple(180.-phi) ; two telescopes at antipodal sites
;error_antipodal = abs(airmass_limit - sqrt(2.)*airmass_simple(180.-phi)) ; two telescopes at antipodal sites
noise =  zeta*error/sqrt(float(samples)) ; penalized only by the zeropoint
noise_0 = zeta*error ; fixed at a single observation
;print, max(noise), max(noise_0)
noise(where(noise ge max(noise))) = max_noise
noise_0(where(noise ge max(noise_0))) = max_noise
;noise(where(phi gt 90.)) = max_noise
;noise_0(where(phi gt 90.)) = max_noise
noise(where(noise gt max_noise)) = max_noise
noise_0(where(noise_0 gt max_noise)) = max_noise
noise_single =(error_single*(alpha_sky + error_single^2.*beta_sky*gamma_sky^2.) + epsilon)/sqrt(float(samples)) + zeta ; for two telescopes at a single site, with a penalty for independent zeropoints
;noise_single(where(phi gt 90.)) = max_noise
noise_single(where(noise_single gt max_noise)) = max_noise
noise_single_0 = error_single*(alpha_sky + error_single^2.*beta_sky*gamma_sky^2.) + epsilon + zeta ; that for a single site, fixed at a single observation
;noise_single_0(where(phi gt 90.)) = max_noise
noise_single_0(where(noise_single_0 gt max_noise)) = max_noise
noise_twin = (error_twin*(alpha_sky + error_twin^2.*beta_sky*gamma_sky^2.) + epsilon)/sqrt(float(samples)) + zeta ; assuming a penalty for two independent zeropoints
noise_twin_0 = error_twin*(alpha_sky + error_twin^2.*beta_sky*gamma_sky^2.) + epsilon + zeta ; fixed at a single observation, assuming a penalty for two independent zeropoints
noise_twin_0(where(noise_0 gt max_noise)) = max_noise
noise_antipodal = (error_antipodal*(alpha_sky + error_antipodal^2.*beta_sky*gamma_sky^2.) + epsilon)/sqrt(float(samples)) + sqrt(2.)*zeta ; assuming a penalty for two independent zeropoints
noise_antipodal_0 = error_antipodal*(alpha_sky + error_antipodal^2.*beta_sky*gamma_sky^2.) + epsilon + sqrt(2.)*zeta ; fixed at a single observation, assuming a penalty for two independent zeropoints
noise_antipodal_0(where(noise_antipodal_0 gt max_noise)) = max_noise

; plot noise
plot, phi, noise, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Observational noise (mag)', charsize=charsize, ystyle=1, yrange=[zeta/2., max_noise], ymargin=[3.5, 0], /ylog, thick=3, color=0
; shading
loadct, 0, /silent ; greyscale
for i = 0, n_elements(phi)-1 do begin
; oplot, [phi(i), phi(i)], [zeta/2., noise_antipodal_0(i)], thick=2, color=245
; oplot, [phi(i), phi(i)], [zeta/2., noise_antipodal(i)], thick=2, color=240
; oplot, [phi(i), phi(i)], [zeta/2., noise_0(i)], thick=2, color=235
; oplot, [phi(i), phi(i)], [zeta/2., noise(i)], thick=2, color=230
; oplot, [phi(i), phi(i)], [zeta/2., noise_single_0(i)], thick=2, color=225
 oplot, [phi(i), phi(i)], [zeta/2., noise_single(i)], thick=2, color=220
 oplot, [phi(i), phi(i)], [zeta/2., noise_twin_0(i)], thick=2, color=150
 oplot, [phi(i), phi(i)], [zeta/2., noise_twin(i)], thick=2, color=100
endfor
; replot
oplot, phi, noise_twin, thick=3, color=100
;oplot, phi, noise_twin, color=10
oplot, phi, noise_twin_0, thick=2, color=150
oplot, phi, noise_single, thick=3, color=220
;oplot, phi, noise_single_0, thick=2, color=225
;oplot, phi, noise, thick=3, color=230
oplot, phi, noise_0, thick=2, color=235
;oplot, phi, noise_antipodal, thick=3, color=240
;oplot, phi, noise_antipodal_0, thick=2, color=245
oplot, phi, noise_antipodal, thick=3, color=255 ; color=0
; angle limits, with coloured delinating lines for those regions
loadct, 39, /silent ; colour, angle limit
oplot, [separation_limit_truly_acausal, separation_limit_truly_acausal], [zeta/2., max_noise], thick=5, color=220 ; red
oplot, [separation_limit_truly_acausal-1., separation_limit_truly_acausal-1.], [zeta/2., max_noise], thick=5, color=150 ; green
oplot, [90., 90.], [zeta/2., max_noise], thick=5, color=50 ; blue
;arrow, 90., 0.1, 85., 0.1, thick=2, color=50, /solid, /data
; limits
loadct, 0, /silent ; greyscale
oplot, [182., 182.], [zeta/2., max_noise], thick=25, color=75 ; dark grey, beyond horizon
;oplot, [-5., 185.], [epsilon, epsilon], linestyle=3, color=0
oplot, [-5., 185.], [mean(noise), mean(noise)], linestyle=2, color=0
oplot, [-5., 185.], [median(noise), median(noise)], linestyle=3, color=0
oplot, [-5., 185.], [zeta, zeta], linestyle=1, color=0
oplot, [-5., 185.], [zeta, zeta], color=0
; labels
xyouts, 10., 0.5, 'Single site (one or two tels.)', charsize=charsize, charthick=charthick, color=0
xyouts, phi(min(where(noise_single gt 0.5))), 0.5, '|', orientation=90., charsize=charsize, charthick=charthick, color=0
xyouts, 100., 0.015, 'Twin sites (90-deg sep.)', charsize=charsize, charthick=charthick, color=0
;xyouts, 100., 0.015, '|', orientation=-90., charsize=charsize, charthick=charthick, color=0
xyouts, 107., 0.5, 'Antipodal sites', charsize=charsize, charthick=charthick, color=0
xyouts, phi(max(where(noise_antipodal gt 0.55))), 0.55, '|', orientation=-90., charsize=charsize, charthick=charthick, color=0
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[zeta/2., max_noise], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[zeta/2., max_noise], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_noise_and_signal.jpg', screen_output, quality=100, true=1

; signal to noise ratio
;ratio = (1./zeta) ; by default, as well as the zeropoint
;ratio_0 = 0.5 ; otherwise
;ratio = sqrt(samples)*abs(signal)/noise
ratio = abs(signal)/noise ; instead, restricted to a single observation
ratio_0 = abs(signal)/noise_0 ; single observation
ratio_prime = abs(signal_prime)/noise
ratio_prime_0 = abs(signal_prime)/noise_0 ; single observation
ratio_single_flat = sqrt(samples)*abs(signal_flat)/noise_single ; linear case
;ratio_single_flat = samples*abs(signal_flat)/noise_single ; linear case
ratio_single = sqrt(samples)*abs(signal)/noise_single
;ratio_single = samples*abs(signal)/noise_single
ratio_single_0 = abs(signal)/noise_single_0 ; single observation
ratio_twin_flat = sqrt(samples)*abs(signal_flat)/noise_twin ; linear case
;ratio_twin_flat = samples*abs(signal_flat)/noise_twin ; linear case
ratio_twin = sqrt(samples)*abs(signal)/noise_twin
;ratio_twin = samples*abs(signal)/noise_twin
ratio_twin_prime = sqrt(samples)*abs(signal_prime)/noise_twin
;ratio_twin_prime = samples*abs(signal_prime)/noise_twin
ratio_twin_0 = abs(signal)/noise_twin_0 ; single observation 
ratio_antipodal_flat = sqrt(samples)*abs(signal_flat)/noise_antipodal ; linear case
;ratio_antipodal_flat = samples*abs(signal_flat)/noise_antipodal ; linear case
ratio_antipodal = sqrt(samples)*abs(signal)/noise_antipodal
;ratio_antipodal = samples*abs(signal)/noise_antipodal
ratio_antipodal_0 = abs(signal)/noise_antipodal_0 ; single observation
ratio_calibration = sqrt(samples)*abs(calibration*signal_flat)/noise_twin
;ratio_calibration = samples*abs(calibration*signal_flat)/noise_twin
ratio_calibration_0 = abs(calibration*signal_flat)/noise_twin_0 ; single observation
;ratio_calibration = sqrt(samples)*abs(calibration*signal_fair)/noise_twin
;;ratio_calibration = samples*abs(calibration*signal_fair)/noise_twin
;ratio_calibration_0 = abs(calibration*signal_fair)/noise_twin_0 ; single observation

; lower limit wanted
;ratio(where(ratio lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_0(where(ratio_0 lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
ratio_single_flat(where(ratio_single_flat lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
ratio_single(where(ratio_single lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_single_0(where(ratio_single_0 lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_twin_flat(where(ratio_twin_flat lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_twin(where(ratio_twin lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_twin_0(where(ratio_twin_0 lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_antipodal_flat(where(ratio_antipodal_flat lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_antipodal(where(ratio_antipodal lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
;ratio_antipodal_0(where(ratio_antipodal_0 lt ratio_wanted)) = 1. ; effectively not useful if lower than wanted
; formal lower limit
ratio(where(ratio lt 1.)) = 1. ; formally, cannot be less than unity
ratio_0(where(ratio_0 lt 1.)) = 1. ; formally, cannot be less than unity
ratio_single_flat(where(ratio_single_flat lt 1.)) = 1. ; formally, cannot be less than unity
ratio_single(where(ratio_single lt 1.)) = 1. ; formally, cannot be less than unity
ratio_single_0(where(ratio_single_0 lt 1.)) = 1. ; formally, cannot be less than unity
ratio_twin_flat(where(ratio_twin_flat lt 1.)) = 1. ; formally, cannot be less than unity
ratio_twin(where(ratio_twin lt 1.)) = 1. ; formally, cannot be less than unity
ratio_twin_0(where(ratio_twin_0 lt 1.)) = 1. ; formally, cannot be less than unity
ratio_antipodal_flat(where(ratio_antipodal_flat lt 1.)) = 1. ; formally, cannot be less than unity
ratio_antipodal(where(ratio_antipodal lt 1.)) = 1. ; formally, cannot be less than unity
ratio_antipodal_0(where(ratio_antipodal_0 lt 1.)) = 1. ; formally, cannot be less than unity
ratio_calibration(where(ratio_calibration lt 1.)) = 1. ; formally, cannot be less than unity
ratio_calibration_0(where(ratio_calibration_0 lt 1.)) = 1. ; formally, cannot be less than unity

; signal-to-noise ratio, ideal
; ideal case
ratio_limit_ideal = sqrt(samples)*omega/zeta
print, 'Idealized S/N, as in space: ', ratio_limit_ideal

; maximum signal-to-noise
max_ratio = max([ratio, ratio_prime, ratio_single, ratio_single_flat, ratio_twin, ratio_twin_flat, ratio_antipodal, ratio_antipodal_flat, ratio_calibration])
print, 'Maximum S/N: ', max_ratio

; read in the real test data, or make some
;readcol, 'correlations_test_hold', format='F, F, F, F', colours_test_min, colours_test_max, correlations_test_min, correlations_test_max, /silent
;readcol, 'signals_test_hold', format='F, F, F, F', snrs_test_min, snrs_test_median, snrs_test_mean, snrs_test_max, /silent
readcol, 'correlations_test', format='F, F, F, F', colours_test_min, colours_test_max, correlations_test_min, correlations_test_max, /silent
readcol, 'signals_test', format='F, F, F, F', snrs_test_min, snrs_test_median, snrs_test_mean, snrs_test_max, /silent
groups = 4 ; there are 4 groups of observations
group = findgen(groups) ; a number
group_label = strarr(groups) ; numbering or labelling them
observations = 3 ; there are 3 observations per group
colour_test_min = fltarr(groups)
colour_test_max = colour_test_min
correlation_test_min = colour_test_min
correlation_test_max = colour_test_min
snr_test_min = colour_test_min
snr_test_max = colour_test_min
snr_test_median = colour_test_min
snr_test_mean = colour_test_min
for j = 0, groups-1 do begin ; now parse the observations into groups, in order to average them
 group(j) = j+1
 group_label(j) = strtrim(string(j+1), 2) ;'Obs. '+strtrim(string(j+1), 2)
 colour_test_min(j) = mean(colours_test_min(observations*j: observations*j + observations-1))
 colour_test_max(j) = mean(colours_test_max(observations*j: observations*j + observations-1))
 correlation_test_min(j) = mean(correlations_test_min(observations*j: observations*j + observations-1))
 correlation_test_max(j) = mean(correlations_test_max(observations*j: observations*j + observations-1))
 snr_test_min(j) = mean(snrs_test_min(observations*j: observations*j + observations-1))
 snr_test_max(j) = mean(snrs_test_max(observations*j: observations*j + observations-1))
 snr_test_median(j) = mean(snrs_test_median(observations*j: observations*j + observations-1))
 snr_test_mean(j) = mean(snrs_test_mean(observations*j: observations*j + observations-1)) 
 print, group(j), colour_test_min(j), colour_test_max(j), correlation_test_min(j), correlation_test_max(j), snr_test_min(j), snr_test_max(j), snr_test_median(j), snr_test_mean(j)
endfor
readcol, 'list_test', format='A, A, A, A, A, A, A, A, F', test_a, test_b, ra_test_a, dec_test_a, ra_test_b, dec_test_b, test_date, test_time, ratio_test, /silent ;, skipline=1
tests = n_elements(test_a)
separation_test = fltarr(tests) ; degrees
airmass_test_a = separation_test
airmass_test_b = separation_test
zenith_test_a = separation_test
zenith_test_b = separation_test
orientation_test_a = separation_test
orientation_test_b = separation_test
epoch_test = double(separation_test) ; MJD
;ratio_test = snr_test_min
;ratio_test = snr_test_max
;ratio_test = snr_test_median
ratio_test = snr_test_mean
ratio_test = sqrt(samples)*snr_test_mean
; reset the test ratio to be a measured one, based on correlation measurements
;ratio_test = 1./colour_test_min
;ratio_test = 1./colour_test_max
;ratio_test = 1./correlation_test_min
;ratio_test = 1./correlation_test_max
for i = 0, tests-1 do begin
 ; object names
 print, 'Star A and B (North/South): ', test_a(i), ' ', test_b(i)

 ; convert input r.a. (hh:mm:ss) and dec. (deg:mm:ss) to decimal degrees
 ; A, North
 ra_deg = floor(360.*float(strmid(ra_test_a(i), 0, 2))/24.) ; degrees
 ra_mm = strmid(ra_test_a(i), 3, 2)
 ra_ss = strmid(ra_test_a(i), 6, 5)
 ra_test_a(i) =  tenv(ra_deg, ra_mm, ra_ss) ; degrees
 dec_deg = float(strmid(dec_test_a(i), 0, 3)) ; degrees, noting that the + or - is also converted, if done this way
 dec_mm = strmid(dec_test_a(i), 4, 2)
 dec_ss = strmid(dec_test_a(i), 7, 4)
 dec_test_a(i) =  tenv(dec_deg, dec_mm, dec_ss) ; degrees
 print, 'A, RA and dec (degrees): ', ra_test_a(i), dec_test_a(i)
 ; B, South
 ra_deg = floor(360.*float(strmid(ra_test_b(i), 0, 2))/24.) ; degrees
 ra_mm = strmid(ra_test_b(i), 3, 2)
 ra_ss = strmid(ra_test_b(i), 6, 5)
 ra_test_b(i) =  tenv(ra_deg, ra_mm, ra_ss) ; degrees
 dec_deg = float(strmid(dec_test_b(i), 0, 3)) ; degrees, noting that the + or - is also converted, if done this way
 dec_mm = strmid(dec_test_b(i), 4, 2)
 dec_ss = strmid(dec_test_b(i), 7, 4)
 dec_test_b(i) =  tenv(dec_deg, dec_mm, dec_ss) ; degrees
 print, 'B, RA and dec (degrees): ', ra_test_b(i), dec_test_b(i)
 
 ; convert test time to components
 test_year = strmid(test_date(i), 0, 4)
 test_month = strmid(test_date(i), 5, 2)
 test_day = strmid(test_date(i), 8, 2)
 test_hour = strmid(test_time(i), 0, 2)
 test_minute = strmid(test_time(i), 3, 2)
 test_second = strmid(test_time(i), 6, 2)
; print, test_year, ' ', test_month, ' ', test_day, ' ', test_hour, ' ', test_minute, ' ', test_second
 
 ; and report the MJD epoch of the test
 epoch = julday(test_month, test_day, test_year, test_hour, test_minute, test_second) - 2400000.5 ; MJD
 epoch_test(i) = epoch
 print, 'Epoch of test (MJD): ', epoch_test(i)
 caldat, epoch_test(i) + 2400000.5, month, day, year, hour, minute, second ; confirm, by converting the Julian night back to month and day UT format
 print, 'On calendar date: ', day, month, year, ' at ', hour, ':', minute, ':', second, ' UT time'

 ; from North
 eq2hor, ra_test_a(i), dec_test_a(i), epoch + 2400000.5, alt_test, az_test, lat=latitude_north, lon=longitude_north, altitude=altitude_north ; epoch in Julian date
 airmass_test_a(i) = airmass_simple(90.-alt_test)
 zenith_test_a(i) = 90.-alt_test ; zenith angle in degrees
 orientation_test_a(i) = az_test ; azimuth angle in degrees East of North
 ; from South
; eq2hor, ra_test_b(i), dec_test_b(i), epoch + 2400000.5, alt_test, az_test, lat=latitude_south, lon=longitude_south, altitude=altitude_south ;, /ws ; epoch in Julian date
 eq2hor, ra_test_b(i), dec_test_b(i), epoch + 2400000.5, alt_test, az_test, lat=latitude_south, lon=-longitude_south, altitude=altitude_south ;, /ws ; epoch in Julian date; kludge: needs to be set to East long.
 airmass_test_b(i) = airmass_simple(90.-alt_test)
 zenith_test_b(i) = 90.-alt_test ; zenith angle in degrees
 orientation_test_b(i) = az_test ; azimuth angle in degrees East of North

 ; and find their separation on the sky
 gcirc, 2, ra_test_a(i), dec_test_a(i), ra_test_b(i), dec_test_b(i), separation ; arcsec
 separation_test(i) = separation/60./60. ; degrees
 print, 'Airmass of A and B: ', airmass_test_a(i), airmass_test_b(i)
 print, 'Zenith angle A and B (degrees): ', zenith_test_a(i), zenith_test_b(i) 
 print, 'Orientation A and B (degrees): ', orientation_test_a(i), orientation_test_b(i) 
 print, 'Separation of stars A and B (degrees): ', separation_test(i)
endfor

; set up the display
window, 13, xsize = 1000, ysize = 900, title='Expected Signal-to-Noise Ratio'
wset, 13
!p.multi = [0, 1, 1]
charsize=2.5

; plot the expected signal-to-noise ratio
plot, phi, ratio, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Expected signal-to-noise ratio', charsize=charsize, ystyle=1, yrange=[min(ratio)/2., 1.5*max_ratio], /ylog, color=0
;plot, phi, ratio, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Expected signal-to-noise ratio', charsize=charsize, ystyle=1, yrange=[ratio_wanted, 1.5*max_ratio], /ylog, color=0
;plot, phi, 0.*ratio, xstyle=1, xrange=[-5., 185.], xtitle='Target separation angle (deg)', ytitle='Expected signal-to-noise ratio', charsize=charsize, ystyle=1, yrange=[0., 1.1*max_ratio], color=0
; shading
loadct, 0, /silent ; greyscale
for i = 0, n_elements(phi)-1 do begin
 oplot, [phi(i), phi(i)], [ratio_calibration(i), max_ratio], thick=2, color=220 ; grey
; oplot, [phi(i), phi(i)], [1.e-4, ratio_prime(i)], thick=2, color=220 ; light grey
endfor
; colour shading
loadct, 39, /silent ; colour, angle limit
for i = 0, n_elements(phi)-1 do begin
 oplot, [phi(i), phi(i)], [1.e-4, ratio_antipodal_flat(i)], thick=2, color=190 ; yellow
; oplot, [phi(i), phi(i)], [1.e-4, ratio_antipodal(i)], thick=2, color=190 ; yellow
 oplot, [phi(i), phi(i)], [1.e-4, ratio_twin_flat(i)], thick=2, color=150 ; green
; oplot, [phi(i), phi(i)], [1.e-4, ratio_twin(i)], thick=2, color=150 ; green
 oplot, [phi(i), phi(i)], [1.e-4, ratio_single_flat(i)], thick=2, color=50 ; blue
; oplot, [phi(i), phi(i)], [1.e-4, ratio_single(i)], thick=2, color=50 ; blue
 if phi(i) ge separation_limit_truly_acausal then oplot, [phi(i), phi(i)], [1.e-4, ratio_antipodal_flat(i)], thick=2, color=220 ; red
 if phi(i) ge separation_limit_truly_acausal then oplot, [phi(i), phi(i)], [1.e-4, ratio_antipodal(i)], thick=2, color=220 ; red 
endfor
;oplot, [182., 182.], [1.e-4, 2*max_ratio], thick=25, color=220 ; red
; shading
loadct, 0, /silent ; greyscale
for i = 0, n_elements(phi)-1 do begin
; oplot, [phi(i), phi(i)], [ratio_calibration(i), max_ratio], thick=2, color=220 ; grey
 oplot, [phi(i), phi(i)], [1.e-4, ratio_prime(i)], thick=2, color=220 ; light grey
endfor
; replot
oplot, phi, ratio_calibration, linestyle=2, color=0
oplot, phi, ratio_antipodal, thick=3, color=255
;oplot, phi, ratio_single, thick=3, color=220
oplot, phi, ratio_twin, thick=3, color=100
oplot, phi, ratio_prime, color=0
;oplot, phi, ratio_prime_0, color=0
;oplot, phi, ratio, color=0
;oplot, phi, ratio_0, color=0
;oplot, phi, ratio_single_flat, color=0
;oplot, phi, ratio_single, color=0
;oplot, phi, ratio_single_0, color=0
;oplot, phi, ratio_twin_flat, color=0
;oplot, phi, ratio_twin, color=0
;oplot, phi, ratio_twin_prime, color=0
;oplot, phi, ratio_twin_0, color=0
;oplot, phi, ratio_antipodal_flat, color=0
;oplot, phi, ratio_antipodal, color=0
;oplot, phi, ratio_antipodal_0, color=0
; limits
; re-shading the lower S/N limit
for i = 0, n_elements(phi)-1 do begin
 oplot, [phi(i), phi(i)], [1.e-4, 1.], thick=2, color=100 ; light grey, below S/N threshold
endfor
oplot, [-5., 185.], [1., 1.], color=0
oplot, [182., 182.], [1.e-4, 2.*max_ratio], thick=25, color=75 ; dark grey, beyond horizon
;oplot, [-5., 185.], [ratio_limit_ideal, ratio_limit_ideal], color=0
oplot, [-5., 185.], [max_ratio, max_ratio], color=0
;oplot, [-5., 185.], [max_ratio/2., max_ratio/2.], linestyle=2, color=0 ; shown as the nominal mean
;oplot, [-5., 185.], [max_ratio/2., max_ratio/2.], linestyle=3, color=0 ; shown as the nominal median
;oplot, [-5., 185.], [mean(ratio), mean(ratio)], linestyle=2, color=0
;oplot, [-5., 185.], [median(ratio), median(ratio)], linestyle=3, color=0
;oplot, [-5., 185.], [mean(ratio_twin), mean(ratio_twin)], linestyle=2, color=0
oplot, [-5., 185.], [median(ratio_twin), median(ratio_twin)], linestyle=3, color=0
oplot, [-5., 185.], [ratio_wanted, ratio_wanted], linestyle=1, color=0
;oplot, [-5., 185.], [1., 1.], linestyle=1, color=0
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 185.], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 185.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[min(ratio)/2., 1.5*max_ratio], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[min(ratio)/2., 1.5*max_ratio], ytickformat='(A1)'
;axis, yaxis=0, color=0, ystyle=1, yrange=[ratio_wanted, 1.5*max_ratio], charsize=charsize
;axis, yaxis=1, color=0, ystyle=1, yrange=[ratio_wanted, 1.5*max_ratio], ytickformat='(A1)'
;axis, yaxis=0, color=0, ystyle=1, yrange=[0., 1.1*max_ratio], charsize=charsize
;axis, yaxis=1, color=0, ystyle=1, yrange=[0., 1.1*max_ratio], ytickformat='(A1)'

; test out the plotting symbols
if just_model eq 'yes' then begin

; default, generic or same-field case
phi_peak = 0.
ratio_peak = max(ratio(where(phi le phi_peak)))
plotsym, 0, 2., /fill, color=200 ; filled grey circle
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0
;plotsym, 0, 2., /fill, color=0 ; peak, effect; best case, filled black circle
;oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0

; for Gemini
plotsym, 0, 2., /fill, color=0 ; filled black circle, generic case
phi_peak = 0.
ratio_peak = max(ratio_twin(where(phi le phi_peak)))
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0
plotsym, 0, 3., thick=2, color=200 ; swappable, grey open circle
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0
plotsym, 8, 3., thick=2, color=0 ; an open square, for same-field case
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0
phi_peak = 90.
ratio_peak = max(ratio_twin(where(phi ge phi_peak)))
plotsym, 4, 3., thick=2, color=200 ; generic case; peak beyond horizon, up-pointing triangle
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0
plotsym, 3, 3., thick=2, color=0 ; test case; a star
oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0

; truly acausal limit
phi_peak = phi(where(phi gt separation_limit_truly_acausal))
ratio_peak = ratio_twin(where(phi gt separation_limit_truly_acausal))
plotsym, 5, 3., thick=2, color=200 ; truly acausal, down-pointing triangle
oplot, [phi_peak(where(ratio_peak eq max(ratio_peak))), phi_peak(where(ratio_peak eq max(ratio_peak)))], [max(ratio_peak), max(ratio_peak)], psym=8, color=0

; just stop
print, 'Just testing the symbols, and stopping.'
stop
endif ; if just testing symbols

; go through the desired plots
for plots = 0, n_elements(plottings)-1 do begin
plotting = plottings(plots)
plotsym, 0, 2., /fill, color=0 ; generic case, black

; plotting target positions on the whole sky
if plotting eq 'sky' then begin

; set up the display
window, 14, xsize = 1200, ysize = 700, xpos = 350, ypos = 200, title='Both Objects on the  Sky'
wset, 14
!p.multi = [0, 1, 1]
charsize=2.5

; create some random, fake points
ra_test = 360.*(randomu(seed, 100)) ; degrees, from zero to 180 degrees in right ascension
dec_test = 90.*(1. - 2.*randomu(seed, 100)) ; degrees, from -90 to 90 degrees in declination
ra_plot = 360.*float(findgen(100))/100.+2.

; plot both objects on the same sky; note that the plot retains the scaling of the previous two, to allow polar overplotting to work
;plotsym, 0, 2., /fill, color=0 ; filled black circle, generic case
plotsym, 0, 2., /fill, color=255 ; kludge, filled white circle, test case made invisible in plot
plot, ra_test, dec_test, xstyle=1, ystyle=1, xrange=[-5., 365.], yrange=[-95, 95.], xtitle='Right Ascension (degrees)', ytitle='Declination (degrees)', charsize=charsize, psym=8, color=0
; shading
for i = 0, 99 do begin
oplot, [ra_plot(i), ra_plot(i)], [latitude_north + zenith_limit, latitude_north + zenith_limit + 90.], thick=35, color=220
oplot, [ra_plot(i), ra_plot(i)], [latitude_south - zenith_limit, latitude_south - zenith_limit - 90.], thick=35, color=220
endfor
; limit
oplot, [-5., 365.], [0., 0.], linestyle=3, color=0 ; equator
oplot, [180., 180.], [-95., 95.], linestyle=1, color=0 ; a dotted line delineating right ascension of 180 degrees
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[-5., 365.], xtickformat='(A1)' 
axis, xaxis=1, color=0, xstyle=1, xrange=[-5., 365.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-95., 95.], ytickformat='(A1)'
axis, yaxis=1, color=0, ystyle=1, yrange=[-95., 95.], ytickformat='(A1)'

endif ; if plotting both objects on the sky

; plotting the orientation of sources
if plotting eq 'polar' then begin

; orientation plots setup
airmass_circle = replicate(airmass_limit-1., 1000)
zenith_circle = replicate(90., 1000)
limit_circle = replicate(zenith_limit, 1000)
orientation_circle = 180.*float(findgen(1000))/1000.
;airmass_test = (airmass_limit-1.)*randomu(seed, 18)
;orientation_test = 10.*findgen(18) ; degrees
;airmass_test = [1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2., 2.1, 2.2, 2.3, 2.4, 2.5]-1. ; in airmass, starting at 1
zenith_test = [0., 10., 20., 30., 40., 50., 60., 70., 80., 90. ] ; in degrees
orientation_zero = 90. ; 270. ; 180. ; 90. ; 0. ; degrees, where 0 is North, and increasing towards the East; so due-East is 90 degrees, south is 180 degrees, and west is 270 degrees
;orientation_test = replicate(orientation_zero, n_elements(airmass_test))+90. ; in degrees, rotated 90 degrees, as the plotting function has zero-angle off towards the right
orientation_test = replicate(orientation_zero, n_elements(zenith_test))+90. ; in degrees, rotated 90 degrees, as the plotting function has zero-angle off towards the right
orientation_test = (!pi/180.)*orientation_test ; in radians

; set up the display
window, 15, xsize = 750, ysize = 700, xpos = 1200, ypos = 500, title='Orientation of Target On Sky: North'
wset, 15
!p.multi = [0, 1, 1]
charsize=2.5

; plot orientation of target in North
;plot, airmass_circle, orientation_circle, /polar, xstyle=1, ystyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], yrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)', ytickformat='(A1)', xtitle='Target of Pair as Viewed from Gemini North', charsize=charsize, color=0
plot, zenith_circle, orientation_circle, /polar, xstyle=1, ystyle=1, xrange=[-90., 90.], yrange=[-90., 90.], xtickformat='(A1)', ytickformat='(A1)', xtitle='Target of Pair as Viewed from Gemini North', charsize=charsize, color=0
; shading
for i = 0, 999 do begin
 oplot, (1.+float(i)/500.)*limit_circle, orientation_circle, /polar, color=220 ; outside the airmass limit
endfor
for i = 0, 999 do begin
; oplot, (1.+float(i)/500.)*airmass_circle, orientation_circle, /polar, color=100 ; outside the airmass limit
 oplot, (1.+float(i)/500.)*zenith_circle, orientation_circle, /polar, color=100 ; outside the horizon limit
endfor
;oplot, airmass_circle, orientation_circle, /polar, color=0 ; a circle at highest possible airmass
;oplot, zenith_circle, orientation_circle, /polar, color=0 ; a circle at highest possible zenith angle
;oplot, limit_circle, orientation_circle, /polar, color=0 ; a circle at airmass limit
;oplot, airmass_test, orientation_test, /polar, color=0 ; a line towards the East
oplot, zenith_test, orientation_test, /polar, color=0 ; a line towards the East

; label
;xyouts, -(airmass_limit-1.), (!pi/180.)*180., 'East', charsize=charsize, color=0, /device
xyouts, 50, 375, 'Eastern', charsize=charsize, color=0, /device
xyouts, 50, 355, 'horizon', charsize=charsize, color=0, /device
xyouts, 400, 600, 'North', charsize=charsize, color=0, /device
; clean up
;axis, xaxis=0, color=0, xstyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)' 
;axis, xaxis=1, color=0, xstyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)'
;axis, yaxis=0, color=0, ystyle=1, yrange=[-(airmass_limit-1.), airmass_limit-1.], ytickformat='(A1)'
;axis, yaxis=1, color=0, ystyle=1, yrange=[-(airmass_limit-1.), airmass_limit-1.], ytickformat='(A1)'
axis, xaxis=0, color=0, xstyle=1, xrange=[-90., 90.], xtickformat='(A1)' 
axis, xaxis=1, color=0, xstyle=1, xrange=[-90., 90.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-90., 90.], ytickformat='(A1)'
axis, yaxis=1, color=0, ystyle=1, yrange=[-90., 90.], ytickformat='(A1)'

; set up the display
window, 16, xsize = 750, ysize = 700, xpos = 0, ypos = 0, title='Orientation of Target On Sky: South'
wset, 16
!p.multi = [0, 1, 1]
charsize=2.5

; plot orientation of target in South
;plot, airmass_circle, orientation_circle, /polar, xstyle=1, ystyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], yrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)', ytickformat='(A1)', xtitle='Target of Pair as Viewed from Gemini South', charsize=charsize, color=0
plot, zenith_circle, orientation_circle, /polar, xstyle=1, ystyle=1, xrange=[-90., 90.], yrange=[-90., 90.], xtickformat='(A1)', ytickformat='(A1)', xtitle='Target of Pair as Viewed from Gemini South', charsize=charsize, color=0
; shading
for i = 0, 999 do begin
 oplot, (1.+float(i)/500.)*limit_circle, orientation_circle, /polar, color=220 ; outside the airmass limit
endfor
for i = 0, 999 do begin
; oplot, (1.+float(i)/500.)*airmass_circle, orientation_circle, /polar, color=100 ; outside the airmass limit
 oplot, (1.+float(i)/500.)*zenith_circle, orientation_circle, /polar, color=100 ; outside the horizon limit
endfor
;oplot, airmass_circle, orientation_circle, /polar, color=0 ; a circle at highest possible airmass
;oplot, zenith_circle, orientation_circle, /polar, color=0 ; a circle at highest possible zenith angle
;oplot, limit_circle, orientation_circle, /polar, color=0 ; a circle at airmass limit
;oplot, airmass_test, orientation_test, /polar, color=0 ; a line towards the East
oplot, zenith_test, orientation_test, /polar, color=0 ; a line towards the East
; label
;xyouts, -(airmass_limit-1.), (!pi/180.)*180., 'East', charsize=charsize, color=0, /device
xyouts, 50, 375, 'Eastern', charsize=charsize, color=0, /device
xyouts, 50, 355, 'horizon', charsize=charsize, color=0, /device
;xyouts, 400, 600, 'North', charsize=charsize, color=0, /device
xyouts, 400, 135, 'South', charsize=charsize, color=0, /device
; clean up
;axis, xaxis=0, color=0, xstyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)' 
;axis, xaxis=1, color=0, xstyle=1, xrange=[-(airmass_limit-1.), airmass_limit-1.], xtickformat='(A1)'
;axis, yaxis=0, color=0, ystyle=1, yrange=[-(airmass_limit-1.), airmass_limit-1.], ytickformat='(A1)'
;axis, yaxis=1, color=0, ystyle=1, yrange=[-(airmass_limit-1.), airmass_limit-1.], ytickformat='(A1)'
axis, xaxis=0, color=0, xstyle=1, xrange=[-90., 90.], xtickformat='(A1)' 
axis, xaxis=1, color=0, xstyle=1, xrange=[-90., 90.], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-90., 90.], ytickformat='(A1)'
axis, yaxis=1, color=0, ystyle=1, yrange=[-90., 90.], ytickformat='(A1)'

endif ; if plotting orientations and directions on sky

; start by plotting the default, generic or same-field case
phi_peak = 0.
ratio_peak = 2.*1./zeta ; max(ratio(where(phi le phi_peak))) ; the best possible SNR
;print, ratio_peak
;plotsym, 0, 2., /fill, color=200 ; filled grey circle
;oplot, [phi_peak, phi_peak], [ratio_peak, ratio_peak], psym=8, color=0

; go through all scenarios
for n = 0, n_elements(experiments) - 1 do begin
experiment = experiments(n)

; settings
if experiment eq 'generic' then begin
 separation_limit = separation_limit_generic
 field_limit = field_limit_generic
 magnitude_limit = magnitude_limit_generic
 colour_limit = colour_limit_generic
 redshift_limit = redshift_limit_generic
 separation_star = separation_star_generic
 print, '****************************************************************************************************************************************************************************'
 print, 'Generic, sky-wide case.'
endif
if experiment eq 'peak' then begin
 separation_limit = separation_limit_peak
 field_limit = field_limit_peak
 magnitude_limit = magnitude_limit_peak
 colour_limit = colour_limit_peak
 redshift_limit = redshift_limit_peak
 separation_star = separation_star_peak
 print, '****************************************************************************************************************************************************************************'
 print, 'Checking near expected peak.'
endif
if experiment eq 'truly_acausal' then begin
 separation_limit = separation_limit_truly_acausal
 field_limit = field_limit_truly_acausal
 magnitude_limit = magnitude_limit_truly_acausal
 colour_limit = colour_limit_truly_acausal
 redshift_limit = redshift_limit_truly_acausal
 separation_star = separation_star_truly_acausal
 print, '****************************************************************************************************************************************************************************'
 print, 'Demanding truly acausal photons.'
endif
if experiment eq 'same_field' then begin
 separation_limit = separation_limit_same_field
 field_limit = field_limit_same_field
 magnitude_limit = magnitude_limit_same_field
 colour_limit = colour_limit_same_field
 redshift_limit = redshift_limit_same_field
 separation_star = separation_star_same_field
 print, '****************************************************************************************************************************************************************************'
 print, 'Requiring pair visible in the same field.'
endif
if semester eq 'A' then begin
 night_0 = julday(2, 1, calendar_year, 0, 0, 0) - 2400000.5 ; MJD of 1 February of a given year
 nights = floor(365.25/2.)
endif
if semester eq 'B' then begin
 night_0 = julday(8, 1, calendar_year, 0, 0, 0) - 2400000.5 ; MJD of 1 August of a given year
 nights = floor(365.25/2.)
endif
if semester eq 'AB' then begin
 night_0 = julday(1, 1, calendar_year, 0, 0, 0) - 2400000.5 ; MJD of 1 second after midnight UTC of 1 January of a given year
 nights = floor(365.25)
endif
print, 'Magnitude limit: ', magnitude_limit
print, 'Colour limit (mag): ', colour_limit
print, 'Redshift limit: ', redshift_limit
print, 'Minimum nights visible: ', minimum_nights

; positions of R-Johnson<18.5 z>3.65 QSOs in MILLIQUAS
;readcol, 'list_quasars_test', format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift, /silent ;, skipline=1
;readcol, 'list_quasars_zgt3.65_rlt18.5', format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift, /silent ;, skipline=1
;readcol, 'list_quasars_zgt' + strmid(strtrim(string(redshift_limit), 2), 0, 4) + '_rlt' + strmid(strtrim(string(magnitude_limit), 2), 0, 4), format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift
readcol, './list_quasars_z_r/list_quasars_zgt' + strmid(strtrim(string(redshift_limit), 2), 0, 4) + '_rlt' + strmid(strtrim(string(magnitude_limit), 2), 0, 4), format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift, /silent ;, skipline=1
;identification = strtrim(identification, 2)
identifications = n_elements(identification)

; convert input r.a. (hh:mm:ss) and dec. (deg:mm:ss) to decimal degrees, and find quasar colour
;print, ra
ra_deg = floor(360.*float(strmid(ra, 0, 2))/24.) ; degrees
ra_mm = strmid(ra, 3, 2)
ra_ss = strmid(ra, 6, 5)
;print, ra_deg
;print, ra_mm
;print, ra_ss
ra =  tenv(ra_deg, ra_mm, ra_ss) ; degrees
;print, ra
;print, dec
dec_deg = float(strmid(dec, 0, 3)) ; degrees, noting that the + or - is also converted, if done this way
dec_mm = strmid(dec, 4, 2)
dec_ss = strmid(dec, 7, 4)
;print, dec_deg
;print, dec_mm
;print, dec_ss
dec =  tenv(dec_deg, dec_mm, dec_ss) ; degrees
;print, dec
magnitude = red
colour = blue - red

;; test of converting back to sexagesimal
;radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc, /hours
;print, ihr, imin, xsec
;print, ideg, imn, xsc

;; test of reformating back to sexagesimal as a string
;ra_dec = adstring(ra(0), dec(0))
;strput, ra_dec, ':', 3
;strput, ra_dec, ':', 6
;strput, ra_dec, ':', 16
;strput, ra_dec, ':', 19
;print, ra_dec

; go through and find angular separations
identification_a = strarr(identifications^2)
identification_b = identification_a
magnitude_a = fltarr(identifications^2)
magnitude_b = magnitude_a
redshift_a = magnitude_a
redshift_b = magnitude_a
ra_a = magnitude_a
ra_b = magnitude_a
dec_a = magnitude_a
dec_b = magnitude_a
separation = magnitude_a
candidate_pair = 0 ; counter, to tally the number of candidate pairs
for i = 0, identifications - 2 do begin

  ; separation of the calibration and A
  gcirc, 2, ra_c, dec_c, ra(i), dec(i), separation_ac ; arcsec
  separation_ac = separation_ac/60./60. ; degrees
;  print, 'Separation of A and C (degrees): ', separation_ac

  ; for each B target
  for j = i+1, identifications - 1 do begin 

  ; separation of the calibration and B
  gcirc, 2, ra_c, dec_c, ra(j), dec(j), separation_bc ; arcsec
  separation_bc = separation_bc/60./60. ; degrees
;  print, 'Separation of B and C (degrees): ', separation_bc

  ; separation of A and B
  gcirc, 2, ra(i), dec(i), ra(j), dec(j), separation_ab ; arcsec
  separation_ab = separation_ab/60./60. ; degrees
;  print, 'Separation of A and B (degrees): ', separation_ab
 
  ; colour difference
  colour_difference = abs(colour(i) - colour(j))

  ; find out if individual target redshift, brightness, their relative colour and those angles work
  if separation_ac le calibration_limit and separation_bc le calibration_limit and separation_ab gt separation_limit and separation_ab le separation_limit+field_limit and magnitude(i) le magnitude_limit and magnitude(j) le magnitude_limit and colour_difference le colour_limit and redshift(i) gt redshift_limit and redshift(j) gt redshift_limit then begin
  
   ; report
;   print, 'Separation of A and C (degrees): ', separation_ac
;   print, 'Separation of B and C (degrees): ', separation_bc
;   print, 'Separation of A and B (degrees): ', separation_ab
;   print, 'Sep. from C, A to B (degrees): ', separation_ac, separation_bc, separation_ab

   ; fix a bug
   if candidate_pair le 0 then candidate_pair = floor(-1.*candidate_pair) ; needed, as there is a bug which can reset this counter to negative numbers

   ; relabel
   candidate_pair = candidate_pair + 1
   identification_a(candidate_pair) = identification(i)
   identification_b(candidate_pair) = identification(j)
   magnitude_a(candidate_pair) = magnitude(i)
   magnitude_b(candidate_pair) = magnitude(j)
   redshift_a(candidate_pair) = redshift(i)
   redshift_b(candidate_pair) = redshift(j)
   ra_a(candidate_pair) = ra(i)
   ra_b(candidate_pair) = ra(j)
   dec_a(candidate_pair) = dec(i)
   dec_b(candidate_pair) = dec(j)
   separation(candidate_pair) = separation_ab

  endif

 endfor

endfor
candidate_pairs = candidate_pair ; total candidate pairs
print, 'There were ', candidate_pairs, ' candidate pairs found on the sky'

; find the unique target names
target = [identification_a, identification_b]
target = target[uniq(target, sort(target))]
targets = n_elements(target)
print, 'Comprising ', targets-1, ' unique targets'

; check if target name is blank
for k = 1, targets-1 do begin
 if strlen(target(k)) le 1 then begin
  print, 'Target name blank', target(k)
  stop
 endif
endfor

; reorder by separation
order = reverse(sort(separation))
identification_a = identification_a(order)
identification_b = identification_b(order)
magnitude_a = magnitude_a(order)
magnitude_b = magnitude_b(order)
redshift_a = redshift_a(order)
redshift_b = redshift_b(order)
ra_a = ra_a(order)
ra_b = ra_b(order)
dec_a = dec_a(order)
dec_b = dec_b(order)
separation = separation(order)

; kludge: avoiding spurious zero-separation cases for antipodal separations
if experiment eq 'truly_acausal' then separation(where(separation eq 0.)) = !values.f_nan ; delete zeroes, in this case
;print, 'All the separations (degrees): '
;print, separation

; report candidate pairs
print, 'Separation (degrees), magnitudes, redshifts, positions, identifications:'
for k = 0, candidate_pairs-1 do begin
 print, separation(k), magnitude_a(k), magnitude_b(k), redshift_a(k), redshift_b(k), ra_a(k), dec_a(k), ra_b(k), dec_b(k), ' ', identification_a(k), ' ', identification_b(k)
endfor

; now, go through these candidate pairs, and find the optimal night to observe each viable pair, minimizing both their airmasses, and plotting them
print, 'First night of semester A or B, or year (if AB), in MJD: ', night_0 + night_offset
caldat, night_0 + night_offset + 2400000.5, month, day, year, hour, minute, second ; confirm, by converting the Julian night back to month and day UT format
print, 'On calendar date: ', day, month, year, ' at ', hour, ':', minute, ' UT time'
night = night_0 + night_offset + float(findgen(nights)) + 2400000.5 ; this is the Julian date at sunset each night, offset from midnight UTC to be sunset in Hawai'i
airmass_a_north = fltarr(nights) + airmass_limit ; set by default to be just beyond the airmasss limit
airmass_b_north = airmass_a_north
airmass_a_south = airmass_a_north
airmass_b_south = airmass_a_north
best_airmass_north = airmass_limit
best_airmass_south = best_airmass_north
pair = 0
best_ab = 0
best_ba = 0
ratio_direct = fltarr(n_elements(separation))
ratio_causal = ratio_direct
ratio_best = ratio_direct
swapped = 0
print, 'Pair, identifications, redshifts, magnitudes: '
if print_long_info eq 'yes' then print, 'A-B, best airmasses, separation (degrees), month and day, nights visible; star calibrators, r magnitudes and separations in arcseconds:'
for k = 0, candidate_pairs-1 do begin

 ; convert position to sexagesimal format for print output
 ra_dec_a = adstring(ra_a(k), dec_a(k))
 strput, ra_dec_a, ':', 3
 strput, ra_dec_a, ':', 6
 strput, ra_dec_a, ':', 16
 strput, ra_dec_a, ':', 19
 ra_dec_a = strtrim(ra_dec_a, 2)
 ra_dec_b = adstring(ra_b(k), dec_b(k))
 strput, ra_dec_b, ':', 3
 strput, ra_dec_b, ':', 6
 strput, ra_dec_b, ':', 16
 strput, ra_dec_b, ':', 19
 ra_dec_b = strtrim(ra_dec_b, 2)

 ; determine if there are UCAC4 or Gaia stars of "comparable" brightness with which to calibrate photometry, by querying the Vizier catalog
 name_a = 'self-calibrated' ; generic name place-holder, in case no calibration star is wanted
 name_b = 'self-calibrated'
 star_mag_a = 0. ; ' N/A '
 star_mag_b = 0. ; ' N/A '
 offset_a = 0. ; ' 0.0 '
 offset_b = 0. ; ' 0.0 '
 if calibrate eq 'yes' then begin
  ; using UCAC 
  if catalog eq 'ucac' then catname = 'UCAC4'
  if catalog eq 'ucac' then constraint = '' ; 'rmag<20' ; '' ; 'rmag<magnitude_limit' should work, but doesn't ; blank is no limit
  ; using Gaia
  if catalog eq 'gaia' then catname = 'I/345/gaia2' ; GAIA ; note that simply 'GAIA' gets the input catalog, not the actual DR2 data release
  if catalog eq 'gaia' then constraint = '' ; 'Gmag<'+strtrim(string(magnitude_limit), 1) ; 'Gmag>1' ; 'Gmag<20.5' ; '' ; must match the catalog, including capitalization, blank equals no constraint; default is to find a star as bright as the target, or it could be just to find any star
; print, constraint
  info_a = queryvizier(catname, [ra_a(k), dec_a(k)], [separation_star_default/60., separation_star_default/60.], constra=constraint, /silent) ; just finding any star within GMOS field of view
  info_b = queryvizier(catname, [ra_b(k), dec_b(k)], [separation_star_default/60., separation_star_default/60.], constra=constraint, /silent) ; just finding any star within GMOS field of view
;  info_a = queryvizier(catname, [ra_a(k), dec_a(k)], separation_star/60., constra=constraint, /silent) ; in arcminutes
;  info_b = queryvizier(catname, [ra_b(k), dec_b(k)], separation_star/60., constra=constraint, /silent) ; in arcminutes
  record_a = size(info_a)
  record_b = size(info_b)
  if record_a(0) eq 0 or record_b(0) eq 0 then begin ; there is no record for one; so both are unobserved
   ; A star
   name_a = 'unobserved'
   ra_star_a = !values.f_nan
   dec_star_a = !values.f_nan
   star_mag_a = !values.f_nan
   ; B star
   name_b = 'unobserved'
   ra_star_b = !values.f_nan
   dec_star_b = !values.f_nan
   star_mag_b = !values.f_nan
;   print, 'No stars found.'
  endif else begin ; there is a record for both, so record those
;   ; UCAC
;   star_a = min(where(info_a.rmag eq max(info_a.rmag, /nan))) ; faintest
   if catalog eq 'ucac' then star_a = min(where(info_a.rmag eq min(info_a.rmag, /nan))) ; brightest
;   star_b = min(where(info_b.rmag eq max(info_b.rmag, /nan))) ; faintest
   if catalog eq 'ucac' then star_b = min(where(info_b.rmag eq min(info_b.rmag, /nan))) ; brightest
   ; Gaia
;   star_a = min(where(info_a.Gmag eq max(info_a.Gmag, /nan))) ; faintest
   if catalog eq 'gaia' then star_a = min(where(info_a.Gmag eq min(info_a.Gmag, /nan))) ; brightest
;   star_b = min(where(info_b.Gmag eq max(info_b.Gmag, /nan))) ; faintest
   if catalog eq 'gaia' then star_b = min(where(info_b.Gmag eq min(info_b.Gmag, /nan))) ; brightest
   
   ; A star
   if catalog eq 'ucac' then name_a = info_a[star_a].ucac4 ; object name, UCAC
   if catalog eq 'gaia' then name_a = info_a[star_a].source ; object name, Gaia
   name_a = string(name_a)
   if catalog eq 'ucac' then ra_star_a = info_a[star_a].raj2000 ; UCAC
   if catalog eq 'ucac' then dec_star_a = info_a[star_a].dej2000 ; UCAC
   if catalog eq 'gaia' then ra_star_a = info_a[star_a].ra_icrs ; Gaia
   if catalog eq 'gaia' then dec_star_a = info_a[star_a].de_icrs ; Gaia
   if catalog eq 'ucac' then star_mag_a = info_a[star_a].rmag ; UCAC
   if catalog eq 'gaia' then star_mag_a = info_a[star_a].Gmag ; Gaia 
   gcirc, 2, ra_a(k), dec_a(k), ra_star_a, dec_star_a, offset_a ; arcsec
   ; B star
   if catalog eq 'ucac' then name_b = info_b[star_b].ucac4 ; object name, UCAC
   if catalog eq 'gaia' then name_b = info_b[star_b].source ; object name, Gaia
   name_b = string(name_b) 
   if catalog eq 'ucac' then ra_star_b = info_b[star_b].raj2000 ; UCAC
   if catalog eq 'ucac' then dec_star_b = info_b[star_b].dej2000 ; UCAC
   if catalog eq 'gaia' then ra_star_b = info_b[star_b].ra_icrs ; Gaia
   if catalog eq 'gaia' then dec_star_b = info_b[star_b].de_icrs ; Gaia
   if catalog eq 'ucac' then star_mag_b = info_b[star_b].rmag ; UCAC
   if catalog eq 'gaia' then star_mag_b = info_b[star_b].Gmag ; Gaia
   gcirc, 2, ra_b(k), dec_b(k), ra_star_b, dec_star_b, offset_b ; arcsec
  endelse
 endif

 ; if the pair works, carry on
 if name_a ne 'unobserved' and name_b ne 'unobserved' and star_mag_a le magnitude_limit and star_mag_b le magnitude_limit and offset_a le separation_star and offset_b le separation_star then begin ; both stars also as bright as the targets and within their respective fields, noting that separation here is now in arcminutes

  ; calculate airmass of those targets for each night at sunset in Hawai`i, offset from midnight UTC to sunset in Hawai`i on that night
  for l = 0, nights-1 do begin
   ; from North
   eq2hor, ra_a(k), dec_a(k), night(l), alt_a_north, az_a_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north ; night is Julian date
   eq2hor, ra_b(k), dec_b(k), night(l), alt_b_north, az_b_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north ; night is Julian date
   ; from South
   eq2hor, ra_a(k), dec_a(k), night(l), alt_a_south, az_a_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south ; night is Julian date
   eq2hor, ra_b(k), dec_b(k), night(l), alt_b_south, az_b_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south ; night is Julian date
   ; airmass
   airmass_a_north(l) = airmass_simple(90.-alt_a_north)
   airmass_b_north(l) = airmass_simple(90.-alt_b_north)
   airmass_a_south(l) = airmass_simple(90.-alt_a_south)
   airmass_b_south(l) = airmass_simple(90.-alt_b_south)
  endfor

  ; find out if there is at least one night when both targets are below an airmass limit; will return -1 if there is none
;  best_nights_ab = night(where(airmass_a_north lt airmass_limit and airmass_b_south lt airmass_limit))
;  best_nights_ba = night(where(airmass_b_north lt airmass_limit and airmass_a_south lt airmass_limit))
  best_nights_ab = night(where(airmass_a_north lt airmass_limit and airmass_b_south lt airmass_limit and airmass_a_north + airmass_b_south lt 2.*airmass_limit))
  best_nights_ba = night(where(airmass_b_north lt airmass_limit and airmass_a_south lt airmass_limit and airmass_b_north + airmass_a_south lt 2.*airmass_limit))
;  best_nights_ab = night(where(airmass_a_north lt airmass_limit and airmass_a_north + airmass_b_south lt 2.*airmass_limit))
;  best_nights_ba = night(where(airmass_b_north lt airmass_limit and airmass_b_north + airmass_a_south lt 2.*airmass_limit))

  ; start by considering A from North and B from South
  best_airmass_a_north = airmass_limit ; assume right at the visible limit to start
  best_airmass_b_south = best_airmass_a_north
  best_airmass_a_south = best_airmass_a_north
  best_airmass_b_north = best_airmass_a_north 
  if n_elements(best_nights_ab) ge minimum_nights and n_elements(best_nights_ab) ge n_elements(best_nights_ba) then begin ; better, or equal to, observable nights with target A from South and target B from North
   first_night = min(night(where(airmass_a_north + airmass_b_south eq min(airmass_a_north + airmass_b_south)))) ; first MJD of lowest combined airmass for both
   last_night = max(night(where(airmass_a_north + airmass_b_south eq min(airmass_a_north + airmass_b_south)))) ; last MJD of lowest combined airmass for both
;   best_night = first_night ; choose the first night
;   best_night = last_night-minimum_nights ; count back the minimum required number of nights from the last
   best_night = first_night + (last_night - first_night)/2. ; choose the middle
;   best_airmass_a = sigfig(airmass_a_north(where(night eq best_night)), 3) ; airmass at start of that best night
;   best_airmass_b = sigfig(airmass_b_south(where(night eq best_night)), 3) ; airmass at start of that best night  
;   print, best_night
   caldat, best_night, best_month, best_day, year, hour, minute, second ; convert best Julian night to month and day format
   best_ab = best_ab + 1
   pair = pair + 1
   print, pair, ' ', identification_a(k), ' ', redshift_a(k), magnitude_a(k), ' ', identification_b(k), ' ', redshift_b(k), magnitude_b(k)

   ; now, find the best hour within that best night, when combined airmasses are lowest
   night_and_hour = best_night + (0. + findgen(hours))/float(hours)
   airmasses_a_north = fltarr(hours) + airmass_limit ; assuming just at the limit to start
;   airmasses_b_north = airmasses_a_north
;   airmasses_a_south = airmasses_a_north
   airmasses_b_south = airmasses_a_north
   orientations_north = fltarr(hours)
   orientations_south = orientations_north
   zeniths_north = orientations_north
   zeniths_south = orientations_north
    for m = 0, hours - 1 do begin
    ; from North
    eq2hor, ra_a(k), dec_a(k), night_and_hour(m), alt_a_north, az_a_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north
;    eq2hor, ra_b(k), dec_b(k), night_and_hour(m), alt_b_north, az_b_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north
    ; from South
;    eq2hor, ra_a(k), dec_a(k), night_and_hour(m), alt_a_south, az_a_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south
    eq2hor, ra_b(k), dec_b(k), night_and_hour(m), alt_b_south, az_b_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south
    ; airmasses
    airmasses_a_north(m) = airmass_simple(90.-alt_a_north)
;    airmasses_b_north(m) = airmass_simple(90.-alt_b_north)
;    airmasses_a_south(m) = airmass_simple(90.-alt_a_south)
    airmasses_b_south(m) = airmass_simple(90.-alt_b_south)
    ; orientations
    orientations_north(m) = az_a_north
;    orientations_north(m) = az_b_north
;    orientations_south(m) = az_a_south
    orientations_south(m) = az_b_south
    ; zeniths
    zeniths_north(m) = 90. - alt_a_north
;    zeniths_north(m) = 90. - alt_b_north
;    zeniths_south(m) = 90. - alt_a_south
    zeniths_south(m) = 90. - alt_b_south
   endfor
   best_night_and_hour = min(night_and_hour(where(airmasses_a_north + airmasses_b_south eq min(airmasses_a_north + airmasses_b_south)))) ; Julian date of lowest airmass for both, now to the hour
;   best_night_and_hour = min(night_and_hour(where(airmasses_b_north + airmasses_a_south eq min(airmasses_b_north + airmasses_a_south)))) ; Julian date of lowest airmass for both, now to the hour
   best_airmass_a_north = sigfig(airmasses_a_north(where(night_and_hour eq best_night_and_hour)), 3)
   best_airmass_b_south = sigfig(airmasses_b_south(where(night_and_hour eq best_night_and_hour)), 3)
;   best_airmass_a_south = sigfig(airmasses_a_south(where(night_and_hour eq best_night_and_hour)), 3)
;   best_airmass_b_north = sigfig(airmasses_b_north(where(night_and_hour eq best_night_and_hour)), 3)
   orientation_north = -sigfig(orientations_north(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees East of North
   orientation_south = -sigfig(orientations_south(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees East of North
   zenith_north = sigfig(zeniths_north(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees
   zenith_south = sigfig(zeniths_south(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees  
   best_airmass_north = best_airmass_a_north ; simply recording A as the best for north
   best_airmass_south = best_airmass_b_south ; simply recording B as the best for south 
   if print_long_info then begin
    print, ' N-S ', best_airmass_north, ' ', best_airmass_south, separation(k), best_month, best_day,  n_elements(best_nights_ab)
    print, ' IDs ', name_a, star_mag_a, round(offset_a), ' ', name_b, star_mag_b, round(offset_b)
   endif
  endif
  
  ; or instead, for B from North and A from South
  if n_elements(best_nights_ba) ge minimum_nights and n_elements(best_nights_ba) gt n_elements(best_nights_ab) then begin ; otherwise, check if the reverse is better than observable nights with target A from North and target B from South
   first_night = min(night(where(airmass_b_north + airmass_a_south eq min(airmass_b_north + airmass_a_south)))) ; first Julian date of lowest combined airmass for both
   last_night = max(night(where(airmass_b_north + airmass_a_south eq min(airmass_b_north + airmass_a_south)))) ; last Julian date of lowest combined airmass for both  
;   best_night = first_night ; choose the first night
;   best_night = last_night-minimum_nights ; count back the minimum required number of nights from the last
   best_night = first_night + (last_night - first_night)/2. ; choose the middle
;   best_airmass_a = sigfig(airmass_a_south(where(night eq best_night)), 3) ; airmass at start of that best night
;   best_airmass_b = sigfig(airmass_b_north(where(night eq best_night)), 3) ; airmass at start of that best night
;   print, best_night
   caldat, best_night, best_month, best_day, year, hour, minute, second ; convert best Julian night to month and day format
   best_ba = best_ba + 1
   pair = pair + 1
   print, pair, ' ', identification_a(k), redshift_a(k), magnitude_a(k), ' ', identification_b(k), redshift_b(k), magnitude_b(k)

   ; now, find the best hour within that night, when combined airmasses are lowest
   night_and_hour = best_night + (0. + findgen(hours))/float(hours)
   airmasses_a_north = fltarr(hours) + airmass_limit ; assuming just at the limit to start
   airmasses_b_north = airmasses_a_north
   airmasses_a_south = airmasses_a_north
;   airmasses_b_south = airmasses_a_north
   orientations_north = fltarr(hours)
   orientations_south = orientations_north
   zeniths_north = orientations_north
   zeniths_south = orientations_north  
   for m = 0, hours - 1 do begin
    ; from North
;    eq2hor, ra_a(k), dec_a(k), night_and_hour(m), alt_a_north, az_a_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north
    eq2hor, ra_b(k), dec_b(k), night_and_hour(m), alt_b_north, az_b_north, lat=latitude_north, lon=longitude_north, altitude=altitude_north
    ; from South
    eq2hor, ra_a(k), dec_a(k), night_and_hour(m), alt_a_south, az_a_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south
;    eq2hor, ra_b(k), dec_b(k), night_and_hour(m), alt_b_south, az_b_south, lat=latitude_south, lon=longitude_south, altitude=altitude_south
    ; airmasses
;    airmasses_a_north(m) = airmass_simple(90.-alt_a_north)
    airmasses_b_north(m) = airmass_simple(90.-alt_b_north)
    airmasses_a_south(m) = airmass_simple(90.-alt_a_south)
;    airmasses_b_south(m) = airmass_simple(90.-alt_b_south)
    ; orientations
;    orientations_north(m) = az_a_north
    orientations_north(m) = az_b_north
    orientations_south(m) = az_a_south
;    orientations_south(m) = az_b_south
     ; zeniths
;    zeniths_north(m) = 90. - alt_a_north
    zeniths_north(m) = 90. - alt_b_north
    zeniths_south(m) = 90. - alt_a_south
;    zeniths_south(m) = 90. - alt_b_south
   endfor
;   best_night_and_hour = min(night_and_hour(where(airmasses_a_north + airmasses_b_south eq min(airmasses_a_north + airmasses_b_south)))) ; Julian date of lowest airmass for both, now to the hour
   best_night_and_hour = min(night_and_hour(where(airmasses_b_north + airmasses_a_south eq min(airmasses_b_north + airmasses_a_south)))) ; Julian date of lowest airmass for both, now to the hour
;   best_airmass_a_north = sigfig(airmasses_a_north(where(night_and_hour eq best_night_and_hour)), 3)
;   best_airmass_b_south = sigfig(airmasses_b_south(where(night_and_hour eq best_night_and_hour)), 3)
   best_airmass_a_south = sigfig(airmasses_a_south(where(night_and_hour eq best_night_and_hour)), 3)
   best_airmass_b_north = sigfig(airmasses_b_north(where(night_and_hour eq best_night_and_hour)), 3)
   orientation_north = -sigfig(orientations_north(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees East of North
   orientation_south = -sigfig(orientations_south(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees East of North
   zenith_north = sigfig(zeniths_north(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees
   zenith_south = sigfig(zeniths_south(where(night_and_hour eq best_night_and_hour)), 3) ; in degrees 
   best_airmass_north = best_airmass_b_north ; simply recording B as the best for north
   best_airmass_south = best_airmass_a_south ; simply recording A as the best for south  
   if print_long_info then begin
    print, ' S-N ', best_airmass_south, ' ', best_airmass_north, separation(k), best_month, best_day,  n_elements(best_nights_ba)
    print, ' IDs ', name_a, star_mag_a, round(offset_a), ' ', name_b, star_mag_b, round(offset_b)
   endif
  endif

  ; and swappable, that is, visible by either telescope
  if experiment eq 'same_field' then begin ; no need to check, if these targets are known to be in the same field
   swappable = 'yes' ; each target is swappable with the other
   swapped = 'all' ; all can be swapped 
  endif else begin
   swappable = 'no' ; by default, the targets are not swappable, but now that check is made
   swappable_nights = night(where(airmass_a_north lt airmass_limit and airmass_b_north lt airmass_limit and airmass_a_south lt airmass_limit and airmass_b_south lt airmass_limit)) ; nights when both targets are visible from either telescope at sunset
;   swappable_nights = best_night_and_hour(where(best_airmass_a_north lt airmass_limit and best_airmass_b_north lt airmass_limit and best_airmass_a_south lt airmass_limit and best_airmass_b_south lt airmass_limit)) ; nights when both targets are at some point best visible from either telescope, doesn't seem to work
   swaps = n_elements(swappable_nights) ; the number of nights swappable, which may return -1 or zero, so in either case none
   if swaps ge 2 then begin ; if there are at least two real nights that these can be viewed from either site, fails for a single night as a swaps=1 case can also indicate a null result
    swapped = swapped + 1
    swappable = 'yes'
;    swappable_nights = swappable_nights-minimum_nights ; count back the minimum required number of nights
    caldat, swappable_nights, swappable_months, swappable_days, year, hour, minute, second ; convert each swappable Julian night to month and day format
    print, 'Swappable'
;    print, 'Swappable for month and day: '
;    for swap = 0, swaps - 1 do begin
;     print, swappable_months(swap), swappable_days(swap)
;    endfor
   endif
  endelse
  
  ; and, for Gemini on that night for those two best airmasses can, at best satisfy the following condition
  noise_north_best = zeta*best_airmass_north/sqrt(float(samples))
  noise_south_best = zeta*best_airmass_south/sqrt(float(samples))
  ratio_causal(k) = 2.*signal_causal(min(where(phi ge separation(k))))/sqrt(noise_north_best^2. + noise_south_best^2.)
  ratio_direct(k) = 2.*signal(min(where(phi ge separation(k))))/sqrt(noise_north_best^2. + noise_south_best^2.) ; calculated directly from the two best airmass results
;  ratio_best(k) = ratio_twin(min(where(phi ge separation(k)))) ; instead, taken to be that for the input model value at that separation
  plotsym, 0, 2., /fill, color=0 ; generic case, black
  if plotting eq 'ratio' and experiment eq 'generic' then begin
;   ; use the generic symbol, but change the colour to match the case
;   if experiment ne 'generic' then loadct, 39, /silent ; colour, angle limit
;   if experiment eq 'peak' then plotsym, 0, 2., /fill, color=150 ; peak, green
;   if experiment eq 'truly_acausal' then plotsym, 0, 2., /fill, color=220 ; truly acausal, red
;   if experiment eq 'same_field' then plotsym, 0, 2., /fill, color=50 ; same field, blue  
   wset, 13 ; plotting S/N
   if ratio_causal(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_causal(k), ratio_causal(k)], psym=8, color=0
   plotsym, 0, 2., /fill, color=100 ; dark grey
   if ratio_direct(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_direct(k), ratio_direct(k)], psym=8, color=0
;   loadct, 0, /silent ; back to greyscale
;   plotsym, 0, 2., /fill, color=0 ; back to generic case, black
  endif
  if plotting eq 'sky' then begin ; plotting targets on the whole sky
   ; use the generic symbol, but change the colour to match the case
   if experiment ne 'generic' then loadct, 39, /silent ; colour, angle limit
   if experiment eq 'peak' then plotsym, 0, 2., /fill, color=150 ; peak, green
   if experiment eq 'truly_acausal' then plotsym, 0, 2., /fill, color=220 ; truly acausal, red
   if experiment eq 'same_field' then plotsym, 0, 2., /fill, color=50 ; same field, blue  
   wset, 14 ; both objects on the same sky
   oplot, [ra_a(k), ra_a(k)], [dec_a(k), dec_a(k)], psym=8, color=0
   oplot, [ra_b(k), ra_b(k)], [dec_b(k), dec_b(k)], psym=8, color=0
;   oplot, [ra_a(k), ra_b(k)], [dec_a(k), dec_b(k)], psym=-8, color=0 ; connect those with a line between them
   loadct, 0, /silent ; back to greyscale   
  endif
  if plotting eq 'polar' then begin ; plotted with North up and East left
   ; use the generic symbol, but change the colour to match the case
   if experiment ne 'generic' then loadct, 39, /silent ; colour, angle limit
   if experiment eq 'peak' then plotsym, 0, 2., /fill, color=150 ; peak, green
   if experiment eq 'truly_acausal' then plotsym, 0, 2., /fill, color=220 ; truly acausal, red
   if experiment eq 'same_field' then plotsym, 0, 2., /fill, color=50 ; same field, blue
   wset, 15 ; North
;   oplot, [best_airmass_north-1., best_airmass_north-1.], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
   oplot, [zenith_north, zenith_north], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
   wset, 16 ; South
;   oplot, [best_airmass_south-1., best_airmass_south-1.], (!pi/180.)*[orientation_south+90., orientation_south+90.], /polar, psym=8, color=0
   oplot, [zenith_south, zenith_south], (!pi/180.)*[orientation_south+90., orientation_south+90.], /polar, psym=8, color=0  
   loadct, 0, /silent ; back to greyscale
  endif

  ; further, mark special cases
  if experiment eq 'peak' then plotsym, 4, 3., thick=2, color=200 ; 0 ; 200 ; upward triangle, grey or black
  if experiment eq 'truly_acausal' then plotsym, 5, 3., thick=2, color=0 ; 200 ; downward triangle, grey or black
  if experiment eq 'same_field' then plotsym, 8, 5., thick=2, color=0 ; box, black
  if plotting eq 'ratio' and experiment ne 'generic' then begin
   wset, 13 ; S/N
   if experiment ne 'peak' then begin
    if ratio_causal(k) gt ratio_wanted then oplot, [separation(k), separation(k)], [ratio_causal(k), ratio_causal(k)], psym=8, color=0
    if ratio_direct(k) gt ratio_wanted then oplot, [separation(k), separation(k)], [ratio_direct(k), ratio_direct(k)], psym=8, color=0
;   if ratio_causal(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_causal(k), ratio_causal(k)], psym=8, color=0
;   if ratio_direct(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_direct(k), ratio_direct(k)], psym=8, color=0
   endif
   if experiment eq 'peak' then oplot, [separation(k), separation(k)], [ratio_peak, ratio_peak], psym=8, color=0
  endif
  if plotting eq 'sky' then begin ; plotting targets on the whole sky
   wset, 14 ; both objects on the same sky
   oplot, [ra_a(k), ra_a(k)], [dec_a(k), dec_a(k)], psym=8, color=0
   oplot, [ra_b(k), ra_b(k)], [dec_b(k), dec_b(k)], psym=8, color=0
   if experiment eq 'truly_acausal' then oplot, [ra_a(k), ra_b(k)], [dec_a(k), dec_b(k)], psym=-8, color=0 ; connect those with a line between them
  endif
  if plotting eq 'polar' then begin ; plotted with North up and East left
   wset, 15 ; North
;   oplot, [best_airmass_north-1., best_airmass_north-1.], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
   oplot, [zenith_north, zenith_north], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
   wset, 16 ; South
;   oplot, [best_airmass_south-1., best_airmass_south-1.], (!pi/180.)*[orientation_south+90. , orientation_south+90.], /polar, psym=8, color=0
   oplot, [zenith_south, zenith_south], (!pi/180.)*[orientation_south+90., orientation_south+90.], /polar, psym=8, color=0  
  endif

  ; and any swappable ones are marked with a grey circle
  if swappable eq 'yes' then begin
   plotsym, 0, 3., thick=2, color=200
   if plotting eq 'ratio' then begin
    wset, 13 ; S/N
    if ratio_causal(k) gt ratio_wanted then oplot, [separation(k), separation(k)], [ratio_causal(k), ratio_causal(k)], psym=8, color=0
    if ratio_direct(k) gt ratio_wanted then oplot, [separation(k), separation(k)], [ratio_direct(k), ratio_direct(k)], psym=8, color=0
;    if ratio_causal(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_causal(k), ratio_causal(k)], psym=8, color=100
;    if ratio_direct(k) gt 1. then oplot, [separation(k), separation(k)], [ratio_direct(k), ratio_direct(k)], psym=8, color=0
   endif
   if plotting eq 'sky' then begin ; plotting targets on the whole sky
    wset, 14 ; both objects on the same sky
    oplot, [ra_a(k), ra_a(k)], [dec_a(k), dec_a(k)], psym=8, color=0
    oplot, [ra_b(k), ra_b(k)], [dec_b(k), dec_b(k)], psym=8, color=0
;    oplot, [ra_a(k), ra_b(k)], [dec_a(k), dec_b(k)], psym=-8, color=0 ; connect those with a line between them   
   endif   
   if plotting eq 'polar' then begin ; plotted with North up and East left
    wset, 15 ; North
;    oplot, [best_airmass_north-1., best_airmass_north-1.], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
    oplot, [zenith_north, zenith_north], (!pi/180.)*[orientation_north+90., orientation_north+90.], /polar, psym=8, color=0
    ; save a screen capture
    screen_output=tvrd(true=1)
    write_jpeg, 'figure_polar_north.jpg', screen_output, quality=100, true=1
    wset, 16 ; South
;    oplot, [best_airmass_south-1., best_airmass_south-1.], (!pi/180.)*[orientation_south+90., orientation_south+90.], /polar, psym=8, color=0
    oplot, [zenith_south, zenith_south], (!pi/180.)*[orientation_south+90., orientation_south+90.], /polar, psym=8, color=0 
    ; save a screen capture
    screen_output=tvrd(true=1)
    write_jpeg, 'figure_polar_south.jpg', screen_output, quality=100, true=1
   endif
  endif
 
 endif ; looped on whether calibration stars are available for both
 
endfor
pairs = pair ; total good pairs
print, 'There were ', pairs, ' pairs found in ', semester, ' semester, and ', swapped ,' could be swapped.'
print, 'Best observed as:
print, 'AB, or North-South: ', best_ab
print, 'BA, or South-North: ', best_ba

; plotting up the test data
if experiment eq 'test' then begin
 print, 'Also plotting up the test data.'
 plotsym, 3, 3., thick=2, color=0 ; a star
 for test = 0, tests-1 do begin
  if plotting eq 'ratio' then begin
   wset, 13 ; S/N
   oplot, [separation_test(test), separation_test(test)], [ratio_test(test), ratio_test(test)], psym=8, color=0
   plotsym, 0, 1., /fill, color=0 ; and a dot
   oplot, [separation_test(test), separation_test(test)], [ratio_test(test), ratio_test(test)], psym=8, color=0
   plotsym, 3, 3., thick=2, color=0 ; back to a star
;   xyouts, separation_test(test)+12., ratio_test(test), group_label(test), alignment=0.5, charsize=charsize, charthick=charthick, color=0
  endif
  if plotting eq 'sky' then begin ; plotting targets on the whole sky
   wset, 14 ; both objects on the same sky
   oplot, [ra_test_a(test), ra_test_a(test)], [dec_test_a(test), dec_test_a(test)], psym=8, color=0
;   xyouts, ra_test_a(test)-20., dec_test_a(test)+1., group_label(test), alignment=0.5, charsize=charsize, charthick=charthick, color=0
   oplot, [ra_test_b(test), ra_test_b(test)], [dec_test_b(test), dec_test_b(test)], psym=8, color=0
;   xyouts, ra_test_b(test)+20., dec_test_b(test)-5., group_label(test), alignment=0.5, charsize=charsize, charthick=charthick, color=0
   plotsym, 0, 1., /fill, color=0 ; and a dot, to mark South
   oplot, [ra_test_b(test), ra_test_b(test)], [dec_test_b(test), dec_test_b(test)], psym=8, color=0 
   plotsym, 3, 3., thick=2, color=0 ; back to a star
  endif
  if plotting eq 'polar' then begin ; plotted with North up and East left
   wset, 15 ; North
;   oplot, [airmass_test_a(test)-1., airmass_test_a(test)-1.], (!pi/180.)*[orientation_test_a(test)+90., orientation_test_a(test)+90.], /polar, psym=8, color=0
   oplot, [zenith_test_a(test), zenith_test_a(test)], (!pi/180.)*[orientation_test_a(test)+90., orientation_test_a(test)+90.], /polar, psym=8, color=0
;   xyouts, zenith_test_a(test), (!pi/180.)*orientation_test_a(test)+90., group(test), charsize=charsize, charthick=charthick, color=0, /data ;, /polar ; does not work in polar coordinates
   ; save a screen capture
   screen_output=tvrd(true=1)
   write_jpeg, 'figure_polar_north.jpg', screen_output, quality=100, true=1
   wset, 16 ; South
;   oplot, [airmass_test_b(test)-1., airmass_test_b(test)-1.], (!pi/180.)*[orientation_test_b(test)+90., orientation_test_b(test)+90.], /polar, psym=8, color=0
   oplot, [zenith_test_b(test), zenith_test_b(test)], (!pi/180.)*[orientation_test_b(test)+90., orientation_test_b(test)+90.], /polar, psym=8, color=0
;   xyouts, zenith_test_b(test), (!pi/180.)*orientation_test_b(test)+90., group(test), charsize=charsize, charthick=charthick, color=0, /data ;, /polar ; does not work in polar coordinates
   plotsym, 0, 1., /fill, color=0 ; and a dot, to mark South
;   oplot, [airmass_test_b(test)-1., airmass_test_b(test)-1.], (!pi/180.)*[orientation_test_b(test)+90., orientation_test_b(test)+90.], /polar, psym=8, color=0
   oplot, [zenith_test_b(test), zenith_test_b(test)], (!pi/180.)*[orientation_test_b(test)+90., orientation_test_b(test)+90.], /polar, psym=8, color=0
   plotsym, 3, 3., thick=2, color=0 ; back to a star
   ; save a screen capture
   screen_output=tvrd(true=1)
   write_jpeg, 'figure_polar_south.jpg', screen_output, quality=100, true=1
  endif
 endfor
endif

endfor ; going through all experiments

; save a screen capture
screen_output=tvrd(true=1)
if plotting eq 'ratio' then write_jpeg, 'figure_ratio.jpg', screen_output, quality=100, true=1
if plotting eq 'sky' then write_jpeg, 'figure_sky.jpg', screen_output, quality=100, true=1
;if plotting eq 'polar' then write_jpeg, 'figure_polar.jpg', screen_output, quality=100, true=1

endfor ; going through all plots
print, 'Done.'

end

; quantum correlation
function quantum_correlation, theta
; Returns the "classical" two-point quantum correlation of point sources. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*(-1.*cos(!pi*theta/180.)) ; test of mirroring
return, -1.*cos(!pi*theta/180.)
end

; Bell theorem
function bell_theorem, theta
; Returns the two-point quantum correlation of point sources, as in Bell's Theorem. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.) ; test of mirroring
return, abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.)
end

function airmass_simple,sza
;+
; ROUTINE:  airmass
;
; PURPOSE:  compute airmass as a function of angle, sza, including
;           spherical earth effects
;
; USEAGE:   result=airmass(sza)
;
; INPUT:    
;   sza     zenith angle   
;
; KEYWORD INPUT:
;
; OUTPUT:   relative airmass,
;           by definition the relative airmass = 1 for sza=0.
;
; References:
;   Kasten, F 1966: A new table and approximate formula for relative
;   airmass. Arch. Meteor. Geophys. Bioklimatol. Ser. B, 14, 206-223
;
;   Leontieva, E.N., and K.H. Stamnes 1996: Remote sensing of cloud
;   optical properties from fround-based measurements of
;   transmittance: a feasibility study, Journal of Applied
;   Meteorology, 35, 2011-2022
;  
; EXAMPLE:  
;
;    sza=findrng(0,89.,dx=.1)
;    plot,sza,1./cos(sza*!dtor),yran=[1,30]
;    plot,sza,cos(sza*!dtor)*airmass(sza)
;
; AUTHOR:   Paul Ricchiazzi                        09 Jan 97
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;           paul@icess.ucsb.edu
;
;-
;
;return,1./(cos((sza<90.)*!dtor)+0.15*exp(-1.253*alog(93.885-(sza<90.))))
return, 1./cos((sza<90.)*!dtor) ; test, setting to the "extreme" case, that is, a simple inverse cosine

end

;+
; NAME:
;        SIGFIG
;
;
; PURPOSE:
;        Accept a scalar numerical value or an array of numbers and
;        return the numbers as strings with the specified number of
;        significant figures.
;
; CALLING SEQUENCE:
;        RESULT = SigFig(Number, Nfig [, /SCIENTIFIC, /PLUSSES, /NUMERICAL)
;
; INPUTS:
;        Number - Scalar or array of numerical values (float, double, int)
;        Nfig   - Number of signficant figures desired in the output
;
; OUTPUTS:
;        String representation of the input with the specified number
;        of signficant figures.
;
; KEYWORD PARAMTERS:
;        /SCIENTIFIC - return the numbers in scientific notation
;        /PLUSSES    - Include plus signs for positive numbers 
;        /NUMERICAL  - Return numerical, rather than string, values
;
; EXAMPLE:
;        IDL> print, sigfig(-0.0001234, 2)      
;        -0.00012
;        IDL> print, sigfig(1.234, 1)
;        1.
;        IDL> print, sigfig(1234, 1) 
;        1000
;        IDL> print, sigfig(-0.0001234, 2, /sci)
;        -1.2e-4
;        IDL> print, sigfig(1234, 2, /plus)
;        +1200
;        IDL> print, sigfig(1234, 2, /plus, /sci)
;        +1.2e+3
;
; MODIFICATION HISTORY:
; Inspired long ago by Erik Rosolowsky's SIGFIG:
;     http://www.cfa.harvard.edu/~erosolow/idl/lib/lib.html#SIGFIG
;
; This version written by JohnJohn Sept 29, 2005
;
; 24 Oct 2007 - If result is a single number, return scalar value
;               instead of an 1-element array. Thanks Mike Liu.
;  2 Apr 2008 - Fixed 1-element array issue, but for real this time.
;-

;;; SF_STR - The way STRING() should behave by default
function sf_str, stringin, format=format
return, strcompress(string(stringin, format=format), /rem)
end

;;; SF_TRANS_DEC - TRANSlate the DECimal point in a number of order
;;;                unity, round it, and translate back. 
function sf_trans_dec, numin, nsigin, order_inc=order_inc
nel = n_elements(numin)

;;; Double precision can't handle more than 19 sig figs
nsig = nsigin < 19

;;; Gonna have to move the decimal nsig-1 places to the right before rounding
move = nsig-1
len = max(strlen(numin))
move = move < (len-1)

;;; Create a string with just the digits, no decimal
nodec = strmid(numin, 0, 1)+strmid(numin, 2, len)

;;; Move the decimal, so nsig digits are to the left of the new
;;; decimal position
num0 = strmid(nodec,0,1+move)+'.'+strmid(nodec,1+move,len)

;;; Round the new number
num1 = sf_str(round(double(num0),/l64))
len1 = strlen(num1)

;;; If the number increases an order of magnitude after rounding, set
;;; order_inc=1 so the calling routine knows to add one to the order 
;;; of magnitude
order_inc = len1 gt nsig
;;; Move the decimal back and return to sender
num  = strmid(num1, 0, 1)+'.'+strmid(num1, 1, nsig-1)
return, num
end

function sigfig, NumIn, Nfig $
                 , string_return=string_return $
                 , scientific=scientific $
                 , numerical=numerical $
                 , plusses=plusses

Num = double(NumIn)
Nel = n_elements(Num)

;;; Convert the input number to scientific notation
TestString = sf_str(abs(double(Num)), format='(e)')
Epos = strpos(TestString[0], 'e')

;;; Test sign of the order
Osign = intarr(Nel)+1
StrOsign = strmid(TestString, Epos+1, 1)
Wneg = where(strosign eq '-', Nneg) 
if Nneg gt 0 then Osign[Wneg] = -1

;;; Test sign of numbers, form string of minus signs for negative vals
NegSign = strarr(Nel) + (keyword_set(plusses) ? '+' : '')
Negative = where(Num lt 0, Nneg)
if Nneg gt 0 then NegSign[Negative] = '-'

;;; What are the orders of magnitude of the values?
Order = fix(sf_str(strmid(TestString, Epos+2, 2)))

;;; Convert all values to order unity for rounding
NumUnit = strmid(TestString,0,epos)

;;; Use TRANS_DEC to round unit values
NumTrans = sf_trans_dec(NumUnit, Nfig, order_inc=Order_Inc)
Order = order + Osign*order_inc
Len = strlen(NumTrans[0])

;;; Exit early without looping for /NUMERICAL or /SCIENTIFIC
if keyword_set(numerical) then begin
    NumRound = NegSign+NumTrans+'e'+StrOsign+sf_str(Order)
    if n_elements(NumRound) eq 1 then return, double(NumRound[0]) else $
      return, double(NumRound)
endif
if keyword_set(scientific) then begin
    NumRound = NegSign+NumTrans+'e'+StrOsign+sf_str(Order)
    if n_elements(NumRound) eq 1 then return, NumRound[0] else $
      return, NumRound
endif

NumRound = strarr(Nel)
for i = 0, Nel-1 do begin
    if Osign[i]*Order[i]+1 gt Nfig then Format = '(I40)' else begin
        d = sf_str(fix(Nfig-(Osign[i]*Order[i])-1) > 0)
        Format = '(F40.' + d + ')'
    endelse
    New = NumTrans[i] * 10d^(Osign[i] * Order[i])
    NumRound[i] = NegSign[i]+sf_str(New, format=Format)
endfor
if n_elements(NumRound) eq 1 then return, NumRound[0]
return, NumRound
end

PRO gcirc,u,ra1,dc1,ra2,dc2,dis                         
;+
; NAME:
;     GCIRC
; PURPOSE:
;     Computes rigorous great circle arc distances.  
; EXPLANATION:
;     Input position can either be either radians, sexagesimal RA, Dec or
;     degrees.   All computations are double precision. 
;
; CALLING SEQUENCE:
;      GCIRC, U, RA1, DC1, RA2, DC2, DIS
;
; INPUTS:
;      U    -- integer = 0,1, or 2: Describes units of inputs and output:
;              0:  everything radians
;              1:  RAx in decimal hours, DCx in decimal
;                       degrees, DIS in arc seconds 
;              2:  RAx and DCx in degrees, DIS in arc seconds
;      RA1  -- Right ascension or longitude of point 1
;      DC1  -- Declination or latitude of point 1
;      RA2  -- Right ascension or longitude of point 2
;      DC2  -- Declination or latitude of point 2
;
; OUTPUTS:
;      DIS  -- Angular distance on the sky between points 1 and 2
;              See U above for units;  double precision  
;
; PROCEDURE:
;      "Haversine formula" see 
;      http://en.wikipedia.org/wiki/Great-circle_distance
;
; NOTES:
;       (1) If RA1,DC1 are scalars, and RA2,DC2 are vectors, then DIS is a
;       vector giving the distance of each element of RA2,DC2 to RA1,DC1.
;       Similarly, if RA1,DC1 are vectors, and RA2, DC2 are scalars, then DIS
;       is a vector giving the distance of each element of RA1, DC1 to 
;       RA2, DC2.    If both RA1,DC1 and RA2,DC2 are vectors then DIS is a
;       vector giving the distance of each element of RA1,DC1 to the 
;       corresponding element of RA2,DC2.    If the input vectors are not the 
;       same length, then excess elements of the longer ones will be ignored.
;
;       (2) The function SPHDIST provides an alternate method of computing
;        a spherical distance.
;
;       (3) The haversine formula can give rounding errors for antipodal
;       points.
;
; PROCEDURE CALLS:
;      None
;
;   MODIFICATION HISTORY:
;      Written in Fortran by R. Hill -- SASC Technologies -- January 3, 1986
;      Translated from FORTRAN to IDL, RSH, STX, 2/6/87
;      Vector arguments allowed    W. Landsman    April 1989
;      Prints result if last argument not given.  RSH, RSTX, 3 Apr. 1998
;      Remove ISARRAY(), V5.1 version        W. Landsman   August 2000
;      Added option U=2                      W. Landsman   October 2006
;      Use double precision for U=0 as advertised R. McMahon/W.L.  April 2007
;      Use havesine formula, which has less roundoff error in the 
;             milliarcsecond regime      W.L. Mar 2009
;-
 compile_opt idl2
 On_error,2                            ;Return to caller

 npar = N_params()
 IF (npar ne 6) and (npar ne 5) THEN BEGIN
   print,'Calling sequence:  GCIRC,U,RA1,DC1,RA2,DC2[,DIS]'
   print,'   U = 0  ==> Everything in radians'
   print, $
   '   U = 1  ==> RAx decimal hours, DCx decimal degrees, DIS arc sec'
   print,'   U = 2  ==> RAx, DCx decimal degrees, DIS arc sec'
   RETURN
 ENDIF


 d2r    = !DPI/180.0d0
 as2r   = !DPI/648000.0d0
 h2r    = !DPI/12.0d0

; Convert input to double precision radians
 CASE u OF
   0:  BEGIN
          rarad1 = double(ra1)
          rarad2 = double(ra2)
          dcrad1 = double(dc1)
          dcrad2 = double(dc2)
       END
   1:  BEGIN
          rarad1 = ra1*h2r
          rarad2 = ra2*h2r
          dcrad1 = dc1*d2r
          dcrad2 = dc2*d2r
       END
   2:  BEGIN  
          rarad1 = ra1*d2r
          rarad2 = ra2*d2r
          dcrad1 = dc1*d2r
          dcrad2 = dc2*d2r
        END
   ELSE:  MESSAGE, $
                'U must be 0 (radians), 1 ( hours, degrees) or 2 (degrees)'
 ENDCASE

 deldec2 = (dcrad2-dcrad1)/2.0d
 delra2 =  (rarad2-rarad1)/2.0d
 sindis = sqrt( sin(deldec2)*sin(deldec2) + $
	  cos(dcrad1)*cos(dcrad2)*sin(delra2)*sin(delra2) )
 dis = 2.0d*asin(sindis) 

 IF (u ne 0) THEN dis = dis/as2r

 IF (npar eq 5) && (N_elements(dis) EQ 1) THEN BEGIN
    IF (u ne 0) && (dis ge 0.1) && (dis le 1000)  $
       THEN fmt = '(F10.4)' $
       ELSE fmt = '(E15.8)'
    IF (u ne 0) THEN units = ' arcsec' ELSE units = ' radians'
    print,'Angular separation is ' + string(dis,format=fmt) + units
 ENDIF

 RETURN
 END




