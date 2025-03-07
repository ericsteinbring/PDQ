pro simulate_quasar_correlation

; Simulates a test of detecting correlated quasar photons, and reports the signal-to-noise ratio for relative photometry of two sources with equal magnitudes, with inputs for zeropoint of filters, efficiency of detectors, size of telescope and extinction.  The number of samples is an input to the final, combined or stacked-exposures result, allowing it to report a reasonable total integration time to achieve this, assuming Poisson statistics.
; Eric Steinbring, 17 October 2024.

; parameters
diameter       = 8.1 ; 1.12838 ; 4. ; 8.1 ; 30 ; m, of the telescope aperture
obstruction    = 1. ; 0.25 ; 1. ; 2. ; m, diameter of central obstruction
lambda_c       = 0.55 ; 0.64 ; 0.55 ; microns, central wavelength of filter
delta_lambda   = 0.16 ; 0.23 ; 0.16 ; microns, bandwidth of filter
flux_0         = 3631. ; 3080. ; 3640. ; 3631. ; Jy ; zeropoint in flux units, where 3631 is set by the definition of the AB-magnitude zeropoint
exposure       = 0.06 ; 0.001 ; 0.005 ; 0.01 ; 0.015 ; 0.02 ; 0.03 ; 0.06 ;  0.1 ; 1. ; exposure time in seconds
overhead       = 0.25 ; 5. ; 2. ; 1. ; 0.5 ; 0.25 ;  0. ; factor of exposure time in addition, as overhead per exposure, unity means that the overhead is equal to the exposure time
efficiency     = 0.85 ; 0.5 ; 0.8 ; 0.85 ; 0.9 ; 1. ; fractional efficiency of detector
magnitude      = 18.5 ; 15. ; 15.5 ; 16. ; 16.5 ; 17. ; 17.5 ; 18. ; 18.5 ; 19. ; 19.5 ; 20. ; mag, of sources; default is 18.5 for roughly a redshift 4 quasar
psf_fwhm       = 0.8 ; 0.5 ; 0.65 ; 0.8 ; 1. ; arcsec, FWHM of image, assumed stellar
extinction     = 0.5 ; 0. ; 0.1 ; 0.2 ; 0.25 ; 0.5 ; 1. ; 2. ; mag, integrated over filter bandpass
sky            = 19. ; 15. ; 16. ; 17. ; 18. ; 19. ; 20. ; 21. ; 22. ; mag arcsec^-2, sky brightness integrated over filter bandpass
samples        = 1000 ; 1.e6 ; 0.5e6 ; 1.e5 ; 10000 ; 1000 ; 500 ; 100 ; 10 ; samples of each exposure, default is 1000
time_minimum   = 0.5 ; 1. ; 2. ; 3. ; 4. ; 5. ; seconds
shifts_maximum = 100 ; 100 ; 250 ; 500 ; 1000 ; per cube, maximum allowed shift in time, looking for variation in the overall correlation of frames
see_signals    = 'yes' ; 'yes' ; 'no' ; 'flat' ; whether there are exploitable signals, washed out by Gaussian or white noise (flat); default is minimal, exploitable-case signal, "flat" is default, and now defunct
add_noise      = 'yes' ; 'yes' ; 'no' ; whether test is noisy or not; default is to add noise

; calculation
shifts_minimum = time_minimum/exposure
print, 'Minumum shifts: ', shifts_minimum
area = !pi*((diameter-obstruction)/2.)^2. ; m^2
source = flux_0*10.^(-0.4*(magnitude+extinction))*1.51e7*delta_lambda ; photons sec^-1 m^-2
flux = area*source ; photons sec^-1
fluence = flux*exposure ; photons per exposure time
psf = !pi*(psf_fwhm/2.)^2. ; arcsec^2
signal = exposure*efficiency*flux/psf ; fluence in photons arcsec^-2
noise = exposure*flux_0*10.^(-0.4*sky*1.)*1.51e7*delta_lambda/psf ; photons arcsec^-2
zeropoint = -2.5*alog10(flux_0) + 8.90 ; magnitudes ; where flux_0 is in Jy flux units, and zeropoint defined to be zero in the AB system

; report
print, 'Diameter (m): ', diameter
print, 'Obstruction (m): ', obstruction
print, 'Area (m^2): ',  area
;print, 'Zeropoint (mag): ', zeropoint
print, 'Source magnitude: ', magnitude
print, 'Total (with extinction): ', magnitude+extinction
print, 'Source (photons s^-1 m^-2): ', source
print, 'Flux (photons s^-1): ', flux
print, 'Fluence (total photons): ', fluence
print, 'Seeing (arcsec): ', psf_fwhm
print, 'Source area (arcsec^2): ', psf
print, 'Exposure (s): ', exposure
print, 'Overhead factor: ', overhead
print, 'Efficiency (percent): ', 100.*efficiency
print, 'Signal (photons arcsec^-2): ', signal
print, 'Noise (photons arcsec^-2): ', noise
print, 'Individual S/N: ', signal/noise
print, 'Approx. photometric error (mag): ', 1./(signal/noise)
print, 'Relative pairwise S/N: ', signal/(sqrt(2.)*noise)
print, 'Samples: ', samples
print, 'Final pairwise S/N: ', sqrt(float(samples))*signal/(sqrt(2.)*noise)
print, 'Integration (s, m, h): ', float(samples)*exposure*(1.+overhead), float(samples)*exposure*(1.+overhead)/60., float(samples)*exposure*(1.+overhead)/60./60.

; calculations
timesteps = findgen(samples)+1

; generate a simulated, perfect, minimal toggle-switch, up or down at each timestep
toggle = replicate(1., samples) ; initially set positive, i.e. "throw" of switch is "up"
for i = 0, samples - 1 do begin
 if timesteps(i) mod 2 eq 1. then toggle(i) = -1. ; at each odd-numbered step, the throw is "down"
endfor
toggle = -1.*toggle ; reverse the toggle orientation, at the initial step
;toggle = 0.5*toggle ; can also make it a total throw of unity instead
;toggle = toggle - mean(toggle, /nan) ; ensure it is centred

; a perfect sinusoidal signal, akin to the quasars correlated like a coherent "siren"
siren = fltarr(samples)
siren = 1.*cos(!pi*timesteps/(time_minimum/exposure)) ; a sine-wave with period equal to the minimum time
;siren = 0.5*siren ; can also make it a total amplitude of unity instead
;siren = siren - mean(siren, /nan) ; median(siren) ; ensure it is centred

; a sequence of perfect signals, plus systematic noise, allowing relative degradation of the output signal
perfect_signals = toggle ; (toggle+siren)/2. ;toggle ; siren ; toggle ; perfect signal, that is, the toggle will always provide a spoiled-sample result
no_signals = 1.-randomn(seed, samples) ; adding an equal unit of Gaussian noise, reducing signal to zero, on average
flat_signals = 1.-randomu(seed, samples) ; adding an equal unit of white noise, reducing signal to exactly zero
; now, force the case that no signal means white noise
no_signals = flat_signals
fair_signals = (no_signals + perfect_signals)/2. ; test: now setting exploitable S/N to be exactly 2:1, that is just as likely to spoil a sample as not
exploitable_signals = (no_signals + 3.*perfect_signals)/4. ; now setting exploitable S/N to be 3, that is 3 out-of-four samples will be spoiled, which is the lower limit for spoiling a test

; setting the signal to look for
;signals = perfect_signals ; test: perfect signal
signals = no_signals ; test: no signal
;if see_signals eq 'yes' then signals = fair_signals ; test: just-fair signal
if see_signals eq 'yes' then signals = exploitable_signals ; test: just-exploitable signal
;if see_signals eq 'flat' then signals = flat_signals ; test: pure white-noise or "flat" signal

; setup display
window, 1, xsize = 1000, ysize = 500, title='Schematic of an Ideally-Correlated Signal'
wset, 1
!p.multi=[0,1,1]
charsize=2.5
loadct, 0, /silent ; greyscale, for models
;loadct, 1, /silent ; bluescale, for models
;loadct, 39, /silent ; colour, for data
!p.background=255
device, set_font='helvetica bold', /tt_font
!p.font=1
!except=0

; plot a schematic of what a "toggle-switch" signal would look like
plotsym, 0, 2., /fill, color=0 ; default, filled black circle
plot, timesteps, perfect_signals, xstyle=1, xrange=[0, 10], xtitle='Samples', xticks=10, ystyle=1, yrange=[-2., 2.], yticks=4, ytitle='Perfect signal (in std. dev. from mean)', charsize=charsize, charthick=chartick, psym=8, color=0
oplot, timesteps, toggle, thick=7, color=200
;; no exploitable signal
;plotsym, 0, 2., /fill, color=200 ; light filled grey circle
;oplot, timesteps, no_signals, psym=8, thick=5, color=0
;; the ignal to be looked for
;plotsym, 0, 2., /fill, color=100 ; dark filled grey circle
;oplot, timesteps, exploitable_signals, psym=8, thick=3, color=0
; plot the "siren" case
oplot, timesteps, siren, color=0
; replot the perfect case
plotsym, 0, 2., /fill, color=0 ; default, filled black circle
oplot, timesteps, perfect_signals, psym=8, color=0
; limits
oplot, [3, 3], [-2, 2], linestyle=1, color=0
oplot, [0, samples], [0., 0.], linestyle=3, color=0

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_schematic_toggle_switch.jpg', screen_output, quality=100, true=1
;write_gif, 'figure_schematic_toggle_switch.gif', tvrd(screen_output)

; simulate the relative correlated colour-signal, that is the relative photon flux North and South
colour_n = (signal/noise)*signals ;+ noise*randomn(seed, samples) ; adding weighted uncorrelated noise
colour_s = -(signal/noise)*signals ;+ noise*randomn(seed, samples) ; adding weighted uncorrelated noise

; adding weighted uncorrelated noise
if add_noise eq 'yes' then begin
 colour_n = colour_n + noise*randomn(seed, samples)
 colour_s = colour_s + noise*randomn(seed, samples)
endif

; each sample is then equivalent to one frame of photometry
frames = samples
 
; take the colour difference
colour = colour_n - colour_s

; report statistics
colour_mean = mean(colour)
colour_median = median(colour)
colour_stddev = stddev(colour)
print, 'Colour mean, median and std. dev.: ', colour_mean, colour_median, colour_stddev 

; and normalize 
;colour = (colour - colour_median)/colour_mean ; normalized to the mean
colour =  (colour - colour_median)/colour_stddev ; normalized to one standard deviation
 
; look for a trend
frame_number = findgen(frames)
trend = poly_fit(frame_number, colour, 1, sigma=sigma, yfit=colour_fit) ; just linear, to start
;print, 'Slope of colour-change per frame: ', trend(1)

; and subtract it from the data
colour = colour - colour_fit
 
; look for a trend, again
frame_number = findgen(frames)
trend = poly_fit(frame_number, colour, 1, sigma=sigma, yfit=colour_fit) ; just linear, to start
; print, 'Slope of colour-change per frame: ', trend(1)

; set up labeling
label_signals = '_signal' ; that is, it is the default of 'yes'
label_noise = '_noise' ; that is, it is the default of 'yes'
label_samples = '_'+strtrim(string(samples),2)
xyouts_signals = 'Exploitable signal' ; 'Exploitable correl.' ; 'Exploitable signal' ; that is, it is the default of 'yes' is for exploitable
xyouts_noise = ', noise-limited' ; that is, it is the default of 'yes'
xyouts_samples = strtrim(string(samples),2)
if see_signals eq 'no' then label_signals = '_no_signal'
;if see_signals eq 'flat' then label_signals = '_flat_signal'
if add_noise eq 'no' then label_noise = '_no_noise'
if see_signals eq 'no' then xyouts_signals = 'No correlation' ; 'Flat signal' ; 'No signal'
;if see_signals eq 'flat' then xyouts_signals = 'Flat signal'
if add_noise eq 'no' then xyouts_noise = ' ' ; 'with minimal noise' ; 'no noise'

; set up the display
window, 2, xsize = 1500, ysize = 500, xpos=0, ypos=0, title='Relative quasar colours'
;window, 2, xsize = 750, ysize = 500, xpos=0, ypos=0, title='Relative quasar colours'
wset, 2
!p.multi = [0, 1, 1]
charsize=2.5

; plot the relative colours
plotsym, 0, 2., /fill, color=0 ; black filled circle, default
plot, colour, xtitle='Sample', ytitle='Rel. A-B colour (in std. dev. from mean)', ystyle=1, yrange=[-4., 4.], charsize=charsize, charthick=charthick, psym=8, color=0
; shading
for l = 0, frames-1 do begin
 oplot, [l, l], [-1., 1.], thick=5, color=220
endfor
; plot the "siren" case
oplot, siren, color=0
; replot
oplot, colour, psym=8, color=0
; limits
;oplot, frame_number, colour_fit, thick=3, color=100
;oplot, [0, frames], [0., 0.], linestyle=2, color=0
oplot, [0, frames], [0., 0.], linestyle=3, color=0
; labels
xyouts, 25, 3.25, xyouts_signals + xyouts_noise, charsize=charsize, charthick=charthick, color=0
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[0, frames], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[0, frames], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[-4., 4.], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[-4., 4.], ytickformat='(A1)'
 
; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_simulated_relative_colours'+label_signals+label_noise+label_samples+'.jpg', screen_output, quality=100, true=1
;write_gif, 'figure_simulated_relative_colours'+label_signals+label_noise+label_samples+'.gif', tvrd(screen_output)

; setup display
window, 3, xsize = 600, ysize = 500, xpos=1400, ypos=0, title='Distribution of an Ideally-Correlated Signal'
wset, 3
!p.multi=[0,1,1]

; find the distribution of colours, and signal
binsize=0.25 ; 0.25 ; 0.1 ; 0.2 ; 0.25 ; 0.5
distribution_min = -4.
distribution_max = 4.
distribution_siren = histogram(siren, min=distribution_min, max=distribution_max, binsize=binsize)
distribution_toggle = histogram(toggle, min=distribution_min, max=distribution_max, binsize=binsize)
distribution_colour = histogram(colour, min=distribution_min, max=distribution_max, binsize=binsize)
distribution_noise = histogram(noise*randomn(seed, samples), min=distribution_min, max=distribution_max, binsize=binsize)
bins = binsize*findgen(n_elements(distribution_colour)) + distribution_min

; smooth
;distribution_siren = smooth(distribution_siren, 3)
;distribution_toggle = smooth(distribution_toggle, 3)
;distribution_colour = smooth(distribution_colour, 3)
;distribution_noise = smooth(distribution_noise, 3)

; normalize
;max_distribution = float(max([distribution_siren, distribution_toggle, distribution_colour, distribution_noise]))
total_distribution = float(total(distribution_toggle))
distribution_siren = distribution_siren/total_distribution ;/max_distribution
distribution_toggle = distribution_toggle/total_distribution ;max_distribution
distribution_colour = distribution_colour/total_distribution ;max_distribution
distribution_noise = distribution_noise/total_distribution ;max_distribution
max_distribution = float(max(distribution_colour))
distribution_normalized = distribution_noise/max_distribution ; set to the width of the noise, and peak of the colour distribution

; report that of a pure, simple normal distribution
normal = gaussian(bins, [max_distribution, 0., 1.]) ; set to the peak of the colour distribution
;normal = gaussian(bins, [0.5, 0., 1.]) ; set to be a peak of 50%
;normal = gaussian(bins, [0.25, 0., 1.]) ; set to be a peak of 25%

; set these to be simply the same thing
distribution_normalized = normal

; and plot them
plot, distribution_siren, bins+binsize/2., xstyle=1, xtitle='Normalized occurances', xrange=[0.0, 0.5], xticks=5, ystyle=1, ytitle='Distribution of A-B colour', yrange=[distribution_min, distribution_max], charsize=charsize, charthick=charthick, color=0 ;, ytickformat='(A1)'
; shading
for i = 0, n_elements(bins)-1 do begin
 if bins(i) gt -1. and bins(i) lt 1. then oplot, [0., 0.5], [bins(i), bins(i)], thick=20, color=220
; oplot, [0., normal(i)], [bins(i), bins(i)], thick=20, color=220
endfor
; replot
oplot, distribution_toggle, bins, thick=2, color=200
oplot, distribution_siren, bins+binsize/2., color=0
;oplot, distribution_noise, bins, thick=2, color=150
oplot, distribution_colour, bins+binsize/2., thick=3, color=0
; overplot a Gaussian, or normalized distribution
;oplot, normal, bins, thick=3, color=100
oplot, distribution_normalized, bins, thick=3, linestyle=2, color=100
; limits
oplot, [0., 0.5], [0., 0.], linestyle=3, color=0
; labels
xyouts, 0.025, 3.25, xyouts_signals + xyouts_noise, charsize=charsize, charthick=charthick, color=0
if add_noise eq 'yes' then xyouts, 0.025, 2.75, 'Samples: ' + xyouts_samples, charsize=charsize, charthick=charthick, color=0
; clean up
;axis, xaxis=0, color=0, xstyle=1, xrange=[0., 0.5], charsize=charsize
;axis, xaxis=1, color=0, xstyle=1, xrange=[0., 0.5], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[distribution_min, distribution_max], charsize=charsize ;ytickformat='(A1)'
axis, yaxis=1, color=0, ystyle=1, yrange=[distribution_min, distribution_max], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_simulated_distribution_colours'+label_signals+label_noise+label_samples+'.jpg', screen_output, quality=100, true=1
;write_gif, 'figure_simulated_distribution_colours'+label_signals+label_noise+label_samples+'.gif', tvrd(screen_output)

; now shift the relative photometry by time, up to the prescribed limit, and look for a maximum and minimum of correlation
shifts = findgen(shifts_maximum) ; start at none and shift South relative to North
;shifts = -floor(shifts_maximum/2.) + findgen(shifts_maximum) ; start by goning back half the maximum number of shifts
colours = fltarr(n_elements(shifts))
correlations = colours
for k = 0, n_elements(shifts)-1 do begin

 ; shift the photometry
; colour_n = shift(colour_n, shifts(k)) ; shift by number of exposures
 colour_s = shift(colour_s, shifts(k)) ; shift by number of exposures
 
 ; take the colour difference
 colour = colour_n - colour_s
 
; ; report statistics
; colour_mean = mean(colour)
; colour_median = median(colour)
; colour_stddev = stddev(colour)
;; print, 'Sequence and shift rel. colour mean, median and std. dev.: ', j+1, k+1, colour_mean, colour_median, colour_stddev 

 ; normalize 
; colour = (colour - colour_median)/colour_mean ; normalized to the mean
 colour = (colour - colour_median)/colour_stddev ; normalized to one standard deviation
 
 ; and record the degree of correlation, e.g., recording the colour standard deviation for the given shift
 colours(k) = mean(abs(colour))
 correlations(k) = stddev(colour)/2.
 
endfor
 
; take absolutes or normalize
;colours = abs(colours)
;colours = colours-mean(colours)
;colours = abs(colours)
;colours = colours/mean(colours)
;correlations = abs(correlations)
;correlations = correlations-mean(correlations)
;correlations = abs(correlations)
;correlations = correlations/mean(correlations)

; a smoothed fit
;colours_smooth = smooth(colours, floor(shifts_maximum/10.), /edge_wrap)
;correlations_smooth = smooth(correlations, floor(shifts_maximum/10.), /edge_wrap)
colours_smooth = smooth(colours, shifts_minimum, /edge_wrap)
correlations_smooth = smooth(correlations, shifts_minimum, /edge_wrap)
 
; find the means
colours_mean = mean(colours_smooth, /nan)
correlations_mean = mean(correlations_smooth, /nan)
 
;; find the maximum and minimum, overall; if multiple, take the first one
;colours_max = max(colours)
;shifts_colours_max = shifts(min(where(colours eq max(colours))))
;correlations_max = max(correlations)
;shifts_correlations_max = shifts(min(where(correlations eq max(correlations))))
;colours_min = min(colours)
;shifts_colours_min = shifts(min(where(colours eq min(colours))))
;correlations_min = min(correlations)
;shifts_correlations_min = shifts(min(where(correlations eq min(correlations))))

;; find the maximum and minimum, after smoothing; if multiple, take the first one
;colours_max = max(colours_smooth)
;shifts_colours_max = shifts(min(where(colours_smooth eq max(colours_smooth))))
;correlations_max = max(correlations_smooth)
;shifts_correlations_max = shifts(min(where(correlations_smooth eq max(correlations_smooth))))
;colours_min = min(colours_smooth)
;shifts_colours_min = shifts(min(where(colours_smooth eq min(colours_smooth))))
;correlations_min = min(correlations_smooth)
;shifts_correlations_min = shifts(min(where(correlations_smooth eq min(correlations_smooth))))
 
; and further, after smoothing find the maximum and minimum within the minimum time; if multiple, take the first one
colours_smooth_short = colours_smooth
colours_smooth_short(where(shifts gt shifts_minimum)) = !values.f_nan
colours_max = max(colours_smooth_short)
shifts_colours_max = shifts(min(where(colours_smooth_short eq max(colours_smooth_short, /nan))))
correlations_smooth_short = correlations_smooth
correlations_smooth_short(where(shifts gt shifts_minimum)) = !values.f_nan
correlations_max = max(correlations_smooth_short, /nan)
shifts_correlations_max = shifts(min(where(correlations_smooth_short eq max(correlations_smooth_short, /nan))))
colours_min = min(colours_smooth_short, /nan)
shifts_colours_min = shifts(min(where(colours_smooth_short eq min(colours_smooth_short, /nan))))
correlations_min = min(correlations_smooth_short, /nan)
shifts_correlations_min = shifts(min(where(correlations_smooth_short eq min(correlations_smooth_short, /nan))))
  
; set up the display
window, 4, xsize = 1000, ysize = 1000, xpos = 900, ypos = 0, title='Simulated correlation of quasar colours'
wset, 4
!p.multi = [0, 1, 2]
charsize=2.5
 
; plot the colour differences
plotsym, 0, 0.5, /fill, color=0 ; default, small black dot
plot, shifts, colours, xstyle=1, xrange=[min(shifts), max(shifts)], xtickformat='(A1)', ystyle=1, ytitle='Abs. colour diff. (in std. dev. rel. to mean)', yrange=[0.95*min(colours), 1.05*max(colours)], ymargin=[0,2], charsize=charsize, charthick=charthick, psym=8, color=0
;; shading
;for l = 0, shifts_minimum-1 do begin
; oplot, [l, l], [0.95*min(colours), 1.05*max(colours)], thick=10, color=220
;endfor
; plot the "siren" case
oplot, timesteps, (colours_max-colours_mean)*siren + colours_mean, color=0
; replot
oplot, shifts, colours, psym=8, color=0
; overplotting the smoothed curve
oplot, shifts, colours_smooth, thick=3, color=100
; maxima and minima
; plotsym, 0, 2., /fill, color=0 ; black filled circle, max
plotsym, 4, 3., thick=2, color=0 ; up-pointing triangle
oplot, [shifts_colours_max, shifts_colours_max], [colours_max, colours_max], psym=8, color=0
; plotsym, 0, 2., /fill, color=200 ; grey filled circle, min
plotsym, 5, 3., thick=2, color=0 ; down-pointing triangle 
oplot, [shifts_colours_min, shifts_colours_min], [colours_min, colours_min], psym=8, color=0
; limits
oplot, [shifts_minimum, shifts_minimum], [0.95*min(colours), 1.05*max(colours)], linestyle=1, color=0
oplot, [0, shifts_maximum], [colours_mean, colours_mean], linestyle=2, color=0
oplot, [0, shifts_maximum], [0.75, 0.75], linestyle=3, color=0 ; expected for no correlation
; clean up
; axis, xaxis=0, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[0.95*min(colours), 1.05*max(colours)], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[0.95*min(colours), 1.05*max(colours)], ytickformat='(A1)'

; and the correlations
plotsym, 0, 0.5, /fill, color=0 ; default, small black dot
plot, shifts, correlations, xstyle=1, xrange=[min(shifts), max(shifts)], xtitle='Exposure shift (samples)', ystyle=1, ytitle='Correlation (in std. devs.)', yrange=[0.95*min(correlations), 1.05*max(correlations)], ymargin=[3.5, 0], charsize=charsize, charthick=charthick, psym=8, color=0
;; shading
;for l = 0, shifts_minimum-1 do begin
; oplot, [l, l], [0.95*min(correlations), 1.05*max(correlations)], thick=10, color=220
;endfor
; plot the "siren" case
oplot, timesteps, (correlations_max-correlations_mean)*siren + correlations_mean, color=0
; replot
oplot, shifts, correlations, psym=8, color=0
; overplotting the smoothed curve
oplot, shifts, correlations_smooth, thick=3, color=100
; maxima and minima
; plotsym, 0, 2., /fill, color=0 ; black filled circle, max
plotsym, 4, 3., thick=2, color=0 ; up-pointing triangle
oplot, [shifts_correlations_max, shifts_correlations_max], [correlations_max, correlations_max], psym=8, color=0
; plotsym, 0, 2., /fill, color=200 ; grey filled circle, min
plotsym, 5, 3., thick=2, color=0 ; down-pointing triangle 
oplot, [shifts_correlations_min, shifts_correlations_min], [correlations_min, correlations_min], psym=8, color=0
; limits
oplot, [shifts_minimum, shifts_minimum], [0.95*min(correlations), 1.05*max(correlations)], linestyle=1, color=0
oplot, [0, shifts_maximum], [correlations_mean, correlations_mean], linestyle=2, color=0
oplot, [0, shifts_maximum], [0.5, 0.5], linestyle=3, color=0 ; expected for no correlation
; clean up
axis, xaxis=0, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], charsize=charsize
axis, xaxis=1, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[0.95*min(correlations), 1.05*max(correlations)], charsize=charsize
axis, yaxis=1, color=0, ystyle=1, yrange=[0.95*min(correlations), 1.05*max(correlations)], ytickformat='(A1)'
 
; record these
print, 'Min/max colour difference of ', colours_min, ' and ', colours_max
print, 'Min/max correlation diff. of ', correlations_min, ' and ', correlations_max
 
; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_simulated_correlations'+label_signals+label_noise+label_samples+'.jpg', screen_output, quality=100, true=1
;write_gif, 'figure_simulated_correlations'+label_signals+label_noise+label_samples+'.gif', tvrd(screen_output)

end
