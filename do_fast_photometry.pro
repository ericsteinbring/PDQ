pro do_fast_photometry

; This program reads in the available test dataset of Alopeke/Zorro pairs which have near-simultaneous observations. These are standard datacubes of Red and Blue filter samples of 1000 frames at 0.06 seconds integration, and the code performs automated aperture photometry on those. The output is the calculated correlated S/N that is reported in "find_quasar_pairs."
; Eric Steinbring, 27 September 2024

; parameters
diameter             = 8.1 ; m, of the telescope aperture
obstruction          = 1. ; m, diameter of central obstruction
pixel_default        = 0.010 ; arcsec/pix, default pixel scale of detectors
gain_default         = 0.8 ; electrons/count, default gain of detectors
lambda_blue          = 562. ; nm, central wavelength of blue filter
delta_lambda_blue    = 54. ; nm, approximate bandpass of blue filter 
lambda_red           = 832. ; nm, central wavelength of red filter
delta_lambda_red     = 40. ; nm, approximate bandpass of red filter
efficiency_blue      = 0.90 ; quantum efficiency of detectors at blue central wavelength of (and including) narrowband filter
efficiency_red       = 0.70 ; quantum efficiency of detectors at red central wavelength of (and including) narrowband filter
exposure_default     = 60. ; seconds
field_size           = 256 ; pixels
frames               = 1000 ; per cube
photometric_aperture = 5. ; 1. ; 2. ; 3. ; 4. ; 5. ; arcsec
photometric_annulus  = 1. ; 1. ; 2. ; 3. ; 4. ; 5. ; arcsec
time_minimum         = 3. ; 1. ; 2. ; 3. ; 4. ; 5. ; seconds, the largest difference in timestamp between North and South frames
shifts_maximum       = 250 ; 100 ; 250 ; 500 ; 1000 ; per cube, maximum allowed shift in time, looking for variation in the overall correlation of frames

; display setups
display_size         = 400 ; pixels, postage stamp size for the frame display panel
stop_display         = 'no' ; 'no' ; 'yes' ; stop after first observation to show default display setup

; setup display
window, 21, xsize = 1600, ysize = 500, title='Alopeke/Zorro frames'
wset, 21
loadct, 0, /silent
!p.background=255
device, set_font='helvetica bold', /tt_font
!p.font=1
!except=0
charsize=2.5
charthick=2.
!p.multi=[0,1,1]

; definitions
c = 299792458 ; m/s

; calculations
area = 100.^2.*!pi*((diameter-obstruction)/2.)^2. ; cm^2, telescope pupil aperture
exposure_time = exposure_default/float(frames)
print, 'Exposure time of individual frame (seconds): ', exposure_time
shifts_minimum = time_minimum/exposure_default*frames
print, 'Minimum shift in exposures between North and South frames: ', shifts_minimum

; open a file for output colour and correlation measures
get_lun, unit0
openw, unit0, 'signals_test'
get_lun, unit1
openw, unit1, 'correlations_test'

; make an aperture and annulus
field = photometric_aperture + 2.*photometric_annulus ; arcseconds
r = field_size/2. ; 50 ; 100 ; pixels, size of postage stamp
;r = floor(field/(2.*pixel_default)); pixels, set instead from the photometric aperture and annulus
;print, photometric_aperture, field, photometric_aperture/field
blank = fltarr(2*r, 2*r)
aperture = float(shift((dist(2.*r)*2./(2.*r) le photometric_aperture/field and dist(2.*r)*2./(2.*r) ge 0.), float(r), float(r)))
;print, min(aperture), max(aperture), total(aperture)
area_aperture = total(aperture) ; in pixels
tvscl, congrid(aperture, display_size, display_size), 0, 0
annulus = float(shift((dist(2.*r)*2./(2.*r) gt photometric_aperture/field and dist(2.*r)*2./(2.*r) lt (photometric_aperture + 2.*photometric_annulus)/field), float(r), float(r)))
annulus[2*r-1:2*r-1, *] = 0.
annulus[*, 2*r-1:2*r-1] = 0.
;print, min(annulus), max(annulus), total(annulus)
tvscl, congrid(annulus, display_size, display_size), display_size, 0
;stop

; make the frame background white for the display
blank = fltarr(1600, 500)+255
tv, blank

; label the display
xyouts, 10, display_size-20+100, 'Gemini-North `Alopeke', charsize=charsize, charthick=charthick, color=0, /device
xyouts, 10, 10, 'Blue', charsize=charsize, charthick=charthick, color=0, /device
xyouts, display_size+10, 10, 'Red', charsize=charsize, charthick=charthick, color=0, /device
xyouts, 2.*display_size+10, display_size-20+100, 'Gemini-South Zorro', charsize=charsize, charthick=charthick, color=0, /device
xyouts, 2.*display_size+10, 10, 'Blue', charsize=charsize, charthick=charthick, color=0, /device
xyouts, 3.*display_size+10, 10, 'Red', charsize=charsize, charthick=charthick, color=0, /device

; go through the list of frame labels
readcol, './Data/inlist_sorted', format='A', inlist, /silent
;readcol, 'inlist_sorted', format='A', inlist, /silent
sequence = n_elements(inlist)/4 ; data come in sequences of 4 cubes: North/South Blue/Red
print, 'Sequences of data in the list: ', sequence 
for j = 0, sequence-1 do begin

 ; read in the set of frames
 inlist_n_b = inlist(j*4+0)
 inlist_n_r = inlist(j*4+1)
 inlist_s_b = inlist(j*4+2)
 inlist_s_r = inlist(j*4+3) 

 ; read in the frames
; frames_n_b = readfits('./Data/N20211020A0659b.fits', /silent)
; frames_n_r = readfits('./Data/N20211020A0659r.fits', /silent)
; frames_s_b = readfits('./Data/S20211021Z0294b.fits', /silent)
; frames_s_r = readfits('./Data/S20211021Z0294r.fits', /silent)
 frames_n_b = readfits('./Data/' + inlist_n_b, /silent)
 frames_n_r = readfits('./Data/' + inlist_n_r, /silent)
 frames_s_b = readfits('./Data/' + inlist_s_b, /silent)
 frames_s_r = readfits('./Data/' + inlist_s_r, /silent)
; frames_n_b = readfits(inlist_n_b, /silent)
; frames_n_r = readfits(inlist_n_r, /silent)
; frames_s_b = readfits(inlist_s_b, /silent)
; frames_s_r = readfits(inlist_s_r, /silent)
; print, size(frames_n_b)
; print, size(frames_n_r)
; print, size(frames_s_b)
; print, size(frames_s_r)

 ; go through all the frames in a cube
 frame_n_b = fltarr(field_size, field_size)
 frame_n_r = frame_n_b
 frame_s_b = frame_n_b
 frame_s_r = frame_n_b
 flux_n_b = fltarr(frames)
 flux_n_r = flux_n_b
 flux_s_b = flux_n_b
 flux_s_r = flux_n_b
 fluence_n_b = flux_n_b
 fluence_n_r = flux_n_b
 fluence_s_b = flux_n_b
 fluence_s_r = flux_n_b
 mag_n_b = flux_n_b
 mag_n_r = flux_n_b
 mag_s_b = flux_n_b
 mag_s_r = flux_n_b
 for i = 0, frames-1 do begin

  ; take one set of four individual frames, North/South and Blue/Red from the cube
  frame_n_b(*, *) = frames_n_b(*, *, i)
  frame_n_r(*, *) = frames_n_r(*, *, i)
  frame_s_b(*, *) = frames_s_b(*, *, i)
  frame_s_r(*, *) = frames_s_r(*, *, i)
  
  ; flip and rotate to align these to North up, East left; Blue is flipped in original dataset, and both need to be rotated 90 degrees
  frame_n_b = reverse(frame_n_b) ; flip in x-axis, that is, the rows
  frame_s_b = reverse(frame_s_b) ; flip in x-axis, that is, the rows
  ; rotate for East to the left, as displayed
  frame_n_b = rotate(frame_n_b, 3) ; rotate 90 degrees clockwise, or 270 counterclockwise
  frame_n_r = rotate(frame_n_r, 3) ; rotate 90 degrees clockwise, or 270 counterclockwise
  frame_s_b = rotate(frame_s_b, 3) ; rotate 90 degrees clockwise, or 270 counterclockwise
  frame_s_r = rotate(frame_s_r, 3) ; rotate 90 degrees clockwise, or 270 counterclockwise

  ; display those
  wset, 21 ; setting back to the frame display panel
  ; first, inverting it for display
  frame_n_b = 1.-(frame_n_b/max(frame_n_b))
  frame_n_r = 1.-(frame_n_r/max(frame_n_r))
  frame_s_b = 1.-(frame_s_b/max(frame_s_b))
  frame_s_r = 1.-(frame_s_r/max(frame_s_r)) 
  tvscl, congrid(frame_n_b, display_size, display_size), 0, 50
  tvscl, congrid(frame_n_r, display_size, display_size), display_size, 50
  tvscl, congrid(frame_s_b, display_size, display_size), 2.*display_size, 50
  tvscl, congrid(frame_s_r, display_size, display_size), 3.*display_size, 50
 
  ; do photometry
;  flux_n_b(i) = total(aperture*frame_n_b)-mean(annulus*frame_n_b)*area_aperture ; in counts
;  flux_n_r(i) = total(aperture*frame_n_r)-mean(annulus*frame_n_r)*area_aperture ; in counts
;  flux_s_b(i) = total(aperture*frame_s_b)-mean(annulus*frame_s_b)*area_aperture ; in counts
;  flux_s_r(i) = total(aperture*frame_s_r)-mean(annulus*frame_s_r)*area_aperture ; in counts
  flux_n_b(i) = total(aperture*frame_n_b)-median(annulus*frame_n_b)*area_aperture ; in counts
  flux_n_r(i) = total(aperture*frame_n_r)-median(annulus*frame_n_r)*area_aperture ; in counts
  flux_s_b(i) = total(aperture*frame_s_b)-median(annulus*frame_s_b)*area_aperture ; in counts
  flux_s_r(i) = total(aperture*frame_s_r)-median(annulus*frame_s_r)*area_aperture ; in counts
;  print, 'Fluxes (counts): ', flux_n_b(i), flux_n_r(i), flux_s_b(i), flux_s_r(i)

  ; report the fluence, so correcting for efficiency
  fluence_n_b(i) = (1./efficiency_blue)*gain_default*flux_n_b(i) ; in photons per exposure
  fluence_n_r(i) = (1./efficiency_red)*gain_default*flux_n_r(i) ; in photons per exposure
  fluence_s_b(i) = (1./efficiency_blue)*gain_default*flux_s_b(i) ; in photons per exposure
  fluence_s_r(i) = (1./efficiency_red)*gain_default*flux_s_r(i) ; in photons per exposure
;  print, 'Fluence (photons/exposure): ', fluence_n_b(i), fluence_n_r(i), fluence_s_b(i), fluence_r(i)

  ; correct for gain, efficiency, exposure time and telescope aperture
  flux_n_b(i) = (1./efficiency_blue)*gain_default*flux_n_b(i)/(exposure_time*area) ; flux in photons per second per cm^2 at the filter central wavelength
  flux_n_r(i) = (1./efficiency_red)*gain_default*flux_n_r(i)/(exposure_time*area) ; flux in photons per second per cm^2 at the filter central wavelength
  flux_s_b(i) = (1./efficiency_blue)*gain_default*flux_s_b(i)/(exposure_time*area) ; flux in photons per second per cm^2 at the filter central wavelength
  flux_s_r(i) = (1./efficiency_red)*gain_default*flux_s_r(i)/(exposure_time*area) ; flux in photons per second per cm^2 at the filter central wavelength
;  print, 'Fluxes (photons/second/cm^2): ', flux_n_b(i), flux_n_r(i), flux_s_b(i), flux_s_r(i)
  
;  ; convert to flux per unit frequency
;  conversion = 6.6225e-26
;  flux_n_b(i) = conversion*flux_n_b(i)*(delta_lambda_blue/lambda_blue) ; flux per unit frequency
;  flux_n_r(i) = conversion*flux_n_r(i)*(delta_lambda_red/lambda_red) ; flux per unit frequency
;  flux_s_b(i) = conversion*flux_s_b(i)*(delta_lambda_blue/lambda_blue) ; flux per unit frequency
;  flux_s_r(i) = conversion*flux_s_r(i)*(delta_lambda_red/lambda_red) ; flux per unit frequency 
;;  print, 'Fluxes per unit frequency (f_nu): ', flux_n_b(i), flux_n_r(i), flux_s_b(i), flux_s_r(i)  

;  ; convert to AB magnitudes, from the definition
;  mag_n_b(i) = -2.5*alog10(flux_n_b(i))-48.6 ; magnitudes
;  mag_n_r(i) = -2.5*alog10(flux_n_r(i))-48.6 ; magnitudes
;  mag_s_b(i) = -2.5*alog10(flux_s_b(i))-48.6 ; magnitudes
;  mag_s_r(i) = -2.5*alog10(flux_s_r(i))-48.6 ; magnitudes
  ; kludge: bypassing this by taking the input fluxes directly in photons instead, and just using an identity
  mag_n_b(i) = -2.5*alog10((flux_n_b(i)/10.*(lambda_blue/1.e9)^2./c))-48.6 ; magnitudes
  mag_n_r(i) = -2.5*alog10((flux_n_r(i)/10.*(lambda_red/1.e9)^2./c))-48.6 ; magnitudes
  mag_s_b(i) = -2.5*alog10((flux_s_b(i)/10.*(lambda_blue/1.e9)^2./c))-48.6 ; magnitudes
  mag_s_r(i) = -2.5*alog10((flux_s_r(i)/10.*(lambda_red/1.e9)^2./c))-48.6 ; magnitudes
;  print, 'Magnitudes (AB): ', mag_n_b(i), mag_n_r(i), mag_s_b(i), mag_s_r(i)

 endfor
 
 ; save a screen capture
 screen_output=tvrd(true=1)
 if j eq 0 then write_jpeg, 'figure_Alopeke_Zorro_frames.jpg', screen_output, quality=100, true=1
 
 ; and report statistics of S/N ; and report statistics of S/N
 fluxes = [flux_n_b, flux_n_r, flux_s_b, flux_s_r]
 fluxes_mean = mean(fluxes, /nan)
 fluxes_median = median(fluxes)
 fluxes_stddev = stddev(fluxes, /nan)
 fluxes_max = max(fluxes, /nan)
 fluxes_min = min(fluxes, /nan)
 snr_min = fluxes_min/fluxes_stddev
 snr_max = fluxes_max/fluxes_stddev
 snr_median = fluxes_median/fluxes_stddev
 snr_mean = fluxes_mean/fluxes_stddev
 print, 'SNR, overall: ', snr_min, snr_median, snr_mean, snr_max
 
; ; report magnitudes
; print, 'Mean magnitudes, N (Red/Blue) and S (Red/Blue) (AB): ', mean(mag_n_b, /nan), mean(mag_n_r, /nan), mean(mag_s_b, /nan), mean(mag_s_r, /nan)
; print, 'Estimated V magnitudes N and S (mean Red/Blue): ', (mean(mag_n_b, /nan) + mean(mag_n_r, /nan))/2., (mean(mag_s_b, /nan) + mean(mag_s_r, /nan))/2.

 ; calculate the relative colours for North and South
; colour_n = flux_n_b - flux_n_r ; in either counts or photons per second per cm^2
; colour_s = flux_s_b - flux_s_r ; in either counts or photons per second per cm^2
 colour_n = fluence_n_b - fluence_n_r ; in photons per exposure
 colour_s = fluence_s_b - fluence_s_r ; in photons per exposure
; colour_n = mag_n_b - mag_n_r ; magnitudes
; colour_s = mag_s_b - mag_s_r ; magnitudes
 
 ; take the colour difference
 colour = colour_n - colour_s
  
 ; report statistics
 colour_mean = mean(colour)
 colour_median = median(colour)
 colour_stddev = stddev(colour)
;  print, 'Sequence colour mean, median and std. dev.: ', j+1, colour_mean, colour_median, colour_stddev
 
 ; and normalize 
; colour = (colour - colour_median)/colour_mean ; normalized to the mean
 colour = (colour - colour_median)/colour_stddev ; normalized to one standard deviation
 
 ; look for a trend
 frame_number = findgen(frames)
 trend = poly_fit(frame_number, colour, 1, sigma=sigma, yfit=colour_fit) ; just linear, to start
; print, 'Slope of colour-change per frame: ', trend(1)

 ; and subtract it from the data
 colour = colour - colour_fit
 
  ; look for a trend, again
 frame_number = findgen(frames)
 trend = poly_fit(frame_number, colour, 1, sigma=sigma, yfit=colour_fit) ; just linear, to start
; print, 'Slope of colour-change per frame: ', trend(1)

 ; set up the display
 window, 22, xsize = 1500, ysize = 500, title='Relative stellar colours'
 wset, 22
 !p.multi = [0, 1, 1]
 charsize=2.5

 ; plot the relative colours
 plotsym, 0, 2., /fill, color=0 ; black filled circle, default
 plot, colour, xtitle='Frame number', ytitle='Rel. A-B colour (in std. dev. from mean)', ystyle=1, yrange=[-4., 4.], charsize=charsize, charthick=charthick, psym=8, color=0
 ; shading
 for l = 0, frames-1 do begin
  oplot, [l, l], [-1., 1.], thick=5, color=220
 endfor
 ; replot
 oplot, colour, psym=8, color=0
 ; limits
; oplot, frame_number, colour_fit, thick=3, color=100
; oplot, [0, frames], [0., 0.], linestyle=2, color=0
 oplot, [0, frames], [0., 0.], linestyle=3, color=0
 ; clean up
 axis, xaxis=0, color=0, xstyle=1, xrange=[0, frames], charsize=charsize
 axis, xaxis=1, color=0, xstyle=1, xrange=[0, frames], xtickformat='(A1)'
 axis, yaxis=0, color=0, ystyle=1, yrange=[-4., 4.], charsize=charsize
 axis, yaxis=1, color=0, ystyle=1, yrange=[-4., 4.], ytickformat='(A1)'
 
 ; save a screen capture
 screen_output=tvrd(true=1)
 if j eq 0 then write_jpeg, 'figure_Alopeke_Zorro_relative_colour.jpg', screen_output, quality=100, true=1
 
 ; setup display
 window, 23, xsize = 600, ysize = 500, xpos=1400, ypos=550, title='Distribution of Signal'
 wset, 23
 !p.multi=[0,1,1]

 ; find the distribution of colours, and signal
 binsize=0.25 ; 0.1 ; 0.2 ; 0.25 ; 0.5
 distribution_min = -4.
 distribution_max = 4.
 distribution_colour = histogram(colour, min=distribution_min, max=distribution_max, binsize=binsize)
 bins = binsize*findgen(n_elements(distribution_colour)) + distribution_min

 ; smooth
 ;distribution_colour = smooth(distribution_colour, 3)

 ; normalize
 total_distribution = float(total(distribution_colour))
 distribution_colour = distribution_colour/total_distribution
 max_distribution = float(max(distribution_colour))
 
 ; and store it
 if j eq 0 then distribution_colours = distribution_colour ; initializing
 if j ge 1 then distribution_colours = distribution_colours + distribution_colour ; adding them up

 ; report that of a pure, normal distribution
 normal = gaussian(bins, [max_distribution, 0., 1.]) ; set to the peak of the colour distribution
; normal = gaussian(bins, [0.5, 0., 1.]) ; set to be a peak of 50%
; normal = gaussian(bins, [0.25, 0., 1.]) ; set to be a peak of 25%

 ; and plot them
 plot, distribution_colour, bins+binsize/2., xstyle=1, xtitle='Normalized occurances', xrange=[0., 0.25], xticks=5, ystyle=1, yrange=[distribution_min, distribution_max], ytitle='Distribution of A-B colour', charsize=charsize, charthick=charthick, thick=3, color=0 ;, ytickformat='(A1)'
 ; shading
 for i = 0, n_elements(bins)-1 do begin
  if bins(i) gt -1. and bins(i) lt 1. then oplot, [0., 0.5], [bins(i), bins(i)], thick=20, color=220
 endfor
 ; replot
 oplot, distribution_colour, bins+binsize/2., thick=3, color=0
 ; overplot a Gaussian
 oplot, normal, bins, thick=3, linestyle=2, color=100
 ; limits
 oplot, [0., 0.25], [0., 0.], linestyle=3, color=0
 ; labels
 xyouts, 0.025, 3.25, 'Single sequence: 1000 samples', charsize=charsize, charthick=charthick, color=0
 ; clean up
; axis, xaxis=0, color=0, xstyle=1, xrange=[0., 0.5], charsize=charsize
; axis, xaxis=1, color=0, xstyle=1, xrange=[0., 0.5], xtickformat='(A1)'
 axis, yaxis=0, color=0, ystyle=1, yrange=[distribution_min, distribution_max], charsize=charsize ;, ytickformat='(A1)'
 axis, yaxis=1, color=0, ystyle=1, yrange=[distribution_min, distribution_max], ytickformat='(A1)'

 ; save a screen capture
 screen_output=tvrd(true=1)
 if j eq 0 then write_jpeg, 'figure_Alopeke_Zorro_distribution_colour.jpg', screen_output, quality=100, true=1
; if j ge 1 then write_jpeg, 'figure_Alopeke_Zorro_distribution_colour_'+strtrim(string(j),2)+'.jpg', screen_output, quality=100, true=1
 
 ; now shift the relative photometry by time, up to the prescribed limit, and look for a maximum and minimum of correlation
 shifts = findgen(shifts_maximum) ; start at none and shift South relative to North
; shifts = -floor(shifts_maximum/2.) + findgen(shifts_maximum) ; start by goning back half the maximum number of shifts
 colours = fltarr(n_elements(shifts))
 correlations = colours
 for k = 0, n_elements(shifts)-1 do begin
 
  ; shift the photometry
;  colour_n = shift(colour_n, shifts(k)) ; shift by number of exposures
  colour_s = shift(colour_s, shifts(k)) ; shift by number of exposures
 
  ; take the colour difference
  colour = colour_n - colour_s
  
;  ; report statistics
;  colour_mean = mean(colour)
;  colour_median = median(colour)
;  colour_stddev = stddev(colour)
;;  print, 'Sequence and shift rel. colour mean, median and std. dev.: ', j+1, k+1, colour_mean, colour_median, colour_stddev

  ; normalize 
;  colour = (colour - colour_median)/colour_mean ; normalized to the mean
  colour = (colour - colour_median)/colour_stddev ; normalized to one standard deviation
 
  ; and record the degree of correlation, e.g., recording the colour standard deviation for the given shift
  colours(k) = mean(abs(colour))
  correlations(k) = stddev(colour)/2.
 
 endfor
 
; ; take absolutes or normalize
; colours = abs(colours)
; colours = colours-mean(colours)
; colours = abs(colours)
; colours = colours/mean(colours)
; correlations = abs(correlations)
; correlations = correlations-mean(correlations)
; correlations = abs(correlations)
; correlations = correlations/mean(correlations)

 ; a smoothed fit
; colours_smooth = smooth(colours, floor(shifts_maximum/10.), /edge_wrap)
; correlations_smooth = smooth(correlations, floor(shifts_maximum/10.), /edge_wrap)
 colours_smooth = smooth(colours, shifts_minimum, /edge_wrap)
 correlations_smooth = smooth(correlations, shifts_minimum, /edge_wrap)
 
 ; find the means
 colours_mean = mean(colours_smooth, /nan)
 correlations_mean = mean(correlations_smooth, /nan)
 
; ; find the maximum and minimum, overall; if multiple, take the first one
; colours_max = max(colours)
; shifts_colours_max = shifts(min(where(colours eq max(colours))))
; correlations_max = max(correlations)
; shifts_correlations_max = shifts(min(where(correlations eq max(correlations))))
; colours_min = min(colours)
; shifts_colours_min = shifts(min(where(colours eq min(colours))))
; correlations_min = min(correlations)
; shifts_correlations_min = shifts(min(where(correlations eq min(correlations))))

; ; find the maximum and minimum, after smoothing; if multiple, take the first one
; colours_max = max(colours_smooth)
; shifts_colours_max = shifts(min(where(colours_smooth eq max(colours_smooth))))
; correlations_max = max(correlations_smooth)
; shifts_correlations_max = shifts(min(where(correlations_smooth eq max(correlations_smooth))))
; colours_min = min(colours_smooth)
; shifts_colours_min = shifts(min(where(colours_smooth eq min(colours_smooth))))
; correlations_min = min(correlations_smooth)
; shifts_correlations_min = shifts(min(where(correlations_smooth eq min(correlations_smooth))))
 
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
 window, 24, xsize = 1000, ysize = 1000, title='Correlation of stellar colours'
 wset, 24
 !p.multi = [0, 1, 2]
 charsize=2.5
 
 ; plot the colour differences
 plotsym, 0, 0.5, /fill, color=0 ; default, small black dot
 plot, shifts, colours, xstyle=1, xrange=[min(shifts), max(shifts)], xtickformat='(A1)', ystyle=1, ytitle='Abs. colour diff. (in std. dev. rel. to mean)', yrange=[0.95*min(colours), 1.05*max(colours)], ymargin=[0,2], charsize=charsize, charthick=charthick, psym=8, color=0
; ; shading
; for l = 0, shifts_minimum-1 do begin
;  oplot, [l, l], [0.95*min(colours), 1.05*max(colours)], thick=10, color=220
; endfor
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
 oplot, [shifts_minimum, shifts_minimum], [0.75*min(colours), 1.25*max(colours)], linestyle=1, color=0
 oplot, [0, shifts_maximum], [colours_mean, colours_mean], linestyle=2, color=0
 oplot, [0, shifts_maximum], [0.75, 0.75], linestyle=3, color=0 ; expected for no correlation
 ; clean up
; axis, xaxis=0, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], charsize=charsize
 axis, xaxis=1, color=0, xstyle=1, xrange=[min(shifts), max(shifts)], xtickformat='(A1)'
 axis, yaxis=0, color=0, ystyle=1, yrange=[0.95*min(colours), 1.05*max(colours)], charsize=charsize
 axis, yaxis=1, color=0, ystyle=1, yrange=[0.95*min(colours), 1.05*max(colours)], ytickformat='(A1)'

 ; and the correlations
 plotsym, 0, 0.5, /fill, color=0 ; default, small black dot
 plot, shifts, correlations, xstyle=1, xrange=[min(shifts), max(shifts)], xtitle='Exposure shift (frames)', ystyle=1, ytitle='Correlation (in std. devs.)', yrange=[0.95*min(correlations), 1.05*max(correlations)], ymargin=[3.5, 0], charsize=charsize, charthick=charthick, psym=8, color=0
; ; shading
; for l = 0, shifts_minimum-1 do begin
;  oplot, [l, l], [0.95*min(correlations), 1.05*max(correlations)], thick=10, color=220
; endfor
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
 print, 'Sequence ' , j+1, ' min/max colour difference of ', colours_min, ' and ', colours_max
 print, 'Sequence ' , j+1, ' min/max correlation diff. of ', correlations_min, ' and ', correlations_max
 
 ; save a screen capture
 screen_output=tvrd(true=1)
 if j eq 0 then write_jpeg, 'figure_Alopeke_Zorro_correlations.jpg', screen_output, quality=100, true=1
 
 ; write them out to a file for the Gemini-pairs reporting code
 printf, unit0, format='(4f)', snr_min, snr_median, snr_mean, snr_max
 printf, unit1, format='(4f)', colours_min, colours_max, correlations_min, correlations_max
 
 ; show display
 if stop_display eq 'yes' then stop

endfor
close, unit0 ; the file containing measured S/N outputs
close, unit1 ; the file containing measured correlation outputs

; setup display
window, 25, xsize = 600, ysize = 500, xpos=1400, ypos=550, title='Distribution of Signals (Averaged)'
wset, 25
!p.multi=[0,1,1]

; normalize
;distribution_colours = distribution_colours/float(sequence) ; an average over the total sequence of observations
total_distribution = float(total(distribution_colours))
distribution_colours = distribution_colours/total_distribution
max_distribution = float(max(distribution_colours))

; report that of a pure, normal distribution
normal = gaussian(bins, [max_distribution, 0., 1.]) ; set to the peak of the colour distribution
; normal = gaussian(bins, [0.5, 0., 1.]) ; set to be a peak of 50%
; normal = gaussian(bins, [0.25, 0., 1.]) ; set to be a peak of 25%

; display the total distribution
plot, distribution_colours, bins+binsize/2., xstyle=1, xtitle='Normalized occurances', xrange=[0., 0.25], xticks=5, ystyle=1, yrange=[distribution_min, distribution_max], ytitle='Distribution of A-B colour', charsize=charsize, charthick=charthick, thick=3, color=0 ;, ytickformat='(A1)'
; shading
for i = 0, n_elements(bins)-1 do begin
 if bins(i) gt -1. and bins(i) lt 1. then oplot, [0., 0.5], [bins(i), bins(i)], thick=20, color=220
endfor
; replot
oplot, distribution_colours, bins+binsize/2., thick=3, color=0
; overplot a Gaussian
oplot, normal, bins, thick=3, linestyle=2, color=100
; limits
oplot, [0., 0.25], [0., 0.], linestyle=3, color=0
; labels
xyouts, 0.025, 3.25, 'Total of all samples', charsize=charsize, charthick=charthick, color=0
; clean up
;axis, xaxis=0, color=0, xstyle=1, xrange=[0., 0.5], charsize=charsize
;axis, xaxis=1, color=0, xstyle=1, xrange=[0., 0.5], xtickformat='(A1)'
axis, yaxis=0, color=0, ystyle=1, yrange=[distribution_min, distribution_max], charsize=charsize ;, ytickformat='(A1)'
axis, yaxis=1, color=0, ystyle=1, yrange=[distribution_min, distribution_max], ytickformat='(A1)'

; save a screen capture
screen_output=tvrd(true=1)
write_jpeg, 'figure_Alopeke_Zorro_distribution_colours.jpg', screen_output, quality=100, true=1

end

