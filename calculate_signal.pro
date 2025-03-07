pro calculate_signal

; Estimates the pairwise signal-to-noise ratio for relative photometry of two quasars with exactly equal magnitudes, a calculation restricted only by zeropoint of filters, efficiency of detectors, size of telescope and extinction.  The number of samples is an input to the final, combined or stacked-exposures result, allowing it to report a reasonable total integration time to achieve this, assuming Poisson statistics.
; Eric Steinbring, 29 January 2019.

; parameters
diameter     = 8. ; 1.12838 ; 4. ; 8. ; 8.1 ; 30 ; m, of the telescope aperture
obstruction  = 0. ; 0.25 ; 0.5 ; 1. ; 2. ; m, diameter of central obstruction
lambda_c     = 0.55 ; 0.64 ; 0.55 ; microns, central wavelength of filter
delta_lambda = 0.16 ; 0.23 ; 0.16 ; microns, bandwidth of filter
flux_0       = 3631. ; 3080. ; 3640. ; 3631. ; Jy ; zeropoint in flux units, where 3631 is set by the definition of the AB-magnitude zeropoint
exposure     = 0.01 ; 0.001 ; 0.005 ; 0.01 ; 0.015 ; 0.02 ; 0.03 ; 0.06 ;  0.1 ; 1. ; exposure time in seconds
overhead     = 0.25 ; 5. ; 2. ; 1. ; 0.5 ; 0.25 ;  0. ; factor of exposure time in addition, as overhead per exposure, unity means that the overhead is equal to the exposure time
efficiency   = 0.5 ; 0.5 ; 0.8 ; 0.85 ; 0.9 ; 1. ; fractional efficiency of detector
magnitude    = 20. ; 18.5 ; 15. ; 15.5 ; 16. ; 16.5 ; 17. ; 17.5 ; 18. ; 18.5 ; 19. ; 19.5 ; 20. ; mag, of sources; default is 18.5 for roughly a redshift 4 quasar
psf_fwhm     = 0.8 ; 0.5 ; 0.65 ; 0.8 ; 1. ; arcsec, FWHM of image, assumed stellar
extinction   = 1. ; 0. ; 0.1 ; 0.2 ; 0.25 ; 0.5 ; 1. ; 2. ; mag, integrated over filter bandpass
sky          = 20. ; 15. ; 16. ; 17. ; 18. ; 19. ; 20. ; 21. ; 22. ; mag arcsec^-2, sky brightness integrated over filter bandpass
samples      = 1000 ; 1.e6 ; 0.5e6 ; 1.e5 ; 10000 ; 1000 ; 100 ; 10 ; samples of each exposure, default is 1000

; calculations
area = !pi*((diameter-obstruction)/2.)^2. ; m^2
source = flux_0*10.^(-0.4*(magnitude+extinction))*1.51e7*delta_lambda ; photons sec^-1 m^-2
flux = area*source ; photons sec^-1
fluence = flux*exposure ; photons per exposure time
;magnitude_limit = ? ; mag, where fluence is one photon
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
print, 'Individual S/N: ', signal/noise, ' or ', 100.*(1./(signal/noise)), ' percent error'
print, 'Approx. photometric error (mag): ', 1./(signal/noise)
print, 'Relative pairwise S/N: ', signal/(sqrt(2.)*noise), ' or ', 100.*(1./(signal/(sqrt(2.)*noise))), ' percent error'
print, 'Samples: ', samples
print, 'Final pairwise S/N: ', sqrt(float(samples))*signal/(sqrt(2.)*noise)
print, 'Integration (s, m, h): ', float(samples)*exposure*(1.+overhead), float(samples)*exposure*(1.+overhead)/60., float(samples)*exposure*(1.+overhead)/60./60.

end
