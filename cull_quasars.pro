pro cull_quasars

; Read in the MILLIQUAS quasars catalog, and cull to a Red magnitude, Blue-Red colour and redshift limit. Note that none of those three is very precise; perhaps within 0.2 mag, with often merely photometric redshifts.  The output is the catalog listed by magnitude and redshift limit, which is now concatenated to generate the full input catalog requred by report_gemini_pairs.
; Eric Steinbring, 29 January 2019 

; parameters
;magnitude_limit   = 20. ; 20. ; 19.5 ; 19. ; 18.5 ; 18. ; 17.5 ; 17. ; 16.5 ; 16. ; 15.5 ; 15. ; magnitudes, maximum Red mag of targets, which is now automatically spooled through
colour_limit      = 2. ; 0.5 ; 1., 2., magnitudes, upper limit for Blue-Red colour, which can be set to 2 here, and restricted later in the search for finding pairs
;redshift_limit    = 4. ; 0., 0.5, 1. ; 1.48 ; 1.5, 2., 2.5, 3., 3.5, 3.65 ; 4. ; 4.5, 5. ; 5.5, 6. ; minimum redshift limit

; as used previously, while testing the file formating
;readcol, 'list_quasars_test', format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift, /silent ;, skipline=1
;readcol, 'list_quasars_cut', format='A, A, A, F, F, F', identification, ra, dec, blue, red, redshift ; , /silent ;, skipline=1 ; a huge vector, making it too slow to work with
;input_number = n_elements(input_number)
;print, input_number

; read in the quasar list
;input_number = 1277820-1 ; 16-1
;openr, unit_in, 'list_quasars_test', /get_lun
;openr, unit_in, 'list_quasars', /get_lun
;openr, unit_in, 'list_quasars_cut', /get_lun ; same file, but without the header
magnitude_limits = [20., 19.5, 19., 18.5, 18.,  17.5, 17., 16.5, 16., 15.5, 15.] ; instead, a vector holding the magnitude limits for which to generate files
magnitudes = n_elements(magnitude_limits)
redshift_limits = [0., 0.5, 1., 1.48, 1.5, 2., 2.5, 3., 3.5, 3.65, 4., 4.5, 5., 5.5, 6.] ; instead, a vector holding the redshift limits for which to generate files
redshifts = n_elements(redshift_limits)
; loop on the redshift limit of the output file
for k = 0, redshifts - 1 do begin
 redshift_limit = redshift_limits(k)
 print, 'Redshift limit: ', redshift_limit
 ; loop on the the magnitude limit of the output file
 for j = 0, magnitudes  - 1 do begin
  magnitude_limit = magnitude_limits(j)
  print, 'Magnitude limit: ', magnitude_limit
  ; read in the quasar list; note that this must be re-read for each instance of redshift and magnitude, as the internal loop uses an end-of-file statement
  openr, unit_in, 'list_quasars_cut', /get_lun ; file without the header
  openw, unit_out, './list_quasars_z_r/list_quasars_zgt' + strmid(strtrim(string(redshift_limit), 2), 0, 4) + '_rlt' + strmid(strtrim(string(magnitude_limit), 2), 0, 4), /get_lun
  junk = ''
  identification_in = ' poop '
  ra_in = ' this '
  dec_in = ' that '
  blue_in = 0.
  red_in = 0.
  redshift_in = 0.
  number = 0
  stuff = ''
; for i = 0, input_number - 1 do begin ; looping on the total number of original entries
  while not eof(unit_in) do begin
   readf, unit_in, stuff
;  print, stuff
   length = strlen(stuff)
;   print, length
;   readf, unit_in, identification_in, ra_in, dec_in, blue_in, red_in, redshift_in
;   indentification_in = strmid(stuff, 0, 1) ; finding whether the identification is not blank
   blue_in = float(strmid(stuff, length-16, 4))
   red_in = float(strmid(stuff, length-11, 4))
   colour_in = blue_in-red_in
   redshift_in = float(strmid(stuff, length-6, 5)) 
;   print, red_in, redshift_in
   if red_in le magnitude_limit and colour_in lt colour_limit and redshift_in ge redshift_limit then begin
;   if identification_in ne '' and red_in le magnitude_limit and colour_in lt colour_limit and redshift_in ge redshift_limit then begin
;   if n_elements(identification_in) gt 2 and red_in le magnitude_limit and colour_in lt colour_limit and redshift_in ge redshift_limit then begin
    if number le 0 then number = floor(-1.*number) ; needed, as there is a bug which resets this counter to negative numbers, possibly by parsing a label incorrectly
    number = number + 1
    identification = identification_in
    ra = ra_in
    dec = dec_in
    blue = blue_in
    red = red_in
    redshift = redshift_in
;    print, number, identification, ra, dec, blue, red, redshift
;    printf, unit_out, identification, ra, dec, blue, red, redshift
;    print, stuff
    printf, unit_out, stuff
   endif
  endwhile
;  endfor
  close, unit_in
  free_lun, unit_in
  close, unit_out
  free_lun, unit_out
  print, 'Found ', number, ' quasars within limits'
 endfor ; looping on magnitude limit
endfor ; looping on redshift limit
;close, unit_in
;free_lun, unit_in
print, 'Done.'

end
