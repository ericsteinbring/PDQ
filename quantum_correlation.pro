; quantum correlation
function quantum_correlation, theta
; Returns the "classical" two-point quantum correlation of point sources. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*(-1.*cos(!pi*theta/180.)) ; test of mirroring
return, -1.*cos(!pi*theta/180.)
end
