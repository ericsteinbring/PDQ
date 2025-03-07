; Bell theorem
function bell_theorem, theta
; Returns the two-point quantum correlation of point sources, as in Bell's Theorem. Theta is in degrees.
;theta = theta-180. ; test of mirroring
;return, -1.*abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.) ; test of mirroring
return, abs(cos(2.*!pi*theta/180.) - cos(!pi*theta/180.)) + cos(!pi*theta/180.)
end
