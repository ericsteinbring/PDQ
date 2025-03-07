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
