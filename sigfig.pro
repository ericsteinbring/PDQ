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
