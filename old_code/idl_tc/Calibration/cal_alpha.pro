pro cal_alpha

freq_flux=[ $
;           [1.40E+09,	 4.08E+01],$
           [9.60E+08,	 6.54E+01],$
           [7.50E+08,	 8.36E+01],$
           [6.35E+08,	 9.71E+01],$
           [4.68E+08,	 1.15E+02],$
           [4.08E+08,	 1.32E+02],$
;           [3.65E+08,	 7.29E+01],$
           [1.60E+08,	 2.46E+02]]

freq=freq_flux(0,*)
print,freq
flux=freq_flux(1,*)
print,'flux',flux

freq=freq/1e6
flux1=alog10(flux)
freq1=alog10(freq)


coef=linfit(freq1,flux1,/double,yfit=yfit)

print,'const',coef(0)
print,'alpha',coef(1)
print,'flux',reform(flux)
print,'fitted flux',10^(yfit)
print,'fitted',(freq/freq(0))^(coef(1))*flux(0)
print,'fitted',10^(coef(0)+coef(1)*freq1)

;alpha=-0.79949634
;const       19.159453
;alpha     -0.72159868

end

