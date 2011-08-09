pro combine_corrmodesvd,field=field,mask=mask,residual=residual,onlysim=onlysim,sigcut5=sigcut5

if not keyword_set(field) then field='f3'

;smode=['1mode','2mode','3mode','4mode','5mode','6mode','7mode','8mode','9mode','10mode','11mode','12mode','15mode','20mode','25mode']
smode=['1mode','2mode','3mode','4mode','5mode','6mode','10mode','15mode','20mode','25mode']

nmode=n_elements(smode)

openw,5,'field3corr.tot.svdredisual.dat'
openw,6,field+'.residual.tot.signal.dat' 
openw,7,field+'.residual.tot.err.dat' 

;
;
for imode=0,nmode-1 do begin

correlate_gain_svdtest,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f3',mode=smode(imode),/test,mask=mask,residual=residual,onlysim=onlysim,sigcut5=sigcut5
printf,5,smode(imode)
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2,corre3
printf,5,'zero:',corrm1,corrm2
printf,5,''
printf,6,corrv1
printf,7,corre1
endfor

close,5
close,6
close,7

endif



;
if field eq 'f4' then begin

if keyword_set(sigcut5) then openw,5,'field4corr.tot.5sigcut.svdredisual.dat' else openw,5,'field4corr.tot.svdredisual.dat'
;
;
for imode=0,nmode-1 do begin

correlate_gain_svdtest,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,field='f4',mode=smode(imode),/test,mask=mask,residual=residual,onlysim=onlysim,sigcut5=sigcut5
printf,5,smode(imode)
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2,corre3
printf,5,'zero:',corrm1,corrm2
printf,5,''
printf,6,corrv1
printf,7,corre1

endfor

close, 5
close,6
close,7

endif




end
