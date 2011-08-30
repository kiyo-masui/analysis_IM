pro combine_simcorr,field=field,sim=sim,sigcut5=sigcut5

if not keyword_set(field) then field='zcosmos'

;smode=['0mode','1mode','2mode','3mode','4mode','5mode','6mode','7mode','8mode','9mode','10mode','11mode','12mode','15mode','20mode','25mode']
smode=['1mode','2mode','3mode','4mode','5mode','6mode','10mode','15mode','20mode','25mode']
nmode=n_elements(smode)
;
openw,5,'zcosmoscorr.simonly.samew.tot.dat'
openw,6,field+'.simonly.samew.tot.signal.dat'
openw,7,field+'.simonly.samew.tot.err.dat'


for imode=0,nmode-1 do begin

correlate_gain_sim,corrv1,corrv2,corre1,corre2,corrm1,corrm2,field='f3',mode=smode(imode),sim=sim
printf,5,smode(imode)
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2
printf,5,'zero:',corrm1,corrm2
printf,5,''
printf,6,corrv1
printf,7,corre1
endfor

close, 5
close,6
close,7



end
