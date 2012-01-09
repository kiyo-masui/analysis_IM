pro combine_corrmode2,field=field,sim=sim,dname=dname,drift=drift,daisy=daisy,onlysim=onlysim

if not keyword_set(field) then field='zcosmos'
if not keyword_set(dname) then dname='02'
if keyword_set(daisy) then obsmode='daisy' 
if not keyword_set(daisy) then drift=drift
if keyword_set(drift) then obsmode='drift'
if keyword_set(onlysim) then flag='.simonly' else flag=''
;if not keyword_set(onlysim) then gain=gain

;smode=['0mode','1mode','2mode','3mode','4mode','5mode','6mode','7mode','8mode','9mode','10mode','11mode','12mode','15mode','20mode','25mode']
smode=['0mode','1mode','2mode','3mode','4mode','5mode','6mode','10mode','15mode','20mode','25mode']
nmode=n_elements(smode)
;

if keyword_set(sim) then openw,6,field+'.sim.tot.signal.'+obsmode+'.dat' else openw,6,field+'.tot.signal.'+obsmode+flag+'.dat'
if keyword_set(sim) then openw,7,field+'.sim.tot.err.'+obsmode+'.dat' else openw,7,field+'.tot.err.'+obsmode+flag+'.dat'
if keyword_set(sim) then openw,8,field+'.sim.tot.autocorr.'+obsmode+'.dat' else openw,8,field+'.tot.autocorr.'+obsmode+flag+'.dat'
if keyword_set(sim) then openw,5,field+'.corr.sim.tot.'+obsmode+'.dat' else openw,5,field+'.corr.tot.'+obsmode+flag+'.dat'

for imode=0,nmode-1 do begin

correlate_gain,corrv1,corrv2,corre1,corre2,corre3,corrm1,corrm2,autocorr,mode=smode(imode),sim=sim,/gain,dname=dname,drift=drift,daisy=daisy,onlysim=onlysim
printf,5,smode(imode)
printf,5,'corr:',corrv1,corrv2
printf,5,'error:',corre1,corre2,corre3
printf,5,'zero:',corrm1,corrm2
printf,5,'autocorr:',autocorr
printf,5,''
printf,6,corrv1
printf,7,corre1
printf,8,autocorr
endfor

close, 5
close,6
close,7
close,8
;


end
