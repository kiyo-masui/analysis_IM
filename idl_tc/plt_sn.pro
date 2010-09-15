pro plt_sn,daisy=daisy,drift=drift

if keyword_set(daisy) then obsmode='daisy' else obsmode='drift'

svd=[0,1,2,3,4,5,6,10,15,20,25]
field='zcosmos'
;
ns=n_elements(svd)
;
sig=dblarr(ns)
err=sig
sim=sig
simonly=sig
svdsim=dblarr(ns)
;
fn=field+'.tot.signal.'+obsmode+'.dat'
fnerr=field+'.tot.err.'+obsmode+'.dat'
fnsim=field+'.sim.tot.signal.'+obsmode+'.dat'
fnsimonly=field+'.tot.signal.'+obsmode+'.simonly.dat'
;fnsvd=field+'.residual.tot'+sigcut+'.signal.dat' 
;
openr,1,fn
readf,1,sig
close,1
;
openr,1,fnerr
readf,1,err
close,1
;
openr,1,fnsim
readf,1,sim
close,1
;
openr,1,fnsimonly
readf,1,simonly
close,1
;
;openr,1,fnsvd
;readf,1,svdsim
;close,1
;
;f3=sig(1:ns)
;f3err=err(1:ns)
;f3sim=sim(1:ns)
;;f3svdsim=svdsim
;f3simonly=simonly(1:ns)

f3=sig
f3err=err
f3sim=sim
f3simonly=simonly



f3sn=f3/f3err

set_plot,'ps'
device,filename='projection.'+obsmode+'.ps'

plot,svd,(f3sim-f3)/f3simonly,psym=2,xtitle='!6 subtracted SVD modes '+obsmode,ytitle='Fractional residual signal',$
     xran=[0,26],/xstyle,yran=[0,2.5],/ystyle
;oplot,[13.5],[2],psym=2
xyouts, [14],[2],'[Corr(s+d)-Corr(d)]/Corr(s)'
xyouts,[14],[2.2],obsmode+' scan'
;oplot,(f3sim-f3),linestyle=2
;oplot,f3simonly,linestyle=1
;oplot,[13.5,14],[1.5,1.5],linestyle=2
;xyouts,[14],[1.5],'Corr(s+d)-Corr(d)'
;oplot,[13.5,14],[1.,1.],linestyle=1
;xyouts,[14],[1.],'Corr(s)'
device,/close

set_plot,'ps'
device,filename='snratio.'+obsmode+'.ps'

plot,svd,f3sn,psym=2,xtitle='!6 subtracted SVD modes '+obsmode,ytitle='SN ratio',$
     yran=[-3,10],/ystyle,xran=[0,26],/xstyle
oploterror,svd,f3sn,f3err/max(f3err)*4,/nohat
;oplot,svd,f4sn,psym=5
;oploterror,svd,f4sn,f4err/max(f4err)*4,/nohat

device,/close

device,filename='correlation.'+obsmode+'.ps'

print,'f3 correction:',(f3simonly/(f3sim-f3))*(100e-6/f3simonly)
;
;
f3corr=f3*(f3simonly/(f3sim-f3))*(100e-6/f3simonly)
f3err=f3err*(100e-6/(f3sim-f3))
corr=(f3corr)
err=f3err
;
plot,svd,f3corr,psym=2,xtitle='!6 subtracted SVD modes '+obsmode,$
     ytitle='Correlation [K]',$
     /ystyle,xran=[-1,26],/xstyle,yran=[-2e-4,12e-4]
oploterror,svd,f3corr,f3err,/nohat
;
oplot,svd,(f3sim-f3),linestyle=2
oplot,svd,f3simonly,linestyle=1
oplot,[14],[10e-4],psym=2
xyouts,[15],[10e-4],'corrected Corr(d)'
oplot,[13.5,15],[9e-4,9e-4],linestyle=2
xyouts,[15],[9e-4],'Corr(s+d)-Corr(d)'
oplot,[13.5,15],[8e-4,8e-4],linestyle=1
xyouts,[15],[8e-4],'Corr(s)'
xyouts,[13.5],[11e-4],obsmode+' scan'

print,'corr',f3corr
print,'err',f3err
print,'sn',f3corr/f3err

device,/close

;device,filename='correlation_svd'+sigcut+'.ps'
;plot,svd,corr,psym=2,xtitle='!6 subtracted SVD modes',$
;     ytitle='Correlation [K]',$
;     /ystyle,xran=[0,26],/xstyle,yran=[-1e-4,6e-4]
;oploterror,svd,corr,err,/nohat

;device,/close

end
