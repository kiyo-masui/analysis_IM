pro plt_auto_sn,sigcut=sigcut

if keyword_set(sigcut) then sigcut='.5sigcut' else sigcut=''
if keyword_set(sigcut) then sigcut5='(5-sigma cut)' else sigcut5=''
;
field2=['field3','field4'] 
svd=[0,1,2,3,4,5,6,8,10,12]
;svd=[1,2,3,4,5,6,7,8,9,10,11,12,15,20,25]
;
ns=n_elements(svd)
;
sig=dblarr(ns)
auto=sig
err=sig
sim=sig
data=sig
derr=sig
;
for i=0,1 do begin
;
field=field2(i)
; 
fn=field+'auto.crossday.10mode+svd.signal.dat'
fnerr=field+'auto.crossday.10mode+svd.err.dat'
fnsim=field+'corr.sim.10mode+svd.signal.dat'
fndata=field+'corr.10mode+svd.signal.dat'
fnderr=field+'corr.10mode+svd.err.dat'
;fnsimonly=field+'.simonly.samew.tot'+sigcut+'.signal.dat'
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
openr,1,fndata
readf,1,data
close,1
;
openr,1,fnderr
readf,1,derr
close,1
;
;openr,1,fnsvd
;readf,1,svdsim
;close,1
;
if field eq 'field3' then begin
   f3=sig
   f3err=err
   f3sim=sim
   f3data=data
   f3derr=derr
endif 
;
if field eq 'field4' then begin
   f4=sig
   f4err=err
   f4sim=sim
   f4data=data
   f4derr=derr
endif
;
endfor


set_plot,'ps'
;
;print,'f3 correction:',(100e-6/(f3sim-f3data))
;print,'f4 correction:',(100e-6/(f4sim-f4data)) 
;
;
factor3=(100e-6/(f3sim-f3data))
factor4=(100e-6/(f4sim-f4data)) 
f3auto=f3*factor3
f4auto=f4*factor4

print,'f3 factor',factor3
print,'f4 factor',factor4
print,''
f3err=f3err*factor3
f4err=f4err*factor4
f3data=f3data*factor3
f4data=f4data*factor4
f3derr=f3derr*factor3
f4derr=f4derr*factor4
;
print,'f3 cross-corr',f3data
print,'f4 cross-corr',f4data
print,''
print,'f3 auto',f3auto
print,'f4 auto',f4auto
print,''
print,'f3 err',f3err
print,'f4 err',f4err
;
;corr=(f3corr/f3err^2.+f4corr/f4err^2.)/(1./f3err^2.+1./f4err^2.)
;err=sqrt(1./(1./f3err^2.+1./f4err^2.))
;
device,filename='auto-correlation.ps'
plot,svd,f3auto,psym=2,xtitle='!6 Spatial-Redshift space subtracted SVD modes',$
     ytitle='Correlation [K]',$
     /ystyle,xran=[-1,15],/xstyle,yran=[0.,1.5e-3]
oploterror,svd,f3auto,f3err,/nohat
oplot,svd,f3data,linestyle=6
oplot,(svd+0.1),f4auto,psym=5
oploterror,(svd+0.1),f4auto,f4err,/nohat
oplot,svd,f4data,linestyle=1
oplot,[9],[0.0013],psym=2
xyouts, [10],[0.0013],'F3 auto'
oplot,[9,9.5],[0.0012,0.0012],linestyle=6
xyouts, [10],[0.0012],'F3 Corr'
oplot,[9],[0.0011],psym=5
xyouts, [10],[0.0011],'F4 auto'
oplot,[9,9.5],[0.001,0.001],linestyle=1
xyouts, [10],[0.001],'F4 Corr'

device,/close

corr=(f3auto+f4auto)/2.
err=(f3err+f4err)/2.

device,filename='autocorr.ps'
plot,svd,corr,psym=2,xtitle='!6 Spatial-Redshift space subtracted SVD modes',$
     ytitle='Correlation [K]',$
     /ystyle,xran=[-1,15],/xstyle,yran=[0.,1.5e-3]
oploterror,svd,corr,err,/nohat
oplot,svd,f3data,linestyle=6
oplot,svd,f4data,linestyle=1
oplot,[9],[0.0013],psym=2
xyouts, [10],[0.0013],'Averaged auto'
oplot,[9,9.5],[0.0012,0.0012],linestyle=6
xyouts, [10],[0.0012],'F3 Corr'
oplot,[9,9.5],[0.0011,0.0011],linestyle=1
xyouts, [10],[0.0011],'F4 Corr'
device,/close

print,'averaged auto',corr
print,'averaged auto err',err

corr=(f3data/f3derr^2.+f4data/f4derr^2.)/(1./f3derr^2.+1./f4derr^2.)
err=sqrt(1./(1./f3derr^2.+1./f4derr^2.))

print,'corr:',corr
print,'corr err:',err

device,filename='correlation.3svd.ps'
plot,svd,corr,psym=2,xtitle='!6 subtracted SVD modes',$
     ytitle='Correlation [K]',$
     /ystyle,xran=[-1,15],/xstyle,yran=[-1e-4,3e-4]
oploterror,svd,corr,err,/nohat
device,/close

end
