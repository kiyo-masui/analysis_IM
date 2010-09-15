pro plt_radio,radio,weight,fn=fn,field=field

set_plot,'ps'
device,/color, bits_per_pixel=8
device,filename='radio.ps'
;device,xsize=18.,ysize=6.
dir='/cita/d/raid-cita/tchang/GBT/zCOSMOS/Processed/18-26/'
;dir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_100muK/sept06/f4_33-36/'
;fname='f3.radio.processed.combined10modesvd3'
;fwname='f3.weight.processed.combined10modesvd3'
fname='f3.radio.processed.combined0modesvd'
fwname='f3.weight.processed.combined0modesvd'
;fname='f4.radio.processed.combined0modesvd'
;fwname='f4.weight.processed.combined0modesvd'
;fname='f4.radiotot0mode_weighted_aftsvd'
;fwname='f4.weighttot0mode_aftsvd'

nx=9
ny=9
nz=560

radio=dblarr(nx,ny,nz)
weight=radio
openr,1,dir+fname
readf,1,radio
close,1
openr,1,dir+fwname
readf,1,weight
close,1
;
print,'total weight',total(weight)
;radio=radio>(-10^(-4.))<10^(-4.)
print,max(radio),min(radio),mean(radio),median(radio)
;
radio=reform(radio,nx*ny,nz)
weight=reform(weight,nx*ny,nz)
neweight=dblarr(nx*ny,nz)
weightx=dblarr(nx*ny,nz)
;
for iz=0,nz-1 do begin
   ind=where(weight(*,iz) gt 0,cind)
   if cind gt 1 then neweight(*,iz)=(1./variance(radio(ind,iz)))
endfor
;
;ind=where(neweight gt 10^10.,cind)
;if cind gt 0 then neweight(ind)=0.

print,max(neweight),min(neweight),mean(neweight)

for ix=0,nx*ny-1 do begin
   weightx(ix,*)=total(weight(ix,*))
endfor
;
;radio=radio;*neweight*weightx
radio=radio*neweight*weightx*30.
;ind=where(abs(radio) gt 10^12.2,cind)
;if cind gt 0 then radio(ind)=0.
weight=neweight*weightx
dist=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   dist(*,*,iz)=1450.+(2+185.+iz)*2.
endfor
meandist=total(dist*weight)/total(weight)
print,'meandist',meandist
; effective redshift considering the weighting is z=0.815
; 
ind=where(radio ne 0,cind)
var=sqrt(total((radio)^2.)/total(weight^2.))
print,'var:',var
print,'rms:',sqrt(var)
help,radio

;radio1=reform(radio(0:nx-1,*))
;radio2=reform(radio((ny-1)*nx:nx*ny-1,*))
;radio=[radio1,radio2]
;
;weight1=reform(weight(0:nx-1,*))
;weight2=reform(weight((ny-1)*nx:nx*ny-1,*))
;radio=radio1
;weight=[weight1,weight2]
;var=sqrt(total(radio^2.)/total(weight^2.))
;ind=where(radio ne 0,cind)
;var=sqrt(variance(radio(ind)))
;print,'var:',var
;var=sqrt(total(radio(ind)^2.)/double(cind))
;print,'var:',var

; for 10modesvd3:  var=0.0037722016
; for 0 mode:  var=0.081460268 with weight
;for 0 mode:  var=0.12808851

;plt_image,transpose(radio),colbar=[-0.06,0.06],/scalable,frame=[1800,2600,0,360],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'
;plt_image,transpose(radio*weight),/colbar,/scalable,frame=[1800,2600,0,360],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'
;plt_image,transpose(radio*neweight),/colbar,/scalable,frame=[1800,2600,0,360],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'
;plt_image,transpose(radio),/colbar,/scalable,frame=[1800,2600,0,120],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'
plt_image,transpose(radio),/colbar,/scalable,frame=[1400,2600,0,40*9],xtitle='!6 Redshift Distance [Mpc]',ytitle='Spatial Distance [Mpc]'

device,/close
!p.multi=0

end
