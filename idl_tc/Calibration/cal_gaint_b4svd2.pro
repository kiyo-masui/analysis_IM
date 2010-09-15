pro cal_gaint_b4svd2,dname=dname,field=field

; calculate T_cal(nu) as a function of the binned (x,z) matrix
; -- basically take T_cal(nu) and convert it into T_cal(z) 

if not keyword_set(field) then field='zcosmos'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='03'

nx=9
ny=9
;
if dname eq '00' then begin
   fcaldir=['5-8']
   ffdir=['9-17', '18-26'] 
endif
;
if dname eq '02' then begin
   fcaldir=['31-34']
   ffdir=['scan9-17', 'scan18-26','scan27-28']
endif
;
if dname eq '03' then begin
   fcaldir=['5-8']
   ffdir=['62-70']
endif


; parameters
nz=560
npol=4
nfile=n_elements(ffdir)

;read in noise
ffcaldir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed/'+dname+'/'+ffdir(nfile-1)+'/'
noise=dblarr(npol,nx,ny,nz)
openr,1,ffcaldir+field+'.noisetot_weighted_aftsvd'
readf,1,noise
close,1
;
weight=dblarr(npol,nx,ny,nz)
openr,1,ffcaldir+field+'.weighttot_aftsvd'
readf,1,weight
close,1
;
;
; read in tcal
tdir='/cita/h/home-1/tchang/projects/GBT/pros/cosmos/Calibration/'+dname+'/'
tcal=dblarr(2,nz)
openr,1,tdir+'tcal_z.txt'
readf,1,tcal
close,1

; first change dimensions
noise2d=reform(noise,npol,nx*ny,nz)
weight2d=reform(weight,npol,nx*ny,nz)

; first calculate g(z)
; by taking Pcal(z) = g(z) Tcal(z)
; so first calculate Pcal(z)--take time median of noisecal
pcal_z=dblarr(npol,nz)
for i=0,npol-1 do begin
   for j=0,nz-1 do begin
      ind=where(weight2d(i,*,j) gt 0,cind)
      if (cind gt 0) then pcal_z(i,j)=median(noise2d(i,ind,j),/double)
   endfor
endfor
;
; now calculate g(z)=pcal_z/tcal_z
gain_z=dblarr(npol,nz)
for i=0,npol-1 do begin
   ind=where(tcal((i mod 2),*) ne 0.,cind)
   if (cind gt 0) then gain_z(i,ind)=pcal_z(i,ind)/tcal((i mod 2),ind)
endfor

openw,1,tdir+'/zcosmos_gain_z_b4svd.txt'
printf,1,gain_z
close,1

set_plot,'ps'
device,filename=tdir+'/zcosmos_gain_z_b4svd.ps'
plot,reform(gain_z(0,*)),xtitle='!6redshift bins',ytitle='gain(z)'
oplot,reform(gain_z(1,*)),linestyle=2
plot,reform(gain_z(2,*)),xtitle='!6redshift bins',ytitle='gain(z)'
oplot,reform(gain_z(3,*)),linestyle=2
device,/close

; calculate gain(t)
; calculate g_(t) = Int_dz [pcal_z * Pcal(z,t)] / Int_dz [Pcal_z^2]
gain_t=dblarr(npol,nx*ny)
for i=0,npol-1 do begin
   for j=0,nx*ny-1 do begin
      ind=where(weight2d(i,j,*) gt 0,cind)
      if (cind gt 0) then gain_t(i,j)=total(pcal_z(i,*)*noise2d(i,j,*))/total(pcal_z(i,*)^2.)
   endfor
endfor

help,gain_t

openw,1,tdir+'/zcosmos_gain_x_b4svd.txt'
printf,1,gain_t
close,1

set_plot,'ps'
device,filename=tdir+'/zcosmos_gain_x_b4svd.ps'
plot,reform(gain_t(0,*)),xtitle='!6spatial bins',ytitle='gain(x)'
oplot,reform(gain_t(1,*)),linestyle=2
plot,reform(gain_t(2,*)),xtitle='!6spatial bins',ytitle='gain(x)'
oplot,reform(gain_t(3,*)),linestyle=2
device,/close


end
