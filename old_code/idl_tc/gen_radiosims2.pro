pro gen_radiosims2, signal=signal,dname=dname

field=['zcosmos']
mdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed/'
if not keyword_set(signal) then signal='100muK'
if signal eq '100muK' then hitemp=0.0001d
if signal eq '200muK' then hitemp=0.0002d
if signal eq '10muK' then hitemp=0.00001d
if not keyword_set(dname) then dname='02'

ffdir=['/02/30/'];,$
       ;'/02/scan27-28']
ffdir=mdir+ffdir
gdir=['~/projects/GBT/pros/cosmos/Calibration/'+dname+'/'];,'~/projects/GBT/pros/cosmos/Calibration/02/']
fbin=[9,9]
nz=560


nx=fbin(0)
ny=fbin(1)
nf=n_elements(ffdir)
npol=4
nzz=560
gainx=dblarr(npol,nx*ny)
gainz=dblarr(npol,nzz)

fdir='~/projects/GBT/zCOSMOS/zCOSMOS.opticaldensity.dat'
set_plot,'ps'
for ii=0,nf-1 do begin

   optical=dblarr(nx,ny,nz)
   ;
   openr,1,fdir
   readf,1,optical
   close,1
   ;
   ;expand optical
   opticaltot=dblarr(npol,nx,ny,nzz)
   for ipol=0,npol-1 do begin
      opticaltot(ipol,*,*,*)=optical
   endfor
   opticaltot=opticaltot*hitemp
   ;
   device,filename=gdir(ii)+'Toptical.ps',/color, bits_per_pixel=8
   opticaltmp=reform(opticaltot,npol,nx*ny,nzz)
   plt_image,transpose(reform(opticaltmp(0,*,*))),/scalable,/colbar
   device,/close
   ;
   ; change to power unit
   openr,1,gdir(ii)+'zcosmos_gain_x_b4svd.txt'
   readf,1,gainx
   close,1
   ;
   openr,1,gdir(ii)+'zcosmos_gain_z_b4svd.txt'
   readf,1,gainz
   close,1
   ;
   optical=reform(opticaltot,npol,nx*ny,nzz)
   for ipol=0,npol-1 do begin
      for ix=0,nx*ny-1 do begin
         for iz=0,nzz-1 do begin
            optical(ipol,ix,iz)=optical(ipol,ix,iz)*gainx(ipol,ix)*gainz(ipol,iz)
         endfor
      endfor
   endfor
   ;
   help,optical
   optical=reform(optical,npol,nx,ny,nzz)
   help,optical
   ;
   device,filename=gdir(ii)+'Poptical.ps',/color, bits_per_pixel=8
   opticaltmp=reform(optical,npol,nx*ny,nzz)
   plt_image,transpose(reform(opticaltmp(0,*,*))),/scalable,/colbar
   device,/close
   ;

;   openw,1,fdir+field+'.radiosim_'+signal
;   printf,1,optical
;   close,1
   ;
   openw,1,gdir(ii)+'radiosim_'+signal
   printf,1,optical
   close,1
endfor


end




