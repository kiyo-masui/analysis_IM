PRO plt_fg,corrv,corre,radio,weight,optical,field=field,fin=fin,radiocut=radiocut,fdir=fdir,noise=noise,sim=sim,daisy=daisy,outfile=outfile,dname=dname,svdfg=svdfg,mode=mode

loadct,0
set_plot,'ps'
device,filename='test.ps',/color, bits_per_pixel=8
tvscl,dist(300)
device,/close
;epsplot,'test.eps', xsize=16, ysize=16
;tvscl,dist(300)
;epsclose


seed=1
if not keyword_set(radiocut) then radiocut=0.1
if not keyword_set(fdir) then fdir=''
if not keyword_set(field) then field='f3'
if not keyword_set(dname) then dname='sept07'
if not keyword_set(mode) then mode=''

if field eq 'f3' then begin
   ;
   thisdir='/cita/scratch/cottontail/tchang/GBT/Processed/'
;   thisdir='/Users/tzu/Projects/GBT/pros/newsvd_day/'
   if dname eq 'aug29' then fdir=dname+'/f3_53-58/'
   if dname eq 'aug30' then fdir=dname+'/f3_45-50/'
   if dname eq 'sept07' then fdir=dname+'/f3_27-29/'
;   fdir=['aug29/f3_53-58/','aug30/f3_45-50/','sept07/f3_27-29/';,'sept17/f3_24-29/']
   ffdir=thisdir+fdir
   ;
endif

if field eq 'f4' then begin
   ;
   thisdir='/cita/scratch/cottontail/tchang/GBT/Processed/'
;   thisdir='/Users/tzu/Projects/GBT/pros/newsvd_day/'
   if dname eq 'aug29' then fdir=dname+'/f4_69-74/'
   if dname eq 'aug30' then fdir=dname+'/f4_67-72/'
   if dname eq 'sept06' then fdir=dname+'/f4_33-36/'
   if dname eq 'sept17' then fdir=dname+'/f4_48-53/'
   ffdir=thisdir+fdir
   ;
endif
set_plot,'ps'

; cross correlate GBT with the deep2 fields

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq0=1420.405751 ;MHz
beam_sigma=beam/2.35482


deepdir='/cita/scratch/cottontail/tchang/GBT/Deep2/'
;deepdir='/Users/tzu/Projects/GBT/Deep2/'
; field center
if (field eq 'f1') then begin
 ;  f0=[214.250,52.50]
   f0=[96.47,59.75]
   fbin=[2,14]
   fn=deepdir+'deep2_f1.dat'
   readcol,fn,ratmp,dectmp,oz,magi,format='d,d,d,d'
   euler,ratmp,dectmp,ora,odec,1
endif 
if (field eq 'f2') then begin
   f0=[252.43229,34.938242]     ;field center coordinate
   fsize=[1.43867,0.500210]    ;apparent size on the sky
   fbin=[12,4]
   fnopt=deepdir+'deep2_f2.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif
if (field eq 'f3') then begin
   f0=[352.49077,0.13736937]
   fsize=[2.24112, 0.538381]
;   fbin=[18,5]
   fbin=[18,4]
   fnopt=deepdir+'deep2_f3.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[12,4]
   fnopt=deepdir+'deep2_f4.dat'
   readcol,fnopt,ora,odec,oz,magi,format='d,d,d,d'
endif


npol=4
nx=fbin(0)
ny=fbin(1)
;nday=n_elements(ffdir)
;
radio=dblarr(npol,nx,ny,nz)
weight=radio
Da_z,oz,dist
distbin=(floor((dist-d0)/dd)) ;d0 is at bin0
;
; read in data
;   tmp=dblarr(npol,nx,ny,nz)
   openr,1,ffdir+field+'.radiotot_weighted_b4svd'
   readf,1,radio
   close,1
;   radioall(i,*,*,*)=tmp
   ;
   openr,1,ffdir+field+'.weighttot_b4svd'
   readf,1,weight
   close,1
;   weightall(i,*,*,*)=tmp
;
;
; apply calibratoin
; reform the arrays
radio=reform(radio,npol,nx*ny,nz)
weight=reform(weight,npol,nx*ny,nz)
;noiseradio=reform(noiseradio,Npol,nx*ny,nz)
help,weight
help,radio
;
;device,filename='fg.ps'
;
; read in gain_x and gain_z
gainx=dblarr(npol,nx*ny)
gainz=dblarr(npol,nz)
;
gdir='~/projects/GBT/pros/'+dname+'/'+field+'/'
;gdir='/Users/tzu/Projects/GBT/pros/'+dname+'/'
openr,1,gdir+'deep2_gain_x_b4svd.txt'
readf,1,gainx
close,1
;
openr,1,gdir+'deep2_gain_z_b4svd.txt'
readf,1,gainz
close,1
;
radio=reform(radio,npol,nx,ny,nz)
weight=reform(weight,npol,nx,ny,nz)
help,radio
help,weight

;for field 3, chop off the top dec strip and some noisy strips
;radioall=radio
;weightall=weight
;
; reform the gain too
gainz=gainz(*,2:nz-4)
nzz=555
gainz=gainz(*,nzz/3:nzz-1)

; get rid of the empty dec strips
;ny=3
;radio=radio(*,*,0:ny-1,*)
;weight=weight(*,*,0:ny-1,*)
if field eq 'f3' then begin
;
radio(*,*,3,*)=0
weight(*,*,3,*)=0
;
; get rid of the ix=13 for every dec strip
;nx=13
;radio=radio(*,0:nx-1,*,*)
;weight=weight(*,0:nx-1,*,*)
radio(*,13:nx-1,*,*)=0
weight(*,13:nx-1,*,*)=0
;
endif

if field eq 'f4' then begin
;
radio(*,*,0,*)=0
weight(*,*,0,*)=0
;
endif
;
;
;
radioall=radio
weightall=weight

device,filename='fg.ps',/color, bits_per_pixel=8
device,xsize=16, ysize=12
;epsplot,'fg.eps';, xsize=16, ysize=16

for ipol=0,Npol-1 do begin

;if keyword_set(nosim) then begin
   radio=reform(radioall(ipol,*,*,*))
   weight=reform(weightall(ipol,*,*,*))
   print,'max(radio),min(radio),mean(radio),variance(radio)'
   print,max(radio),min(radio),mean(radio),variance(radio)
   ind=where(weight gt 0,cind)

   ; do a svd on the spatial-redsfhit space
   nz=560
   zstart=2
   radio=radio(*,*,zstart:nz-4)
   weight=weight(*,*,zstart:nz-4)
   nz=555
   nnz=3
   if (field ne 'f1' and field ne '3c286') then begin
      help,radio
      help,weight
      radio=radio(*,*,nz/3:nz-1)
;        ;help,radio
      weight=weight(*,*,nz/3:nz-1)
      if ipol eq 0 then begin
         inddist=where(distbin ge (zstart+nz/3),count)
         print,'sources with z>0.7',count
         distbin=distbin(inddist)
         ora=ora(inddist)
         odec=odec(inddist)
         distbin=distbin-nz/3-zstart
      endif
      nz=(nz/3)*2
      nnz=2
   endif else begin
      if ipol eq 0 then begin
         inddist=where(distbin ge zstart,count)
         print,'sources with z>zstart',count
         distbin=distbin(inddist)
         ora=ora(inddist)
         odec=odec(inddist)
         distbin=distbin-zstart
      endif
   endelse
   ;
   if ipol eq 0 then begin
      radioadd=dblarr(nx*ny,nz)
      weightadd=dblarr(nx*ny,nz)
   endif
   ;
   print,'mean HI emission:',mean(radio),min(radio),max(radio)
     ;
;   plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(radio*weight,nx*ny,nz)),/scalable,/colbar
;   plt_image,transpose(reform(weight,nx*ny,nz)),/scalable,/colbar
   ;
   radio=reform(radio,nx*ny,nz)
   weight=reform(weight,nx*ny,nz)
   ;
   ; do -1
   ind=where(weight gt 0,cind)
   radio(ind)=radio(ind)-1.
     ;
     ;
     ; set back to kelvin
   for ix=0,nx*ny-1 do begin
      for iz=0,nz-1 do begin
         if (gainx(ipol,ix)*gainz(ipol,iz) ne 0.) then begin
            radio(ix,iz)=radio(ix,iz)/gainx(ipol,ix)/gainz(ipol,iz) 
         endif else begin
            weight(ix,iz)=0.
         endelse
      endfor
   endfor
   ;
   ;
;   tvscl,transpose(reform(radio,nx*ny,nz))
   plt_image,transpose(reform(radio,nx*ny,nz)),/scalable,/colbar
   ;
   radioadd=radioadd+radio*weight
   weightadd=weightadd+weight
   ;
endfor


; change the dimension of radio and do proper weight
; do weighted average?
ind=where(weightadd gt 0,cind)
radioadd(ind)=radioadd(ind)/weightadd(ind)
;
;
;plt_image,transpose(reform(radioadd,nx*ny,nz)),/scalable,/colbar
;plt_image,transpose(reform(radioadd*weightadd,nx*ny,nz)),/scalable,/colbar
;plt_image,transpose(reform(weightadd,nx*ny,nz)),/scalable,/colbar
;
;read,test
device,/close
;
;epsclose

END

      
      
