pro docalib_svdfill,field=field,ffdir=ffdir,dname=dname

; field:  data field
; fn:  fits file name

;if not keyword_set(field) then field='f3'  ; for the radio/optical cube setup
;if not keyword_set(fn) then fn='/home/scratch/kbandura/friday/field3_daisy_266-270.raw.acs.fits'

if not keyword_set(field) then field='3c218'  ; for the radio/optical cube setup
if not keyword_set(dname) then dname='04'

; more efficient to put the loop over file outside the loop this funciton
; since there are now multiple filed in the same session.

if not keyword_set(ffdir) then begin
  if dname eq '00' then begin
     ffdir=['5-8'];,'34-37']
  endif

  if dname eq '02' then begin
     ffdir=['5-8']
  endif

  if dname eq '03' then begin
     ffdir=['5-8']
  endif

  if dname eq '04' then begin
     ffdir=['5-8']
  endif
endif

outdirroot = '/mnt/raid-project/gmrt/kiyo/wiggleZ/calibration/'
; indirroot = '/cita/d/raid-cita/tchang/wiggleZ/GBT10B_036/'
indirroot = '/mnt/raid-cita/tchang/wiggleZ/GBT10B_036/'

set_plot,'ps'
!p.multi=[0,8,1]

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq_rest=1420.405751 ;MHz
nfreq=2048
Npol=4

; field center
if (field eq '3c286') then begin
   f0= [202.78453, 30.509155] 
endif
;3c67       02:24:12.2973    +27:50:12.016
if (field eq '3c324') then begin
   f0 = [237.45427,21.427460]
endif
if (field eq '3c48') then begin
   f0 = [24.427248,33.157969]
endif
if (field eq '3c348') then begin
   f0 = [252.783945,    4.992588]
endif
; 3c 218  09h18m05.7s -12d05m44s [J2000]
if (field eq '3c218') then begin
   f0= [139.52375, -12.095556]
endif
;

nfile=n_elements(ffdir)
;nn=22976L  ;359*8*
;nn=7680  ;120*8*8
nn=3840 ;60*8*8.
nscan=4
ntime=60
npol=4
nwin=8
nfreq=2048

windowb=[0,62,125,190,257,327,400,476]
windowe=[61,124,189,256,326,399,475,559]


tcaltsys=dblarr(nwin*nscan*nfile,npol)
itcal=0

for ifile=0,nfile-1 do begin
;for ifile=0,0 do begin

   fdir = outdirroot+dname+'/'+field+'_'+ffdir(ifile)+'/'
   spawn,'mkdir -p '+fdir
   fn = indirroot+dname+'_'+field+'_onoff_'+ffdir(ifile)+'.raw.acs.fits'
   
   if dname eq '00' and ffdir(ifile) eq '34-37' then begin
      nn=240*8*8
      ntime=240
   endif


   for scan=0L,nscan-1L do begin

      print,'******'
      print,'** processing scan ***',scan
      print,'******'
      
      scann=strtrim(string(scan),2)
      rec0=nn*scan
      rec1=nn*(scan+1)-1L
   
;   rec0=0
;   rec1=nn
      print,'fdir: ',fdir
      print,'fn: ',fn
      print,'rec0',rec0
      print,'rec1',rec1
      zmax_rec=nz
      zmin_rec=nz
      flagnumon=0L
      flagnumoff=0L
      
      !p.multi=0
  ; device,filename='Tnoise'+scann+'.ps'

      average_freq=dblarr(ntime,npol*2,8*nwin)

      time_freq_tot=dblarr(ntime,npol,nfreq*nwin)
      freqkeep=dblarr(nfreq*nwin)
      time_freq_matrix=dblarr(ntime,npol,nfreq*3/5*nwin)
      time_freq_crosspol=time_freq_matrix
      flagtot=time_freq_matrix
;   flagtot=make_array(ntime,npol,nfreq*3/5*nwin,/double,value=1.)
   ; collect the frequencies,ra,dec as well
      freqtot=dblarr(nfreq*3/5*nwin)
      ratot=dblarr(ntime*nwin)
      dectot=ratot
   ;
      tcaltsys=dblarr(nwin,npol)
      tcal=dblarr(nwin,npol)
      tsys=dblarr(nwin,npol)

      for i=0,nwin-1 do begin
         
         device,filename=fdir+'bandpass'+strtrim(string(scan),2)+strtrim(string(i),2)+'.ps'
         
         calibration_svd,t_off,t_on,noisecal,noise_time,noise_bpass,flag_off,flag_on,ra,dec,freq,pol,parallatic,recn=[rec0,rec1],thisif=i,threshould=0.02,fn=fn,/drift  

         openw,1,fdir+'TcalTsys.scan'+strtrim(string(scan),2)+strtrim(string(i),2)+'.txt'
         printf,1,noisecal
         close,1

         openw,1,fdir+'Pcal.scan'+strtrim(string(scan),2)+strtrim(string(i),2)+'.txt'
         printf,1,noise_time
         close,1

         openw,1,fdir+'Freq.scan'+strtrim(string(scan),2)+strtrim(string(i),2)+'.txt'
         printf,1,freq
         close,1

         
         for ipol=0,npol-1 do begin
            tcaltsys(i,ipol)=median(noisecal(ipol,nfreq/5:nfreq/5*3))
            tcal(i,ipol)=median(noise_time(0,ipol,nfreq/5:nfreq/5*3))
            tsys(i,ipol)=median(noise_bpass(0,ipol,nfreq/5:nfreq/5*3))
         endfor

      endfor

      openw,1,fdir+'TcalTsys.scan'+strtrim(string(scan),2)+'.freqmed.txt'
      printf,1,tcaltsys
      printf,1,''
      printf,1,tcal
      printf,1,''
      printf,1,tsys
      close,1

   endfor

endfor

end



