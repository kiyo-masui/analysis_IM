pro docalib_svdfill_gain_b4svd,radio,weight,t_off,radio_on,weight_on,field=field,fn=fn,sigmacut=sigmacut,dname=dname,fitsfile=fitsfile

; field:  data field
; fn:  fits file name

if not keyword_set(field) then field='f4'  ; for the radio/optical cube setup
if not keyword_set(fn) then fn='/cita/scratch/cottontail/tchang/GBT/sept07/d4_field3_drift_21-26.raw.acs.fits'
if not keyword_set(dname) then dname='sept07'

if dname eq 'sept07' then begin
;   ffdir=['05-08']
   ffdir=['09-14', '15-20', '21-26','27-29'] ;sept07
   dnum='d4'
endif
if dname eq 'aug29' then begin
   ffdir=['05-10', '15-20', '21-26', '27-32', '37-42', '47-52', '53-58'] ;aug29
   dnum='d1'
endif
if dname eq 'aug30' then begin
   ffdir=['9-14', '15-20', '21-26', '27-32','33-38', '39-44', '45-50'] ;aug30
   dnum='d2'
endif
if dname eq 'sept17' then begin
   ffdir=['12-17','18-23','24-29'];,'30-35','36-41','42-47','48-53'] ;sept 17
   dnum='d5'
endif
;
if field eq 'f4' then begin
;
if dname eq 'sept06' then begin
   ffdir=['09-14', '15-20', '21-26','27-32','33-36'] ;sept06
   dnum='d3'
endif
if dname eq 'aug29' then begin
   ffdir=['63-68','69-74'] ;aug29
   dnum='d1'
endif
if dname eq 'aug30' then begin
   ffdir=['55-60','61-66', '67-72'] ;aug30
   dnum='d2'
endif
if dname eq 'sept17' then begin
   ffdir=['30-35','36-41','42-47','48-53'] ;sept 17
   dnum='d5'
endif
;
endif
;
;
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
if (field eq 'nissim') then begin
   f0=[40.801455,-5.8477362]
   fbin=[12,4]
;   radio=dblarr(12,4,nz)
endif
if (field eq 'f1') then begin
;   f0=[214.250,52.50]
   ; use galactic coordiate
   f0=[96.47,59.75]
   fbin=[2,14]
endif 
if (field eq 'f2') then begin
   f0=[252.43229,34.938242]     ;field center coordinate
   fsize=[1.43867,0.500210]    ;apparent size on the sky
   fbin=[12,4]
endif
if (field eq 'f3') then begin
   f0=[352.49077,0.13736937]
   fsize=[2.24112, 0.538381]
   fbin=[18,4]
endif
if (field eq 'f4') then begin
   f0=[37.255199,0.58632926]
   fbin=[12,4]
endif
if (field eq '3c286') then begin
   f0= [202.78453, 30.509155] 
endif
;3c67       02:24:12.2973    +27:50:12.016
if (field eq '3c324') then begin
   f0 = [237.45427,21.427460]
endif
;
;
; parameters and matrices
nx=fbin(0)
ny=fbin(1)
radio=dblarr(Npol,nx,ny,nz)
radiotot=radio
weight=radio
weighttot=weight
noiseradio=dblarr(Npol,nx,ny,nz)
noisetot=noiseradio

; parameters
nfile=n_elements(ffdir)
nn=22976L  ;359*8*8
nscan=6
ntime=359
npol=4
nwin=8
nfreq=2048
if dname eq 'sept17' then begin
   if field eq 'f3' then begin
      nn=28736L
      ntime=449
   endif
endif
if field eq 'f4' then begin
   nn=19136
   ntime=299
endif
;
for ifile=0,nfile-1 do begin
   ;
   if (dname eq 'sept07' and ifile eq 3) then nscan=3
   if (dname eq 'sept06' and ifile eq 4) then nscan=4
   fdir='/cita/scratch/cottontail/tchang/GBT/Processed/'+dname+'/'+field+'_'+ffdir(ifile)+'/'
   spawn,'mkdir '+fdir
   ;fn='/cita/scratch/bear/pen/GBT/DATA/GBT08B_037/'+dname+'/'+dnum+'_field3_drift_'+ffdir(ifile)+'.raw.acs.fits'
   if field eq 'f3' then fn='/cita/scratch/bear/pen/GBT/DATA/GBT08B_037/'+dname+'/'+dnum+'_field3_drift_'+ffdir(ifile)+'.raw.acs.fits'
   if field eq 'f4' then fn='/cita/scratch/bear/pen/GBT/DATA/GBT08B_037/'+dname+'/'+dnum+'_field4_drift_'+ffdir(ifile)+'.raw.acs.fits'

   ;
   radio=dblarr(Npol,nx,ny,nz)
   weight=radio
   noiseradio=dblarr(Npol,nx,ny,nz)
   filenum=strtrim(string(ifile),2)
   ;
   for scan=0L,nscan-1L do begin
      ;
      print,'******'
      print,'** processing scan ***',scan
      print,'******'
      ;
      scann=strtrim(string(scan),2)
      rec0=nn*scan
      rec1=nn*(scan+1)-1L
   
      print,'fdir: ',fdir
      print,'fn: ',fn
      print,'rec0',rec0
      print,'rec1',rec1
      zmax_rec=nz
      zmin_rec=nz
      flagnumon=0L
      flagnumoff=0L

      !p.multi=0

      average_freq=dblarr(ntime,npol*2,8*nwin)
      time_freq_tot=dblarr(ntime,npol,nfreq*nwin)
      freqkeep=dblarr(nfreq*nwin)
      time_freq_matrix=dblarr(ntime,npol,nfreq*3/5*nwin)
      time_freq_crosspol=time_freq_matrix
      time_freq_noise=dblarr(ntime,npol,nfreq*3/5*nwin)
      flagtot=time_freq_matrix
      ; collect the frequencies,ra,dec as well
      freqtot=dblarr(nfreq*3/5*nwin)
      ratot=dblarr(ntime*nwin)
      dectot=ratot
      ra=dblarr(ntime)
      dec=dblarr(ntime)
      ;
      if keyword_set(fitsfile) then begin

      ;
      device,filename=fdir+'bandpass'+strtrim(string(scan),2)+'.ps'

      for i=0,nwin-1 do begin

         calibration_svd,t_off,t_on,noise,flag_off,flag_on,ra,dec,freq,pol,parallatic,recn=[rec0,rec1],thisif=i,threshould=0.02,fn=fn,/drift  

         print,'freq:',freq(0:10)
         print,'freq2:',freq(nfreq/2:nfreq/2+20)
         ; compare the overlap
         tbin=dblarr(ntime,npol,nfreq)
         tbin(*,0,*)=t_off(*,0,*)
         tbin(*,1,*)=t_off(*,3,*)
         tbin(*,2,*)=t_on(*,0,*)
         tbin(*,3,*)=t_on(*,3,*)
         ;
         tbinpol=dblarr(ntime,npol,nfreq)
         tbinpol(*,0,*)=t_off(*,1,*)
         tbinpol(*,1,*)=t_off(*,2,*)
         tbinpol(*,2,*)=t_on(*,1,*)
         tbinpol(*,3,*)=t_on(*,2,*)
         ;
         tbinflag=dblarr(ntime,npol,nfreq)
         tbinflag(*,0,*)=flag_off(*,0,*)
         tbinflag(*,1,*)=flag_off(*,3,*)
         tbinflag(*,2,*)=flag_on(*,0,*)
         tbinflag(*,3,*)=flag_on(*,3,*)
         ;
         ; rebin the noise cal source - keep only XX and YY
         tbinnoise=dblarr(ntime,npol,nfreq)
         tbinnoise(*,0,*)=noise(*,0,*)
         tbinnoise(*,1,*)=noise(*,3,*)
         tbinnoise(*,2,*)=noise(*,0,*)
         tbinnoise(*,3,*)=noise(*,3,*)
         ;
         ind=where(tbinflag eq 0,cind)
         print,'flagged:',cind,double(cind)/n_elements(tbinflag)
         ;
         ind=where(tbinflag(*,0,*)  eq 0,cind)
         print,'flagged:',cind,double(cind)/n_elements(tbinflag(*,0,*))
         print,'tbinflag off:',n_elements(tbinflag(*,0,*))
         ;
         ;
         ; calculate the bandpass (time median)
         bandpass=dblarr(npol,nfreq)
         ; do bandpass correction
         ; calculate the bandpass by taking time median per pol per freq
         for kk=0,nfreq-1 do begin
            for jj=0,npol-1 do begin
               ind=where(tbinflag(*,jj,kk) eq 1,cind)
               if (cind gt 1) then begin
                  band=tbin(ind,jj,kk)
                  fval=median(band,/double)
                  bandpass(jj,kk)=fval
               endif else begin
                  print,'no bandpass value:',jj,kk
                  bandpass(jj,kk)=1.
               endelse
            endfor
         endfor
         ;
         ; apply bandpass - divide data by its time-median
         for kk=0,nfreq-1 do begin 
            for jj=0,npol-1 do begin
               if bandpass(jj,kk) ne 0. then begin
                  tbin(*,jj,kk)=tbin(*,jj,kk)/bandpass(jj,kk)
                  ; divide the noise source by the same bandpass
                  tbinnoise(*,jj,kk)=tbinnoise(*,jj,kk)/bandpass(jj,kk)
               endif else begin
                  tbin(*,jj,kk)=0.
                  tbinflag(*,jj,kk)=0.
                  ; put the same mask 
                  tbinnoise(*,jj,kk)=0.
               endelse
            endfor
            if bandpass(0,kk)*bandpass(1,kk) ne 0. then begin
               tbinpol(*,0,kk)=tbinpol(*,0,kk)/(sqrt(bandpass(0,kk)*bandpass(1,kk)))
               tbinpol(*,1,kk)=tbinpol(*,1,kk)/(sqrt(bandpass(0,kk)*bandpass(1,kk)))
            endif 
            if bandpass(2,kk)*bandpass(3,kk) ne 0. then begin
               tbinpol(*,2,kk)=tbinpol(*,2,kk)/(sqrt(bandpass(2,kk)*bandpass(3,kk)))
               tbinpol(*,3,kk)=tbinpol(*,3,kk)/(sqrt(bandpass(2,kk)*bandpass(3,kk)))
            endif 
         endfor
         ;
         if (i eq 0) then begin
            tbin0=tbin
            freq0=freq
            tbinflag0=tbinflag
         endif else begin
            ;
            ratio=dblarr(ntime,npol)
            for ii=0,ntime-1 do begin
               for jj=0,npol-1 do begin
                                ; cdelta_freq is negative, and it
                                ; starts from the lowest freq window
                                ; 695 MHz
                  tbin0slice=tbin0(ii,jj,(nfreq/10-1):(nfreq/10+nfreq/5-1))
                  ind0=where(tbinflag0(ii,jj,(nfreq/10-1):(nfreq/10+nfreq/5-1)) eq 1, cind0)
                  tbinslice=tbin(ii,jj,(nfreq*3/5+nfreq/10):(nfreq*3/5+nfreq/10+nfreq/5))  
                  ind=where(tbinflag0(ii,jj,(nfreq*3/5+nfreq/10):(nfreq*3/5+nfreq/10+nfreq/5)) eq 1, cind)
                  if ((cind0 gt 1) and (cind gt 1)) then $
                     ratio(ii,jj)=median(tbin0slice(ind0),/double)/median(tbinslice(ind),/double) else ratio(ii,jj)=1.d0
                  ;
                  if (abs(ratio(ii,jj)) gt 2. or ratio(ii,jj) eq 1.d0) then begin
                     print,'ratio:',ratio(ii,jj),ii,jj,i
                  endif
               endfor
            endfor
            ;
            ;
            ; now match the IF windows
            for ii=0,ntime-1 do begin
               for jj=0,npol-1 do begin
                  if finite(ratio(ii,jj)) then begin
                     tbin(ii,jj,*)=tbin(ii,jj,*)*ratio(ii,jj)
                     tbinnoise(ii,jj,*)=tbinnoise(ii,jj,*)*ratio(ii,jj)
                  endif
               endfor
            endfor
            ;
            ;
            tbin0=tbin
            freq0=freq
            tbinflag0=tbinflag
         endelse
         ; write tbin into a big matrix
         nkeep=nfreq*3/5
         time_freq_tot(*,*,i*nfreq:(i+1)*nfreq-1)=tbin
         freqkeep(i*nfreq:(i+1)*nfreq-1)=freq
         time_freq_matrix(*,*,i*nkeep:(i+1)*nkeep-1)=tbin(*,*,nfreq/5+1:nfreq*4/5-1)
         time_freq_crosspol(*,*,i*nkeep:(i+1)*nkeep-1)=tbinpol(*,*,nfreq/5+1:nfreq*4/5-1)
         time_freq_noise(*,*,i*nkeep:(i+1)*nkeep-1)=tbinnoise(*,*,nfreq/5+1:nfreq*4/5-1)
         flagtot(*,*,i*nkeep:(i+1)*nkeep-1)=tbinflag(*,*,nfreq/5+1:nfreq*4/5-1)
         freqtot(i*nkeep:(i+1)*nkeep-1)=freq(nfreq/5+1:nfreq*4/5-1)
         ratot(i*ntime:(i+1)*ntime-1)=ra
         dectot(i*ntime:(i+1)*ntime-1)=dec
         ;
      endfor
      ;
      device,/close
      ;
      device,filename=fdir+'Tsky_spec'+strtrim(string(scan),2)+'.ps'
      plt_image,reform(time_freq_matrix(*,0,*)),/scalable,/colbar
      plt_image,reform(time_freq_matrix(*,1,*)),/scalable,/colbar
      plt_image,reform(time_freq_matrix(*,2,*)),/scalable,/colbar
      plt_image,reform(time_freq_matrix(*,3,*)),/scalable,/colbar
      ;
      ; now do flagging before svd, using cross-polariztion flags
      ; the flagged one is stored in time_freq_flag
      time_freq_flag=time_freq_matrix
      ;
      flagnumon=0L
      flagnumoff=0L
      dnu=50e6/double(nfreq)
      dt=0.5d
      sigma=2./(dnu*dt)         ;where did the 2 come from? ans: 2=1+1: pol^2+pol^2
      sigth=8.                  ; sigma threshould
      for i=0,ntime-1 do begin
      ; noise-off
         r= (time_freq_crosspol(i,0,*)^2+time_freq_crosspol(i,1,*)^2.) / $
            (time_freq_matrix(i,0,*)*time_freq_matrix(i,1,*))
         ;
         ind=where(r gt sigma*(sigth^2),cind)
         if (cind gt 0) then begin
            for j=0,1 do begin
               time_freq_flag(i,j,ind)=dblarr(cind)
               flagtot(i,j,ind)=dblarr(cind)
               time_freq_noise(i,j,ind)=0.
            endfor
            flagnumoff=flagnumoff+cind
         endif
         ;
      ; noise-on
         r= (time_freq_crosspol(i,2,*)^2+time_freq_crosspol(i,3,*)^2.) / $
            (time_freq_matrix(i,2,*)*time_freq_matrix(i,3,*))
         ;
         ind=where(r gt sigma*(sigth^2),cind)
         ;
         if (cind gt 0) then begin
            for j=2,3 do begin
               flagtot(i,j,ind)=dblarr(cind)
               time_freq_flag(i,j,ind)=dblarr(cind)
               time_freq_noise(i,j,ind)=0.
            endfor
            flagnumon=flagnumon+cind
         endif
      endfor
      ;
      plt_image,transpose(reform(time_freq_flag(*,0,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_flag(*,1,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_flag(*,2,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_flag(*,3,*))),/scalable,/colbar
      ;
      device,/close
      ;
      ; write outputs here, per scan, per file, per day
      openw,1,fdir+'radio_tf_b4svd'+scann
      printf,1,time_freq_flag
      close,1
      ;
      openw,1,fdir+'flagtot_tf_b4svd'+scann
      printf,1,flagtot
      close,1
      ;
      openw,1,fdir+'noisecal_tf_b4svd'+scann
      printf,1,time_freq_noise
      close,1
      ;
      if (scan eq 0) then begin
         ;
         openw,1,fdir+'freqtot'
         printf,1,freqtot
         close,1
         ;
         openw,1,fdir+'ra'
         printf,1,ra
         close,1
         ;
         openw,1,fdir+'dec'
         printf,1,dec
         close,1
         ;
      endif
      ;
      ;
      device,filename=fdir+'Tnoise'+strtrim(string(scan),2)+'.ps'
      plt_image,transpose(reform(time_freq_noise(*,0,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,1,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,2,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,3,*))),/scalable,/colbar
      device,/close
      ;
      print,'flagged points off',flagnumoff
      print,'flagged percentage off',double(flagnumoff)/double(n_elements(time_freq_matrix(*,0,*)))
      print,'flagged points on',flagnumon
      print,'flagged percentage on',double(flagnumon)/double(n_elements(time_freq_matrix(*,0,*)))
      ;
      ;
   endif else begin
      
            ; read in the data files
      openr,1,fdir+'radio_tf_b4svd'+scann
      readf,1,time_freq_flag
      close,1
      ;
      openr,1,fdir+'flagtot_tf_b4svd'+scann
      readf,1,flagtot
      close,1
      ;
      openr,1,fdir+'noisecal_tf_b4svd'+scann
      readf,1,time_freq_noise
      close,1
      ;
      openr,1,fdir+'freqtot'+scann
      readf,1,freqtot
      close,1
      ;
      openr,1,fdir+'ra'+scann
      readf,1,ra
      close,1
      ;
      openr,1,fdir+'dec'+scann
      readf,1,dec
      close,1
      ;   
   endelse
      ;

      if keyword_set(svdfg) then begin
         ;
         ; svd
         ; do each pol separately
         ;
         radiosvd=dblarr(ntime,npol,nfreq*3/5*npol*2)
         radiores=radiosvd
         radiosvd1mode=radiosvd
         radiores1mode=radiosvd
         radiosvd2mode=radiosvd
         radiores2mode=radiosvd
         radio2dtot=radiosvd
         radiosvd5mode=radiosvd
         radiores5mode=radiosvd
         radiosvd10mode=radiosvd
         radiores10mode=radiosvd
         ;
         for ii=0,npol-1 do begin
            ;
            radio2d=reform(time_freq_flag(*,ii,*))
            ;
            la_svd,radio2d, s,u,v,/double,status=status
            print,'SVD status:',status
            ;
            help,s
            help,u
            help,v
            ;
            w=dblarr(n_elements(s))
            w2=w
            w5=w
            w10=w
            w(0)=s(0)
            w2(0:1)=s(0:1)
            w5(0:4)=s(0:4)
            w10(0:9)=s(0:9)
            ;
            radiosvd(*,ii,*) = u ## diag_matrix(s) ## transpose(v)
            radiosvd1mode(*,ii,*) =u ## diag_matrix(w) ## transpose(v)
            radiosvd2mode(*,ii,*) = u ## diag_matrix(w2) ## transpose(v)
            radiosvd5mode(*,ii,*) =u ## diag_matrix(w5) ## transpose(v)
            radiosvd10mode(*,ii,*) = u ## diag_matrix(w10) ## transpose(v)
            ;
         endfor
         ;
         radiores = time_freq_flag-radiosvd
         radiores1mode = time_freq_flag-radiosvd1mode
         radiores2mode = time_freq_flag-radiosvd2mode
         radiores5mode = time_freq_flag-radiosvd5mode
         radiores10mode = time_freq_flag-radiosvd10mode
         ;
         device,filename=fdir+'Tsky_svd'+strtrim(string(scan),2)+'.ps'
         vfil=v
         vfil(15:*,*)=0
         plt_image, transpose(reform(time_freq_flag(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd1mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores1mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd2mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores2mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd5mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores5mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd10mode(*,2,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores10mode(*,2,*))),/scalable,/colbar
         ;
         ; do sigma flagging
         dnu=50e6/double(nfreq)
         dt=0.5d
         sigma=1./sqrt(dnu*dt)
         sigma_th=5.
         ;
         ind=where(abs(radiores5mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then begin
            radiores5mode(ind)=0
            flagtot(ind)=0
            time_freq_noise(ind)=0.
         endif
         ;
         ind=where(abs(radiores10mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores10mode(ind)=0
         ind=where(abs(radiores1mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores1mode(ind)=0
         ind=where(abs(radiores2mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores2mode(ind)=0
         plt_image, transpose(reform(radiores1mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores2mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores5mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores10mode(*,0,*))),/scalable,/colbar
         ; fill in the svd values to the holes
         ind=where(flagtot eq 0, cind)
         time_freq_flag(ind)=radiosvd5mode(ind)
         ; do frequency weighting
         freqweight=dblarr(npol,nfreq*3/5*nwin)
         for ii=0,npol-1 do begin
            for kk=0,nfreq*3/5*nwin-1 do begin
               ind=where(flagtot(*,ii,kk) eq 1,cind)
               if cind gt 1 then begin
                  freqweight(ii,kk)=variance(radiores5mode(ind,ii,kk))
                  time_freq_flag(ind,ii,kk)=time_freq_flag(ind,ii,kk)/freqweight(ii,kk)
               endif
            endfor
         endfor
         ;
         plt_image, transpose(reform(time_freq_flag(*,0,*))),/scalable,/colbar
         ;
         for ii=0,npol-1 do begin
            ;
            radio2d=reform(time_freq_flag(*,ii,*))
            la_svd,radio2d, s,u,v,/double,status=status
            print,'SVD status:',status
            help,s
            help,u
            help,v
            w=dblarr(n_elements(s))
            w2=w
            w5=w
            w10=w
            w(0)=s(0)
            w2(0:1)=s(0:1)
            w5(0:4)=s(0:4)
            w10(0:9)=s(0:9)
            ;
            radiosvd(*,ii,*) = u ## diag_matrix(s) ## transpose(v)
            radiosvd1mode(*,ii,*) =u ## diag_matrix(w) ## transpose(v)
            radiosvd2mode(*,ii,*) = u ## diag_matrix(w2) ## transpose(v)
            radiosvd5mode(*,ii,*) =u ## diag_matrix(w5) ## transpose(v)
            radiosvd10mode(*,ii,*) = u ## diag_matrix(w10) ## transpose(v)
            ;
         endfor
         ;
         radiores = time_freq_flag-radiosvd
         radiores1mode = time_freq_flag-radiosvd1mode
         radiores2mode = time_freq_flag-radiosvd2mode
         radiores5mode = time_freq_flag-radiosvd5mode
         radiores10mode = time_freq_flag-radiosvd10mode
         ; de-freq weight
         for ii=0,npol-1 do begin
            for kk=0,nfreq*3/5*nwin-1 do begin
               ind=where(flagtot(*,ii,kk) eq 1,cind)
               if cind gt 1 then begin
                  radiores5mode(ind,ii,kk)=radiores5mode(ind,ii,kk)*freqweight(ii,kk)
                  radiores10mode(ind,ii,kk)=radiores10mode(ind,ii,kk)*freqweight(ii,kk)
                  radiores1mode(ind,ii,kk)=radiores1mode(ind,ii,kk)*freqweight(ii,kk)
                  radiores2mode(ind,ii,kk)=radiores2mode(ind,ii,kk)*freqweight(ii,kk)
               endif
            endfor
         endfor
         ; do sigma flagging
         dnu=50e6/double(nfreq)
         dt=0.5d
         sigma=1./sqrt(dnu*dt)
         sigma_th=5.
         ;
         ind=where(abs(radiores5mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then begin
            radiores5mode(ind)=0
            flagtot(ind)=0
            time_freq_noise(ind)=0.
         endif
         ;
         ind=where(abs(radiores10mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores10mode(ind)=0
         ind=where(abs(radiores1mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores1mode(ind)=0
         ind=where(abs(radiores2mode) gt (sigma_th*sigma),cind)
         if (cind gt 0) then radiores2mode(ind)=0
         ;
         plt_image,v,/scalable,/colbar
         plot,v(0,*)
         plot,v(*,0)
         plt_image,reform(v(0:10,*)),/scalable,/colbar
         plot,w(0:10)
         ;
         plt_image, transpose(reform(radiosvd1mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores1mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd2mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores2mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd5mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores5mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiosvd10mode(*,0,*))),/scalable,/colbar
         plt_image, transpose(reform(radiores10mode(*,0,*))),/scalable,/colbar
         ;
         device,/close
         ;
         device,filename=fdir+'Tnoise_flag'+strtrim(string(scan),2)+'.ps'
         plt_image,transpose(reform(time_freq_noise(*,0,*))),/scalable,/colbar
         plt_image,transpose(reform(time_freq_noise(*,1,*))),/scalable,/colbar
         plt_image,transpose(reform(time_freq_noise(*,2,*))),/scalable,/colbar
         plt_image,transpose(reform(time_freq_noise(*,3,*))),/scalable,/colbar
         device,/close
         ;
         ;
      endif
      ;
      ;
      ; rest frequency of the 8 IF's 
      ;
      ; 695, 725, 755, 785, 815, 845, 875, 905  in MHz
      ;
      ; each 50 MHz bandwidth
      ;
      Da_z,(freq_rest/freqtot-1.),distslice
      zslice=(floor((distslice-d0)/dd))
      print,'zslice range:',min(zslice),max(zslice)
      ;
      indtmp=where(flagtot(*,0,*) eq 0,cind)
      print,'flagtot x,off:',cind
      indtmp=where(flagtot eq 0, cind) 
      print,'flagtot:', cind
      print,'flagtot %',double(cind)/n_elements(flagtot)
      help,flagtot
      ;
      mk_hicube_noise_fg,time_freq_flag,flagtot,time_freq_noise,ra,dec,freqtot,pol,$
                         parallatic,radio=radio,weight=weight,$
                         noiseradio=noiseradio,field=field,flagnum=flagnumoff
;   mk_hicube_noise_fg,radiores5mode,flagtot,time_freq_noise,ra,dec,freqtot,pol,parallatic,radio=radio,weight=weight,$
;             noiseradio=noiseradio,field=field,flagnum=flagnumoff
;   mk_hicube_svd,time_freq_noise,flagnoise,ra,dec,freqtot,pol,parallatic,radio=noiseradio,weight=noiseweight,$
;             field=field,flagnum=flagnumoff,Npol=2

;   mk_hicube,t_on,flagtot_on,ra,dec,freq,pol,parallatic,radio=radio_on,weight=weight_on,$
;             field=field,kstart=kstart,kend=kend,flagnum=flagnumon
      ind=where(radio ne 0,cind)
      print,''
      print,'outside mk_hicube:  non-zero radio elements:',cind
      print,''
      ;
      radio2=radio
      ind=where(weight gt 0.,cind)
      radio2(ind)=radio2(ind)/weight(ind)
      ;
      device,filename=fdir+'mk_hicub_radio_b4svd'+strtrim(string(scan),2)+'.ps'
      radim1=reform(radio2(0,*,*,*))
      radim2=reform(radio2(1,*,*,*))
      radim3=reform(radio2(2,*,*,*))
      radim4=reform(radio2(3,*,*,*))
      radim1=reform(radim1,nx*ny,nz,/overwrite)
      radim2=reform(radim2,nx*ny,nz,/overwrite)
      radim3=reform(radim3,nx*ny,nz,/overwrite)
      radim4=reform(radim4,nx*ny,nz,/overwrite)
      ;
      plt_image, transpose(radim1),/scalable,/colbar
      plt_image, transpose(radim2),/scalable,/colbar
      plt_image, transpose(radim3),/scalable,/colbar
      plt_image, transpose(radim4),/scalable,/colbar
      ;
      radim1=reform(weight(0,*,*,*))
      radim2=reform(weight(1,*,*,*))
      radim3=reform(weight(2,*,*,*))
      radim4=reform(weight(3,*,*,*))
      radim1=reform(radim1,nx*ny,nz,/overwrite)
      radim2=reform(radim2,nx*ny,nz,/overwrite)
      radim3=reform(radim3,nx*ny,nz,/overwrite)
      radim4=reform(radim4,nx*ny,nz,/overwrite)
      ;
      plt_image, transpose(radim1),/scalable,/colbar
      plt_image, transpose(radim2),/scalable,/colbar
      plt_image, transpose(radim3),/scalable,/colbar
      plt_image, transpose(radim4),/scalable,/colbar
      device,/close
      ;
   endfor
   ;
   help,radio
   radiotot=radiotot+radio
   weighttot=weighttot+weight
   noisetot=noisetot+noiseradio
   ;
   openw,1,fdir+field+'.radio_b4svd'
   printf,1,radio
   close,1
   ;
   openw,1,fdir+field+'.weight_b4svd'
   printf,1,weight
   close,1
   ;
   openw,1,fdir+field+'.noisecal_b4svd'
   printf,1,noiseradio
   close,1
   ;
   openw,1,fdir+field+'.ra'
   printf,1,ra
   close,1
   ;
   openw,1,fdir+field+'.dec'
   printf,1,dec
   close,1
   ;
   radio2=radio
   ind=where(weight gt 0.,cind)
   radio2(ind)=radio2(ind)/weight(ind)
   ;
   openw,1,fdir+field+'.radio_weighted_b4svd'
   printf,1,radio2
   close,1
   ;
   device,filename=fdir+'Tsky_b4svd.ps'
   radim1=reform(radio2(0,*,*,*))
   radim2=reform(radio2(1,*,*,*))
   radim3=reform(radio2(2,*,*,*))
   radim4=reform(radio2(3,*,*,*))
   radim1=reform(radim1,nx*ny,nz,/overwrite)
   radim2=reform(radim2,nx*ny,nz,/overwrite)
   radim3=reform(radim3,nx*ny,nz,/overwrite)
   radim4=reform(radim4,nx*ny,nz,/overwrite)
   ;
   plt_image, transpose(radim1),/scalable,/colbar
   plt_image, transpose(radim2),/scalable,/colbar
   plt_image, transpose(radim3),/scalable,/colbar
   plt_image, transpose(radim4),/scalable,/colbar
   device,/close
   ;
   radio2=noiseradio
   ind=where(weight gt 0.,cind)
   radio2(ind)=radio2(ind)/weight(ind)
   ;
   openw,1,fdir+field+'.noisecal_weighted_b4svd'
   printf,1,radio2
   close,1
   ;
   device,filename=fdir+'Tsky_noise_b4svd.ps'
   radim1=reform(radio2(0,*,*,*))
   radim2=reform(radio2(1,*,*,*))
   radim3=reform(radio2(2,*,*,*))
   radim4=reform(radio2(3,*,*,*))
   radim1=reform(radim1,nx*ny,nz,/overwrite)
   radim2=reform(radim2,nx*ny,nz,/overwrite)
   radim3=reform(radim3,nx*ny,nz,/overwrite)
   radim4=reform(radim4,nx*ny,nz,/overwrite)
   ;
   plt_image, transpose(radim1),/scalable,/colbar
   plt_image, transpose(radim2),/scalable,/colbar
   plt_image, transpose(radim3),/scalable,/colbar
   plt_image, transpose(radim4),/scalable,/colbar
   device,/close
   ;
   ;
endfor
;
;
; total file outputs per day
openw,1,fdir+field+'.radiotot_b4svd'
printf,1,radiotot
close,1
;
openw,1,fdir+field+'.weighttot_b4svd'
printf,1,weighttot
close,1
;
ind=where(weighttot gt 0,cind)
radiotot(ind)=radiotot(ind)/weighttot(ind)
;
openw,1,fdir+field+'.radiotot_weighted_b4svd'
printf,1,radiotot
close,1
;
openw,1,fdir+field+'.noisetot_b4svd'
printf,1,noisetot
close,1
;
ind=where(weighttot gt 0,cind)
noisetot(ind)=noisetot(ind)/weighttot(ind)
;
openw,1,fdir+field+'.noisetot_weighted_b4svd'
printf,1,noisetot
close,1
;
;
;
end
