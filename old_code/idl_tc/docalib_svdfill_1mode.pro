pro docalib_svdfill_1mode,radio,weight,t_off,radio_on,weight_on,field=field,fn=fn,sigmacut=sigmacut,dname=dname,svdfg=svdfg,fitsfile=fitsfile,sim=sim,signal=signal,daisy=daisy,drift=drift

; field:  data field
; fn:  fits file name

if not keyword_set(field) then field='zcosmos'  ; for the radio/optical cube setup
if not keyword_set(signal) then signal='100muK'
if not keyword_set(dname) then dname='02'
if keyword_set(daisy) then obsmode='daisy' else obsmode='drift'
if not keyword_set(daisy) then drift=drift

;ffdir=['9-17','18-26']  ;for 00
if dname eq '02' then begin
   if keyword_set(daisy) then ffdir=['30'] 
   if keyword_set(drift) then ffdir=['scan9-17','scan18-26','scan27-28']
;   if keyword_set(drift) then ffdir=['scan18-26','scan27-28']
endif

set_plot,'ps'
device,/color, bits_per_pixel=8

; GBT parameters
beam=15./60.d0 ; beam size
dra=beam/2.d
ddec=beam/2.d
;dra=3./60.d0 ;arcmin
;ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq_rest=1420.405751 ;MHz
nfreq=2048
Npol=4

; field center
;f0=[150.11667,2.2269444]
;fbin=[40,40]
f0=[150.11635,2.2247776]
fbin=[9,9]
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
radiotot20mode=radio
weighttot20mode=weight
radiotot25mode=radio
weighttot25mode=weight
radiotot10mode=radio
weighttot10mode=weight
radiotot5mode=radio
weighttot5mode=weight
radiotot0mode=radio
weighttot0mode=weight

; parameters
nfile=n_elements(ffdir)
nn=16320L  ;255.*8*8  8 freq windows and 4 pol x 2 on-off
nscan=9
;nscan=1
ntime=255
npol=4
nwin=8

;
for ifile=0,nfile-1 do begin
   ;
   fdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed/'+dname+'/'+ffdir(ifile)+'/'
   spawn,'mkdir '+fdir
;   fn='/cita/d/raid-cita/tchang/GBT/GBT09C_075/00_zCOSMOS_drift_'+ffdir(ifile)+'.raw.acs.fits'
   fn='/cita/d/raid-project/gmrt/tchang/GBT/GBT09C_075/'+dname+'_zCOSMOS_'+obsmode+'_'+ffdir(ifile)+'.raw.acs.fits'
   ;
   nscan=9
   if ffdir(ifile) eq 'scan27-28' then nscan=2
   if ffdir(ifile) eq '30' then nscan=1
   ;
   radio=dblarr(Npol,nx,ny,nz)
   weight=radio
   radio20mode=dblarr(Npol,nx,ny,nz)
   weight20mode=radio
   radio25mode=dblarr(Npol,nx,ny,nz)
   weight25mode=radio
   radio10mode=dblarr(Npol,nx,ny,nz)
   weight10mode=radio
   radio5mode=dblarr(Npol,nx,ny,nz)
   weight5mode=radio
   radio0mode=dblarr(Npol,nx,ny,nz)
   weight0mode=radio
   noiseradio=dblarr(Npol,nx,ny,nz)
   noiseradio20mode=dblarr(Npol,nx,ny,nz)
   noiseradio25mode=dblarr(Npol,nx,ny,nz)
   noiseradio10mode=dblarr(Npol,nx,ny,nz)
   noiseradio5mode=dblarr(Npol,nx,ny,nz)
   noiseradio0mode=dblarr(Npol,nx,ny,nz)
   filenum=strtrim(string(ifile),2)
   ;
   ;
   time_freq_tot=dblarr(ntime*nscan,npol,nfreq*nwin)
   freqkeep=dblarr(nfreq*nwin)
   time_freq_matrix=dblarr(ntime*nscan,npol,nfreq*3/5*nwin)
   time_freq_crosspol=time_freq_matrix
   time_freq_noise=dblarr(ntime*nscan,npol,nfreq*3/5*nwin)
   flagtot=time_freq_matrix
   time_freq_flag=time_freq_matrix
      ; collect the frequencies,ra,dec as well
   freqtot=dblarr(nfreq*3/5*nwin)
   ratot=dblarr(ntime*nscan)
   dectot=dblarr(ntime*nscan)
   ;
   flagnumon=0L
   flagnumoff=0L
   flagnumoff20mode=0L
   flagnumoff25mode=0L
   flagnumoff10mode=0L
   flagnumoff5mode=0L
   flagnumoff0mode=0L
   ;
   ;
   if keyword_set(fitsfile) then begin
      ;
      for scan=0L,nscan-1L do begin
         
         print,'******'
         print,'** processing scan ***',scan
         print,'******'
         ;
         scann=strtrim(string(scan),2)
         rec0=nn*scan
         rec1=nn*(scan+1)-1L
         ;
         fdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed/'+dname+'/'+ffdir(ifile)+'/'
         print,'fdir: ',fdir
         print,'fn: ',fn
         print,'rec0',rec0
         print,'rec1',rec1
         zmax_rec=nz
         zmin_rec=nz
         ;
         for i=0,nwin-1 do begin
            ;
            calibration_svd,t_off,t_on,noise,flag_off,flag_on,ra,dec,freq,pol,parallatic,recn=[rec0,rec1],thisif=i,threshould=0.02,fn=fn,drift=drift,daisy=daisy  
            ;
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
            time_freq_tot((scan*ntime):(scan+1)*ntime-1,*,i*nfreq:(i+1)*nfreq-1)=tbin
            freqkeep(i*nfreq:(i+1)*nfreq-1)=freq
            time_freq_matrix((scan*ntime):(scan+1)*ntime-1,*,i*nkeep:(i+1)*nkeep-1)=tbin(*,*,nfreq/5+1:nfreq*4/5-1)
            time_freq_crosspol((scan*ntime):(scan+1)*ntime-1,*,i*nkeep:(i+1)*nkeep-1)=tbinpol(*,*,nfreq/5+1:nfreq*4/5-1)
            time_freq_noise((scan*ntime):(scan+1)*ntime-1,*,i*nkeep:(i+1)*nkeep-1)=tbinnoise(*,*,nfreq/5+1:nfreq*4/5-1)
            flagtot((scan*ntime):(scan+1)*ntime-1,*,i*nkeep:(i+1)*nkeep-1)=tbinflag(*,*,nfreq/5+1:nfreq*4/5-1)
            freqtot(i*nkeep:(i+1)*nkeep-1)=freq(nfreq/5+1:nfreq*4/5-1)
            if (i eq 0) then begin
               ratot((scan*ntime):(scan+1)*ntime-1)=ra
               dectot((scan*ntime):(scan+1)*ntime-1)=dec
            endif
            ;
         endfor
         ;
      endfor
      ;
      ; now have collected all 6 scans into one big file
      ;
      device,filename=fdir+'Psky_freq_matched.mode.ps'
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
      sigth=5.                  ; sigma threshould
      for i=0L,ntime*nscan-1 do begin
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
               time_freq_crosspol(i,j,ind)=0.
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
               time_freq_crosspol(i,j,ind)=0.
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
      ; write outputs here, per 6 scan, per file, per day
      openw,1,fdir+'radio_tf_b4svd'
      printf,1,time_freq_flag
      close,1
      ;
      openw,1,fdir+'flagtot_tf_b4svd'
      printf,1,flagtot
      close,1
      ;
      openw,1,fdir+'noisecal_tf_b4svd'
      printf,1,time_freq_noise
      close,1
      ;
      openw,1,fdir+'crosspol_tf_b4svd'
      printf,1,time_freq_crosspol
      close,1
      ;
      openw,1,fdir+'freqtot'
      printf,1,freqtot
      close,1
      ;
      openw,1,fdir+'ra'
      printf,1,ratot
      close,1
      ;
      openw,1,fdir+'dec'
      printf,1,dectot
      close,1
      ;
      ;
      ;
      device,filename=fdir+'Pnoise_freq_matched.mode.ps'
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
   endif else begin
      ;
      ; read in the data files
      openr,1,fdir+'radio_tf_b4svd'
      readf,1,time_freq_flag
      close,1
      ; 
      openr,1,fdir+'flagtot_tf_b4svd'
      readf,1,flagtot
      close,1
      ;
      openr,1,fdir+'noisecal_tf_b4svd'
      readf,1,time_freq_noise
      close,1
      ;
      openr,1,fdir+'freqtot'
      readf,1,freqtot
      close,1
      ;
      openr,1,fdir+'ra'
      readf,1,ratot
      close,1
      ;
      openr,1,fdir+'dec'
      readf,1,dectot
      close,1
      ;   
   endelse
      ;
   if keyword_set(sim) then begin
      ;
;      for scan=0,nscan-1 do begin
      radiotfsim=dblarr(ntime*nscan,npol,nfreq*3/5*nwin)
      openr,1,fdir+'simradio_tf_b4svd_'+signal
      readf,1,radiotfsim
      close,1
      help,time_freq_flag
print,min(time_freq_flag),max(time_freq_flag),mean(time_freq_flag)
      time_freq_flag=time_freq_flag+radiotfsim
      print,min(time_freq_flag),max(time_freq_flag),mean(time_freq_flag)
;read,test
         ;
;;         for ipol=0,npol-1 do begin
;         time_freq_flag((scan*ntime):(scan+1)*ntime-1,*,*)=time_freq_flag((scan*ntime):(scan+1)*ntime-1,*,*)+radiotfsim
;;         endfor
         ;
;      endfor
      ;
      fdir='/cita/d/raid-project/gmrt/tchang/GBT/zCOSMOS/Processed_sim_'+signal+'/'+dname+'/'+ffdir(ifile)+'/'
      spawn,'mkdir '+fdir
      ;
   endif
   ;
   if keyword_set(svdfg) then begin
      ;
      ind=where(flagtot eq 1, cind)
      time_freq_flag(ind)=time_freq_flag(ind)-1.
      ;
                                ; do frequency rebin to reduce SVD
                                ; workload -- bin into 2Mpc redshift
                                ;             bins
      time_z_flag=dblarr(ntime*nscan,npol,nz)
      flagtotz=dblarr(ntime*nscan,npol,nz)
      time_z_noise=time_z_flag
 ;     freqnum=dblarr(ntime*nscan,npol,nz)
      Da_z,(freq_rest/freqtot-1.),distslice
      distbin=(floor((distslice-d0)/dd)) ;d0 is at bin0
      for it=0,ntime*nscan-1 do begin
         for ipol=0,npol-1 do begin
            flagslice=reform(flagtot(it,ipol,*))
            for iz=0,nz-1 do begin
               ind=where(distbin eq iz and flagslice eq 1,cind)
               if cind gt 0 then begin
                  time_z_flag(it,ipol,iz)=total(time_freq_flag(it,ipol,ind))/double(cind)
                  time_z_noise(it,ipol,iz)=total(time_freq_noise(it,ipol,ind))/double(cind)
                  flagtotz(it,ipol,iz)=cind
               endif
  ;             freqnum(it,ipol,iz)=cind
            endfor
         endfor
      endfor
      ;
      time_freq_flag=time_z_flag
      time_freq_noise=time_z_noise
      flagtot=flagtotz
      ;
      ; produce copies of flagtot and time_freq_flag and time_freq_noise
      flagtot20mode=flagtot
      flagtot25mode=flagtot
      flagtot10mode=flagtot
      flagtot5mode=flagtot
      flagtot0mode=flagtot
      ;
      time_freq_flag20mode=time_freq_flag
      time_freq_flag25mode=time_freq_flag
      time_freq_flag10mode=time_freq_flag
      time_freq_flag5mode=time_freq_flag
      ;
      ;
      ; svd
      ; do each pol separately
      ;
      radiosvd5mode=dblarr(ntime*nscan,npol,nz)
      radiosvd10mode=radiosvd5mode
      radiosvd15mode=radiosvd5mode
      radiosvd20mode=radiosvd5mode
      radiosvd25mode=radiosvd5mode
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
         w5=dblarr(n_elements(s))
         w10=w5
         w15=w5
         w20=w5
         w25=w5
         w5(0:0)=s(0:0)
         w10(0:1)=s(0:1)
         w15(0:2)=s(0:2)
         w20(0:3)=s(0:3)
         w25(0:5)=s(0:5)
         ;
         radiosvd5mode(*,ii,*) = u ## diag_matrix(w5) ## transpose(v)
         radiosvd10mode(*,ii,*) =u ## diag_matrix(w10) ## transpose(v)
         radiosvd15mode(*,ii,*) = u ## diag_matrix(w15) ## transpose(v)
         radiosvd20mode(*,ii,*) =u ## diag_matrix(w20) ## transpose(v)
         radiosvd25mode(*,ii,*) = u ## diag_matrix(w25) ## transpose(v)
         ;
      endfor
      ;
      radiores5mode = time_freq_flag-radiosvd5mode
      radiores10mode = time_freq_flag-radiosvd10mode
      radiores15mode = time_freq_flag-radiosvd15mode
      radiores20mode = time_freq_flag-radiosvd20mode
      radiores25mode = time_freq_flag-radiosvd25mode
      radiores0mode = time_freq_flag
      ;
      ; calculate thermal
      dnu=50e6/double(nfreq)
      dt=0.5d
 ;     freqnum=freqnum>1
      thermal=(dnu*dt)*(flagtot)
      ;
      device,filename=fdir+'Psky_svd1.mode.ps'
      ;calculate the weight per frequency and apply
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores0mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores0mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores5mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores5mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores10mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores10mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
       endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores15mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores15mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores20mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores20mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores25mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores25mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      device,/close
      ;
      ; do sigma flagging
      dnu=50e6/double(nfreq)
      dt=0.5d
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot>1)
      sigma_th=5.
      ;
      ;
      ind=where(abs(radiores0mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores0mode(ind)=0
         flagtot0mode(ind)=0
      endif
      ;
      ind=where(abs(radiores15mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores15mode(ind)=0
         flagtot(ind)=0
         time_freq_noise(ind)=0.
      endif
      ;
      ind=where(abs(radiores20mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores20mode(ind)=0
         flagtot20mode(ind)=0
      endif
      ;
      ind=where(abs(radiores25mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores25mode(ind)=0
         flagtot25mode(ind)=0
      endif
      ;
      ind=where(abs(radiores5mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores5mode(ind)=0
         flagtot5mode(ind)=0
      endif
      ;
      ind=where(abs(radiores10mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores10mode(ind)=0
         flagtot10mode(ind)=0
      endif
      ;
      ;
      device,filename=fdir+'Psky_svd1.scut.mode.ps'
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot0mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores0mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores0mode(*,0,ifreq)*tfweight*(flagtot0mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot5mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores5mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores5mode(*,0,ifreq)*tfweight*(flagtot5mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot10mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores10mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores10mode(*,0,ifreq)*tfweight*(flagtot10mode(*,0,ifreq)<1)
       endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores15mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores15mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot20mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores20mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores20mode(*,0,ifreq)*tfweight*(flagtot20mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot25mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores25mode(*,0,ifreq)))<(thermal(*,0,ifreq))
         tfradio(*,ifreq)=radiores25mode(*,0,ifreq)*tfweight*(flagtot25mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      ;
      device,/close
      ;
      ; fill in the svd values to the holes
      ind=where(flagtot eq 0, cind)
      time_freq_flag(ind)=radiosvd15mode(ind)
      ;   
      ind=where(flagtot20mode eq 0, cind)
      time_freq_flag20mode(ind)=radiosvd20mode(ind)
      ;
      ind=where(flagtot25mode eq 0, cind)
      time_freq_flag25mode(ind)=radiosvd25mode(ind)
      ;
      ind=where(flagtot10mode eq 0, cind)
      time_freq_flag10mode(ind)=radiosvd10mode(ind)
      ;
      ind=where(flagtot5mode eq 0, cind)
      time_freq_flag5mode(ind)=radiosvd5mode(ind)
      ;
      ;
      ; do frequency weighting
      freqweight=dblarr(npol,nz)
      freqweight20mode=dblarr(npol,nz)
      freqweight25mode=dblarr(npol,nz)
      freqweight10mode=dblarr(npol,nz)
      freqweight5mode=dblarr(npol,nz)
      ;
      for ii=0,npol-1 do begin
         for kk=0,nz-1 do begin
            ind=where(flagtot(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               freqweight(ii,kk)=(1./variance(radiores15mode(ind,ii,kk)))<mean((dnu*dt)*(flagtot(ind,ii,kk)>1))
               time_freq_flag(ind,ii,kk)=time_freq_flag(ind,ii,kk)*freqweight(ii,kk)
            endif
            ;
            ind=where(flagtot20mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               freqweight20mode(ii,kk)=(1./variance(radiores20mode(ind,ii,kk)))<mean((dnu*dt)*(flagtot20mode(ind,ii,kk)>1))
               time_freq_flag20mode(ind,ii,kk)=time_freq_flag20mode(ind,ii,kk)*freqweight20mode(ii,kk)
            endif
            ;
            ind=where(flagtot25mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               freqweight25mode(ii,kk)=(1./variance(radiores25mode(ind,ii,kk)))<mean((dnu*dt)*(flagtot25mode(ind,ii,kk)>1))
               time_freq_flag25mode(ind,ii,kk)=time_freq_flag25mode(ind,ii,kk)*freqweight25mode(ii,kk)
            endif
            ;
            ind=where(flagtot10mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               freqweight10mode(ii,kk)=(1./variance(radiores10mode(ind,ii,kk)))<mean((dnu*dt)*(flagtot10mode(ind,ii,kk)>1))
               time_freq_flag10mode(ind,ii,kk)=time_freq_flag10mode(ind,ii,kk)*freqweight10mode(ii,kk)
            endif
            ;
            ind=where(flagtot5mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               freqweight5mode(ii,kk)=(1./variance(radiores5mode(ind,ii,kk)))<mean((dnu*dt)*(flagtot5mode(ind,ii,kk)>1))
               time_freq_flag5mode(ind,ii,kk)=time_freq_flag5mode(ind,ii,kk)*freqweight5mode(ii,kk)
            endif
            ;
         endfor
      endfor
      ;
      ;
      openw,2,fdir+field+'.timevec_1mode'
      openw,3,fdir+field+'.timevec_2mode'
      openw,4,fdir+field+'.timevec_3mode'
      openw,5,fdir+field+'.timevec_4mode'
      openw,6,fdir+field+'.timevec_6mode'
      openw,12,fdir+field+'.freqvec_1mode'
      openw,13,fdir+field+'.freqvec_2mode'
      openw,14,fdir+field+'.freqvec_3mode'
      openw,15,fdir+field+'.freqvec_4mode'
      openw,16,fdir+field+'.freqvec_6mode'
      openw,22,fdir+field+'.eigenvalue_1mode'
      openw,23,fdir+field+'.eigenvalue_2mode'
      openw,24,fdir+field+'.eigenvalue_3mode'
      openw,25,fdir+field+'.eigenvalue_4mode'
      openw,26,fdir+field+'.eigenvalue_6mode'
      for ii=0,npol-1 do begin
         ;
         radio2d=reform(time_freq_flag(*,ii,*))
         la_svd,radio2d, s,u,v,/double,status=status
         print,'SVD status:',status
         help,s
         help,u
         help,v
         w5=dblarr(n_elements(s))
         w10=w5
         w15=w5
         w20=w5
         w25=w5
         w15(0:2)=s(0:2)
         ;
         printf,4,u
         printf,14,v
         printf,24,s
         ;
         radiosvd15mode(*,ii,*) = u ## diag_matrix(w15) ## transpose(v)
         ;
         radio2d=reform(time_freq_flag20mode(*,ii,*))
         la_svd,radio2d, s,u,v,/double,status=status
         w20(0:3)=s(0:3)
         radiosvd20mode(*,ii,*) = u ## diag_matrix(w20) ## transpose(v)
         ;
         printf,5,u
         printf,15,v
         printf,25,s
         ;
         radio2d=reform(time_freq_flag25mode(*,ii,*))
         la_svd,radio2d, s,u,v,/double,status=status
         w25(0:5)=s(0:5)
         radiosvd25mode(*,ii,*) = u ## diag_matrix(w25) ## transpose(v)
         ;
         printf,6,u
         printf,16,v
         printf,26,s
         ;
         radio2d=reform(time_freq_flag10mode(*,ii,*))
         la_svd,radio2d, s,u,v,/double,status=status
         w10(0:1)=s(0:1)
         radiosvd10mode(*,ii,*) = u ## diag_matrix(w10) ## transpose(v)
         ;
         printf,3,u
         printf,13,v
         printf,23,s
         ;
         radio2d=reform(time_freq_flag5mode(*,ii,*))
         la_svd,radio2d, s,u,v,/double,status=status
         w5(0:0)=s(0:0)
         radiosvd5mode(*,ii,*) = u ## diag_matrix(w5) ## transpose(v)
         ;
         printf,2,u
         printf,12,v
         printf,22,s
         ;
      endfor
      ;
      close,2
      close,3
      close,4
      close,5
      close,6
      close,12
      close,13
      close,14
      close,15
      close,16
      close,22
      close,23
      close,24
      close,25
      close,26
      ;
      radiores5mode = time_freq_flag5mode-radiosvd5mode
      radiores10mode = time_freq_flag10mode-radiosvd10mode
      radiores15mode = time_freq_flag-radiosvd15mode
      radiores20mode = time_freq_flag20mode-radiosvd20mode
      radiores25mode = time_freq_flag25mode-radiosvd25mode
      ; de-freq weight
      for ii=0,npol-1 do begin
         for kk=0,nz-1 do begin
            ind=where(flagtot(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               radiores15mode(ind,ii,kk)=radiores15mode(ind,ii,kk)/freqweight(ii,kk)
            endif
            ind=where(flagtot20mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               radiores20mode(ind,ii,kk)=radiores20mode(ind,ii,kk)/freqweight20mode(ii,kk)
            endif
            ind=where(flagtot25mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               radiores25mode(ind,ii,kk)=radiores25mode(ind,ii,kk)/freqweight25mode(ii,kk)
            endif
            ind=where(flagtot10mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               radiores10mode(ind,ii,kk)=radiores10mode(ind,ii,kk)/freqweight10mode(ii,kk)
            endif
            ind=where(flagtot5mode(*,ii,kk) gt 0,cind)
            if cind gt 1 then begin
               radiores5mode(ind,ii,kk)=radiores5mode(ind,ii,kk)/freqweight5mode(ii,kk)
            endif
         endfor
      endfor
      ; do sigma flagging
      dnu=50e6/double(nfreq)
      dt=0.5d
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot>1)
      sigma_th=5.
      ;
      ind=where(abs(radiores15mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores15mode(ind)=0
         flagtot(ind)=0
         time_freq_noise(ind)=0.
      endif
      ;
      ;
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot20mode>1)
      ind=where(abs(radiores20mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores20mode(ind)=0
         flagtot20mode(ind)=0
      endif
      ;
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot25mode>1)
      ind=where(abs(radiores25mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores25mode(ind)=0
         flagtot25mode(ind)=0
      endif
      ;
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot10mode>1)
      ind=where(abs(radiores10mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores10mode(ind)=0
         flagtot10mode(ind)=0
      endif
      ;
      sigma=1./sqrt(dnu*dt)/sqrt(flagtot5mode>1)
      ind=where(abs(radiores5mode) gt (sigma_th*sigma),cind)
      if (cind gt 0) then begin
         radiores5mode(ind)=0
         flagtot5mode(ind)=0
      endif
      ;
      ;
      ;
      device,filename=fdir+'Psky_svd2.mode.ps'
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot0mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores0mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores0mode(*,0,ifreq)*tfweight*(flagtot0mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_0mode.dat'
      printf,1,tfradio
      close,1

      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot5mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores5mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores5mode(*,0,ifreq)*tfweight*(flagtot5mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_1mode.dat'
      printf,1,tfradio
      close,1

      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot10mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores10mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores10mode(*,0,ifreq)*tfweight*(flagtot10mode(*,0,ifreq)<1)
       endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_2mode.dat'
      printf,1,tfradio
      close,1
      ;
      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores15mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores15mode(*,0,ifreq)*tfweight*(flagtot(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_3mode.dat'
      printf,1,tfradio
      close,1

      ;   
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot20mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores20mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores20mode(*,0,ifreq)*tfweight*(flagtot20mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_4mode.dat'
      printf,1,tfradio
      close,1

      ;
      tfradio=dblarr(ntime*nscan,nz)
      thermal=(dnu*dt)*(flagtot25mode)
      for ifreq=0,nz-1 do begin
         tfweight=(1./variance(radiores25mode(*,0,ifreq)))<thermal(*,0,ifreq)
         tfradio(*,ifreq)=radiores25mode(*,0,ifreq)*tfweight*(flagtot25mode(*,0,ifreq)<1)
      endfor
      plt_image, transpose(tfradio),/scalable,/colbar
      openw,1,fdir+field+'.Psky_svd2_6mode.dat'
      printf,1,tfradio
      close,1
      ;
      device,/close
      ;
      ;
      device,filename=fdir+field+'.Pnoise_flag.mode.ps'
      plt_image,transpose(reform(time_freq_noise(*,0,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,1,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,2,*))),/scalable,/colbar
      plt_image,transpose(reform(time_freq_noise(*,3,*))),/scalable,/colbar
      device,/close
      ;
   endif
      ;
      ; rest frequency of the 8 IF's 
      ;
      ; 695, 725, 755, 785, 815, 845, 875, 905  in MHz
      ;
      ; each 50 MHz bandwidth
      ;
   if not keyword_set(svdfg) then begin
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
      mk_hicube_noise,radiores15mode,flagtot,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio,weight=weight,$
                      noiseradio=noiseradio,field=field,flagnum=flagnumoff
      mk_hicube_noise,radiores20mode,flagtot20mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio20mode,weight=weight20mode,$
                      field=field,noiseradio=noiseradio20mode,flagnum=flagnumoff20mode
      mk_hicube_noise,radiores25mode,flagtot25mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio25mode,weight=weight25mode,$
                      field=field,noiseradio=noiseradio25mode,flagnum=flagnumoff25mode
      ;
      ind=where(radio ne 0,cind)
      print,''
      print,'outside mk_hicube:  non-zero radio elements:',cind
      print,''
      ;
      ;
   endif else begin
      ;
      mk_hicube_half,radiores15mode,flagtot,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio,weight=weight,$
                      noiseradio=noiseradio,field=field,flagnum=flagnumoff
      mk_hicube_half,radiores20mode,flagtot20mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio20mode,weight=weight20mode,$
                      field=field,noiseradio=noiseradio20mode,flagnum=flagnumoff20mode
      mk_hicube_half,radiores25mode,flagtot25mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio25mode,weight=weight25mode,$
                      field=field,noiseradio=noiseradio25mode,flagnum=flagnumoff25mode
      
      mk_hicube_half,radiores10mode,flagtot10mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio10mode,weight=weight10mode,$
                      field=field,noiseradio=noiseradio10mode,flagnum=flagnumoff10mode

      mk_hicube_half,radiores5mode,flagtot5mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio5mode,weight=weight5mode,$
                      field=field,noiseradio=noiseradio5mode,flagnum=flagnumoff5mode

      mk_hicube_half,radiores0mode,flagtot0mode,time_freq_noise,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio0mode,weight=weight0mode,$
                      field=field,noiseradio=noiseradio0mode,flagnum=flagnumoff0mode
      ;
   endelse
;
   device,filename=fdir+'Psky_xz_aftsvd.mode.ps'
   ;
   ;
   radio2=radio0mode
   ind=where(weight0mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight0mode(ind)
   ;
   thermal=(dnu*dt)*(weight0mode)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight0mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight0mode(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   radio2=radio5mode
   ind=where(weight5mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight5mode(ind)
   ;
   thermal=(dnu*dt)*(weight5mode)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight5mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight5mode(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   radio2=radio10mode
   ind=where(weight10mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight10mode(ind)
   ;
   thermal=(dnu*dt)*(weight10mode>1)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight10mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight10mode(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   ;
   radio2=radio
   ind=where(weight gt 0.,cind)
   radio2(ind)=radio2(ind)/weight(ind)
   ;
   thermal=(dnu*dt)*(weight>1)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   radio2=radio20mode
   ind=where(weight20mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight20mode(ind)
   ;
   thermal=(dnu*dt)*(weight20mode>1)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight20mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight20mode(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   radio2=radio25mode
   ind=where(weight25mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight25mode(ind)
   ;
   thermal=(dnu*dt)*(weight25mode>1)
   radiowt=dblarr(nx,ny,nz)
   for iz=0,nz-1 do begin
      ind=where(weight25mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<thermal(0,*,*,iz)
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weight25mode(0,*,*,iz)<1)
      endif
   endfor
   radiowt=reform(radiowt,nx*ny,nz)
   plt_image, transpose(radiowt),/scalable,/colbar
   ;
   device,/close
   ;
   radiotot=radiotot+radio
   weighttot=weighttot+weight
   noisetot=noisetot+noiseradio
   radiotot20mode=radiotot20mode+radio20mode
   weighttot20mode=weighttot20mode+weight20mode
   radiotot25mode=radiotot25mode+radio25mode
   weighttot25mode=weighttot25mode+weight25mode
   radiotot10mode=radiotot10mode+radio10mode
   weighttot10mode=weighttot10mode+weight10mode
   radiotot5mode=radiotot5mode+radio5mode
   weighttot5mode=weighttot5mode+weight5mode
   radiotot0mode=radiotot0mode+radio0mode
   weighttot0mode=weighttot0mode+weight0mode
   ;
   ;
   openw,1,fdir+field+'.radio3mode_aftsvd'
   printf,1,radio
   close,1
   ;
   openw,1,fdir+field+'.weight3mode_aftsvd'
   printf,1,weight
   close,1
   ;
   openw,1,fdir+field+'.radio4mode_aftsvd'
   printf,1,radio20mode
   close,1
   ;
   openw,1,fdir+field+'.weight4mode_aftsvd'
   printf,1,weight20mode
   close,1
   ;
   openw,1,fdir+field+'.radio6mode_aftsvd'
   printf,1,radio25mode
   close,1
   ;
   openw,1,fdir+field+'.weight6mode_aftsvd'
   printf,1,weight25mode
   close,1
   ;
   openw,1,fdir+field+'.radio2mode_aftsvd'
   printf,1,radio10mode
   close,1
   ;
   openw,1,fdir+field+'.weight2mode_aftsvd'
   printf,1,weight10mode
   close,1
   ;
   openw,1,fdir+field+'.radio1mode_aftsvd'
   printf,1,radio5mode
   close,1
   ;
   openw,1,fdir+field+'.weight1mode_aftsvd'
   printf,1,weight5mode
   close,1
   ;
   openw,1,fdir+field+'.radio0mode_aftsvd'
   printf,1,radio0mode
   close,1
   ;
   openw,1,fdir+field+'.weight0mode_aftsvd'
   printf,1,weight0mode
   close,1
   ;
   openw,1,fdir+field+'.noisecal_aftsvd'
   printf,1,noiseradio
   close,1
   ;
   openw,1,fdir+field+'.ra'
   printf,1,ratot
   close,1
   ;
   openw,1,fdir+field+'.dec'
   printf,1,dectot
   close,1
   ;
   radio2=radio
   ind=where(weight gt 0.,cind)
   radio2(ind)=radio2(ind)/weight(ind)
   ;
   openw,1,fdir+field+'.radio3mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=radio20mode
   ind=where(weight20mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight20mode(ind)
   ;
   openw,1,fdir+field+'.radio4mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=radio25mode
   ind=where(weight25mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight25mode(ind)
   ;
   openw,1,fdir+field+'.radio6mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=radio10mode
   ind=where(weight10mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight10mode(ind)
   ;
   openw,1,fdir+field+'.radio2mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=radio5mode
   ind=where(weight5mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight5mode(ind)
   ;
   openw,1,fdir+field+'.radio1mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=radio0mode
   ind=where(weight0mode gt 0.,cind)
   radio2(ind)=radio2(ind)/weight0mode(ind)
   ;
   openw,1,fdir+field+'.radio0mode_weighted_aftsvd'
   printf,1,radio2
   close,1
   ;
   radio2=noiseradio
   ind=where(weight gt 0.,cind)
   radio2(ind)=radio2(ind)/weight(ind)
   ;
;   openw,1,fdir+field+'.noisecal_weighted_aftsvd'
;   printf,1,radio2
;   close,1
   ;
;   device,filename=fdir+'Psky_noise_aftsvd.ps'
;   radim1=reform(radio2(0,*,*,*))
;   radim2=reform(radio2(1,*,*,*))
;   radim3=reform(radio2(2,*,*,*))
;   radim4=reform(radio2(3,*,*,*))
;   radim1=reform(radim1,nx*ny,nz,/overwrite)
;   radim2=reform(radim2,nx*ny,nz,/overwrite)
;   radim3=reform(radim3,nx*ny,nz,/overwrite)
;   radim4=reform(radim4,nx*ny,nz,/overwrite)
   ;
;   plt_image, transpose(radim1),/scalable,/colbar
;   plt_image, transpose(radim2),/scalable,/colbar
;   plt_image, transpose(radim3),/scalable,/colbar
;   plt_image, transpose(radim4),/scalable,/colbar
;   device,/close
   ;
   ;
endfor
;
;
; total file outputs per day
openw,1,fdir+field+'.radiotot3mode_aftsvd'
printf,1,radiotot
close,1
;
openw,1,fdir+field+'.weighttot3mode_aftsvd'
printf,1,weighttot
close,1
;
ind=where(weighttot gt 0,cind)
radiotot(ind)=radiotot(ind)/weighttot(ind)
;
openw,1,fdir+field+'.radiotot3mode_weighted_aftsvd'
printf,1,radiotot
close,1
;
;

openw,1,fdir+field+'.radiotot4mode_aftsvd'
printf,1,radiotot20mode
close,1
;
openw,1,fdir+field+'.weighttot4mode_aftsvd'
printf,1,weighttot20mode
close,1
;
ind=where(weighttot20mode gt 0,cind)
radiotot20mode(ind)=radiotot20mode(ind)/weighttot20mode(ind)
;
openw,1,fdir+field+'.radiotot4mode_weighted_aftsvd'
printf,1,radiotot20mode
close,1
;
openw,1,fdir+field+'.radiotot6mode_aftsvd'
printf,1,radiotot25mode
close,1
;
openw,1,fdir+field+'.weighttot6mode_aftsvd'
printf,1,weighttot25mode
close,1
;
ind=where(weighttot25mode gt 0,cind)
radiotot25mode(ind)=radiotot25mode(ind)/weighttot25mode(ind)
;
openw,1,fdir+field+'.radiotot6mode_weighted_aftsvd'
printf,1,radiotot25mode
close,1
;
;
openw,1,fdir+field+'.radiotot2mode_aftsvd'
printf,1,radiotot10mode
close,1
;
openw,1,fdir+field+'.weighttot2mode_aftsvd'
printf,1,weighttot10mode
close,1
;
ind=where(weighttot10mode gt 0,cind)
radiotot10mode(ind)=radiotot10mode(ind)/weighttot10mode(ind)
;
openw,1,fdir+field+'.radiotot2mode_weighted_aftsvd'
printf,1,radiotot10mode
close,1
;
;
openw,1,fdir+field+'.radiotot1mode_aftsvd'
printf,1,radiotot5mode
close,1
;
openw,1,fdir+field+'.weighttot1mode_aftsvd'
printf,1,weighttot5mode
close,1
;
ind=where(weighttot5mode gt 0,cind)
radiotot5mode(ind)=radiotot5mode(ind)/weighttot5mode(ind)
;
openw,1,fdir+field+'.radiotot1mode_weighted_aftsvd'
printf,1,radiotot5mode
close,1
;
;
openw,1,fdir+field+'.radiotot0mode_aftsvd'
printf,1,radiotot0mode
close,1
;
openw,1,fdir+field+'.weighttot0mode_aftsvd'
printf,1,weighttot0mode
close,1
;
ind=where(weighttot0mode gt 0,cind)
radiotot0mode(ind)=radiotot0mode(ind)/weighttot0mode(ind)
;
openw,1,fdir+field+'.radiotot0mode_weighted_aftsvd'
printf,1,radiotot0mode
close,1
;
;
device,filename=fdir+'Pskytot_aftsvd.mode.ps'
;
radio2=radiotot0mode
;
thermal=(dnu*dt)*(weighttot0mode)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   ind=where(weighttot0mode(0,*,*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2(0,*,*,iz)
      zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
      radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot0mode(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
;
radio2=radiotot5mode
;
thermal=(dnu*dt)*(weighttot5mode>1)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   ind=where(weighttot5mode(0,*,*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2(0,*,*,iz)
      zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
      radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot5mode(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
radio2=radiotot10mode
;
thermal=(dnu*dt)*(weighttot10mode>1)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   ind=where(weighttot10mode(0,*,*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2(0,*,*,iz)
      zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
      radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot10mode(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
radio2=radiotot
;
thermal=(dnu*dt)*(weighttot>1)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   ind=where(weighttot(0,*,*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2(0,*,*,iz)
      zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
      radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
radio2=radiotot20mode
;
thermal=(dnu*dt)*(weighttot20mode>1)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
      ind=where(weighttot20mode(0,*,*,iz) gt 0,cind)
      if cind gt 1 then begin
         radiotmp=radio2(0,*,*,iz)
         zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
         radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot20mode(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
radio2=radiotot25mode
;
thermal=(dnu*dt)*(weighttot25mode>1)
radiowt=dblarr(nx,ny,nz)
for iz=0,nz-1 do begin
   ind=where(weighttot25mode(0,*,*,iz) gt 0,cind)
   if cind gt 1 then begin
      radiotmp=radio2(0,*,*,iz)
      zweight=(1./variance(radiotmp(ind)))<(thermal(0,*,*,iz))
      radiowt(*,*,iz)=radio2(0,*,*,iz)*zweight*(weighttot25mode(0,*,*,iz)<1)
   endif
endfor
radiowt=reform(radiowt,nx*ny,nz)
plt_image, transpose(radiowt),/scalable,/colbar
;
device,/close
;
;openw,1,fdir+field+'.noisetot_aftsvd'
;printf,1,noisetot
;close,1
;
;ind=where(weighttot gt 0,cind)
;noisetot(ind)=noisetot(ind)/weighttot(ind)
;
;openw,1,fdir+field+'.noisetot_weighted_aftsvd'
;printf,1,noisetot
;close,1
;
;
;
end

