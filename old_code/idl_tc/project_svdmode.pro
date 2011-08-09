pro project_svdmode,iradio,field=field,sim=sim,rebin=rebin

; svd dimension time vector = [z,t]=[560, 2154] = [560, 2240]

if not keyword_set(field) then field='f3'
if not keyword_set(signal) then signal='100muK'
if keyword_set(sigcut) then sigcut='_5sigcut' else sigcut=''

;mode=[1]
;smode=['1mode']
mode=[1,2,3,4,5,6,7,8,9,10,11,12,15,20,25]
smode=['1mode','2mode','3mode','4mode','5mode','6mode','7mode','8mode','9mode','10mode','11mode','12mode','15mode','20mode','25mode']


nmode=n_elements(smode)
;
if field eq 'f3' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then
;   thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed'+sigcut+'/'
   print,'thisdir: ',thisdir
   fdir=['/f3_53-58/','/f3_45-50/','/f3_27-29/']
   subdir=['/aug29/f3_05-10/', '/aug29/f3_15-20/',  '/aug29/f3_21-26/',  '/aug29/f3_27-32/',  $
           '/aug29/f3_37-42/',  '/aug29/f3_47-52/',  '/aug29/f3_53-58/', $
           '/aug30/f3_15-20/',  '/aug30/f3_21-26/',  '/aug30/f3_27-32/', $ 
           '/aug30/f3_33-38/',  '/aug30/f3_39-44/',  '/aug30/f3_45-50/', $ 
           '/aug30/f3_9-14/', '/sept07/f3_09-14/', '/sept07/f3_15-20/', $
           '/sept07/f3_21-26/', '/sept07/f3_27-29/']
   if not keyword_set(dname) then dname=['aug29','aug30','sept07']
   if n_elements(dname) eq 1 then begin
      if dname eq 'aug29' then begin
         fdir='/f3_53-58/'
         subdir=subdir(0:6)
      endif
      if dname eq 'aug30' then begin
         fdir='/f3_45-50/'
         subdir=subdir(7:13)
      endif
      if dname eq 'sept07' then begin
         fdir='/f3_27-29/'
         subdir=subdir(14:17)
      endif
   endif
   ffdir=thisdir+dname+fdir
   submapdir=thisdir+subdir
;   subfmdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'+subdir
   subfmdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'+subdir
   nx=40
   ny=6
   ;
endif


if field eq 'f4' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then
;   thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/';   if keyword_set(sim) then
;   thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'

   print,'thisdir: ',thisdir
   fdir=['/f4_69-74/','/f4_67-72/','/f4_33-36/']
   subdir=['/aug29/f4_63-68/', '/aug29/f4_69-74/','/aug30/f4_55-60/', $
           '/aug30/f4_61-66/', '/aug30/f4_67-72/',$
           '/sept06/f4_09-14/',  '/sept06/f4_15-20/',  '/sept06/f4_21-26/', $
           '/sept06/f4_27-32/',  '/sept06/f4_33-36/']
   if not keyword_set(dname) then dname=['aug29','aug30','sept06']
   if n_elements(dname) eq 1 then begin
      if dname eq 'aug29' then begin
         fdir='/f4_69-74/'
         subdir=subdir(0:1)
      endif
      if dname eq 'aug30' then begin
         fdir='/f4_67-72/'
         subdir=subdir(2:4)
      endif
      if dname eq 'sept06' then begin
         fdir='/f4_33-36/'
         subdir=subdir(5:9)
      endif
   endif
   ffdir=thisdir+dname+fdir
   submapdir=thisdir+subdir
;   subfmdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'+subdir
   subfmdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'+subdir
   nx=30
   ny=6
   ;
endif
set_plot,'ps'
device,/color, bits_per_pixel=8

; GBT parameters
beam=15./60.d0 ; beam size
dra=3./60.d0 ;arcmin
ddec=3./60.d0
d0=1450.  ;2Mpc per redshift bin starting at D=1450 h^-1 Mpc
dd=2.
nz=1120/dd
freq_rest=1420.405751 ;MHz

; parameters
nfile=n_elements(submapdir)
nn=22976L  ;359*8*8
nscan=6
ntime=359
npol=4
nwin=8
nfreq=2048
if field eq 'f4' then begin
   nn=19136
   ntime=299
endif
;
nz=560

!p.multi=0
;!p.multi=[0,1,5]
;device,filename=field+'.residual.ps'

radiotot=dblarr(nmode,npol,nx,ny,nz)
weighttot=radiotot
;
;for ifile=0,nfile-1 do begin
for ifile=0,0 do begin
   ;
   device,filename=field+'.residual.'+strtrim(string(ifile),2)+'.ps'
   ;fdir='/cita/d/raid-cita/tchang/GBT/Processed/'+dname+'/'+field+'_'+ffdir(ifile)+'/'
   fdir=submapdir(ifile)
;   fmdir=subfmdir(ifile)
   fmdir=fdir
   if fdir eq '/cita/d/raid-cita/tchang/GBT/Processed/sept07/f3_27-29/' then nscan=3
   if fdir eq '/cita/d/raid-cita/tchang/GBT/Processed/sept06/f4_33-36/' then nscan=4
   nt=ntime*nscan
   ;
   ;
   ; read in the simulation and rebin
   freqtot=dblarr(nfreq*3/5*nwin)
   ratot=dblarr(ntime*nscan)
   dectot=dblarr(ntime*nscan)
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
   ; read in the data files for test only
   if keyword_set(rebin) then begin

      radiotf=dblarr(nt,npol,nfreq*3/5*nwin)
      flagtot=radiotf
;   openr,1,fdir+'radio_tf_b4svd'
;   readf,1,radiotf
;   close,1
   ; 
   ; 
   openr,1,fdir+'simradio_tf_b4svd_'+signal
   readf,1,radiotf
   close,1
   ;
   flagtot=make_array(ntime*nscan,npol,nfreq*3/5*nwin,value=1.)
   openr,1,fdir+'flagtot_tf_b4svd'
   readf,1,flagtot
   close,1
   ;
   ;
   radiotz=dblarr(nt,npol,nz)
   flagtotz=dblarr(nt,npol,nz)
   Da_z,(freq_rest/freqtot-1.),distslice
   distbin=(floor((distslice-d0)/dd)) ;d0 is at bin0                               
   for it=0,nt-1 do begin
      for ipol=0,npol-1 do begin
         flagslice=reform(flagtot(it,ipol,*))
         for iz=0,nz-1 do begin
            ind=where(distbin eq iz and flagslice eq 1,cind)
            if cind gt 0 then begin
               radiotz(it,ipol,iz)=total(radiotf(it,ipol,ind))/double(cind)
               flagtotz(it,ipol,iz)=cind
            endif
         endfor
      endfor
   endfor
   flagtot=flagtotz
   ;
   openw,1,fdir+field+'simradio_tz_b4svd_'+signal
   printf,1,radiotz
   close,1
   ;
   openw,1,fdir+field+'flagtot_tz_b4svd_'+signal
   printf,1,flagtot
   close,1
   ;
   print,'finish reading'

   ;
endif else begin
   ;
   radiotz=dblarr(nt,npol,nz)
   flagtotz=dblarr(nt,npol,nz)
   openr,1,fdir+'simradio_tz_b4svd_'+signal
   readf,1,radiotz
   close,1
   ;
   flagtot=make_array(ntime*nscan,npol,nfreq*3/5*nwin,value=1.)
   openr,1,fdir+'flagtot_tz_b4svd'
   readf,1,flagtot
   close,1

endelse

   ;
   for imode=0,nmode-1 do begin
      ;
      timevec=dblarr(nz,nt,npol)  
      if (mode(imode) lt 5 or mode(imode) eq 6) then begin
         fmmdir=subfmdir(ifile)  
         openr,1,fmmdir+field+'.freqvec_'+smode(imode)
      endif else begin
         fmmdir=fmdir 
         if (mode(imode) gt 12 or mode(imode) eq 10 or mode(imode) eq 5) then begin
            openr,1,fmmdir+'freqvec_'+smode(imode) 
         endif else openr,1,fmmdir+field+'.freqvec_'+smode(imode)
      endelse
      ;
      readf,1,timevec
      close,1
      ;
      ;
   ; calculate the projected sim values using the data time-eigenmodes
      radiotmp=dblarr(nt,npol,nz)
      difftmp=radiotmp
;      for ipol=0,npol-1 do begin
      for ipol=0,0 do begin
         iradio=reform(radiotz(*,ipol,*))
;         la_svd,iradio,s,u,v,/double,status=status
         uvec=reform(timevec(*,*,ipol))
;         help,s
;         help,u
;         help,v
;         uvec=v
         ;
         matrix=iradio ## (uvec)
         help,matrix
         ;
         nm=n_elements(matrix(*,0))
         sm=dblarr(nm)
         vvec=dblarr(nm,nm)
         for im=0,nm-1 do begin
            vec=reform(matrix(im,*))
            sm(im)=sqrt(vec ## transpose(vec))
            if sm(im) gt 0 then vvec(im,*)=vec/sm(im)
         endfor
         ;
;         help,sm
;         print,'sm:',sm
;         help,s
;         print,'s',s
;         ;
;         print,'s-sm',s-sm
;         read,test
 ;        ;
 ;        help,vvec
;         print,'vvec:',vvec(0,*)
;         print,''
;         print,'u',u(0,*)
;         print,''
;         read,test
;         print,'vvec:',vvec(*,0)
;         print,''
;         print,'u',u(*,0)
;         print,''
         ;
;         read,test
         ;
         ss=dblarr(nm)
         ss(0:mode(imode)-1)=1.
         ;
         sm2=make_array(nm,value=1.)
;         ss(0:mode(imode)-1)=sm(0:mode(imode)-1)
         origin0=uvec ## transpose(uvec)
      
         help,origin0
         origin01=iradio ## uvec ## transpose(uvec)
         origin02=iradio ## origin0
         help,origin01
         help,origin02
         origin=matrix ## transpose(uvec)
         origin3=matrix ## diag_matrix(sm2) ## transpose(uvec)
         origin2=vvec ## diag_matrix(sm) ## transpose(uvec)
         help,origin
         help,iradio
         diff=origin-iradio
         ind=where(diff ne 0.,cind)
         print,'not zero elements',cind
         print,'max ,min diff',max(diff),min(diff),mean(diff),median(diff)
         print,'max, min origin',max(origin),min(origin),mean(origin),median(origin)
         print,'max, min iradio',max(iradio),min(iradio),mean(iradio),median(iradio)
         if ipol eq 0 then begin
            plt_image,transpose(iradio),/scalable,/colbar
      ;      plt_image,transpose(origin0),/scalable,/colbar
;            plt_image,transpose(origin01),/scalable,/colbar
;            plt_image,transpose(origin02),/scalable,/colbar
            plt_image,transpose(origin),/scalable,/colbar
            plt_image,transpose(origin2),/scalable,/colbar
;            plt_image,transpose(origin3),/scalable,/colbar
            plt_image,transpose(diff),/scalable,/colbar
         endif
         if mode(imode) gt 0 then sm2(0:mode(imode)-1)=0.
         residual=matrix ## diag_matrix(sm2) ## transpose(uvec)
         if ipol eq 0 then plt_image,transpose(residual),/scalable,/colbar
          ;
         diff=matrix ## diag_matrix(ss) ## transpose(uvec)
         radiotmp(*,ipol,*)=residual
         difftmp(*,ipol,*)=diff
      endfor
      ;
      ; rebin
      radio=dblarr(npol,nx,ny,nz)
      weight=radio
      noiseradio=radio
      radiodiff=radio
      weightdiff=weight
      mk_hicube_half,radiotmp,flagtot,flagtot,ratot,dectot,freqtot,pol,$
                      parallatic,radio=radio,weight=weight,$
                      field=field,noiseradio=noiseradio
      mk_hicube_half,difftmp,flagtot,flagtot,ratot,dectot,freqtot,pol,$
                     parallatic,radio=radiodiff,weight=weightdiff,$
                     field=field,noiseradio=noiseradio
      ;
      radio2=radio
      ind=where(weight gt 0.,cind)
      radio2(ind)=radio2(ind)/weight(ind)
      ;
      openw,1,fmdir+field+'.svdresidual_weighted_'+smode(imode)
      printf,1,radio2
      close,1
      ;
      openw,1,fmdir+field+'.svdresidual_weight_'+smode(imode)
      printf,1,weight
      close,1
      ;
      radiotot(imode,*,*,*,*)=radiotot(imode,*,*,*,*)+radio
      weighttot(imode,*,*,*,*)=weighttot(imode,*,*,*,*)+weight
      ;
      radio2=radiodiff
      ind=where(weightdiff gt 0.,cind)
      radio2(ind)=radio2(ind)/weightdiff(ind)
      ;
      openw,1,fmdir+field+'.svddiff_weighted_'+smode(imode)
      printf,1,radio2
      close,1
   endfor
   device,/close
endfor


;radio2=radiotot
ind=where(weighttot gt 0.,cind)
radiotot(ind)=radiotot(ind)/weighttot(ind)
;
;device,filename=field+'.residual.ps'
;
for imode=0,nmode-1 do begin
   openw,1,fmdir+field+'.svdresidualtot_weighted_'+smode(imode)
   printf,1,reform(radiotot(imode,*,*,*))
   close,1
   ;
   openw,1,fmdir+field+'.svdresidualtot_weight_'+smode(imode)
   printf,1,reform(weighttot(imode,*,*,*))
   close,1
   ;
endfor

device,/close
!p.multi=0   
   
end
