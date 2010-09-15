pro plt_svdmode,field=field,sim=sim

; svd dimension time vector = [z,t]=[560, 2154] = [560, 2240]

if not keyword_set(field) then field='f3'
if not keyword_set(signal) then signal='100muK'

if field eq 'f3' then begin
   ;
;   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_simonly_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'
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
   ;
endif


if field eq 'f4' then begin
   ;
   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed/'
;   if keyword_set(sim) then thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode_sim_'+signal+'/' else thisdir='/cita/d/raid-cita/tchang/GBT/Processed_mode/'
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
   ;
endif
set_plot,'ps'
device,/color, bits_per_pixel=8

; parameters
nfile=n_elements(submapdir)
nfile=1
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
device,filename=field+'.timevec.ps'
;for ifile=0,nfile-1 do begin
for ifile=8,8 do begin
   ;
;   if (dname eq 'sept07' and ifile eq 3) then nscan=3
;   if (dname eq 'sept06' and ifile eq 4) then nscan=4
   ;fdir='/cita/d/raid-cita/tchang/GBT/Processed/'+dname+'/'+field+'_'+ffdir(ifile)+'/'
   fdir=submapdir(ifile)
   if fdir eq '/cita/d/raid-cita/tchang/GBT/Processed/sept07/f3_27-29/' then nscan=3
   if fdir eq '/cita/d/raid-cita/tchang/GBT/Processed/sept06/f4_33-36/' then nscan=4
   nt=ntime*nscan
   ;
   timevec=dblarr(nz,nt)   
   openr,1,fdir+field+'.freqvec_1mode'
   readf,1,timevec
   close,1
   ;
   for imode=0,9 do begin
      plot,reform(timevec(imode,*))
   endfor
   ;
;   timevec=dblarr(nz,nt)   
;   openr,1,fdir+field+'.freqvec_2mode'
;   readf,1,timevec
;   close,1
   ;
;   for imode=0,6 do begin
;      plot,reform(timevec(imode,*))
;   endfor
   ;
;   timevec=dblarr(nz,nt)   
;   openr,1,fdir+field+'.freqvec_6mode'
;   readf,1,timevec
;   close,1
;   ;
;   for imode=0,6 do begin
;      plot,reform(timevec(imode,*))
;   endfor

endfor
device,/close
!p.multi=0   
   
end
