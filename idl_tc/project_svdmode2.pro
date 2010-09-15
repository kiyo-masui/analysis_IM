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
                                ; calculate the projected sim values
                                ; using the data time-eigenmodes         
   radio=dblarr(nt,npol,nz)
   for ipol=0,npol-1 do begin
      iradio=reform(radiotz(*,ipol,*))
      uvec=timevec(ipol*ntime:(ipol+1)*ntime-1)
      ;
      matrix=iradio ## transpose(uvec)
      help,matrix
      ;
      nm=n_elements(matrix(*,0))
      sm=dblarr(nm)
      vvec=dblarr(nm,nm)
      for im=0,nm-1 do begin
         vec=refrom(matrix(*,im))
         sm(im)=sqrt(vec ## transpose(vec))
         vvec(*,im)=vec/sm(im)
      endfor
      ;
      residual=
   endfor
      
