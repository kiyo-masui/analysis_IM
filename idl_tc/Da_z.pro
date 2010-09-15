function dadz,z
  return,1./sqrt(0.24*(1.+z)^3.+0.76)
end


pro Da_z,z,da
; given a redshift calculate the comoving angular diameter distance 
;
; parameter
c0=2.9998e5
Ho=100.
c0Ho=c0/Ho

da=qromb('dadz',0,z,/double)*c0Ho

end
