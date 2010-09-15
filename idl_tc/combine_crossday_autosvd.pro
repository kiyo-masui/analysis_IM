pro combine_crossday_autosvd,field=field

;if not keyword_set(field) then field='f4'

;
;if field eq 'f3' then begin

; field3
openw,5,'field3auto.crossday.10mode+svd.dat'
openw,6,'field3auto.crossday.10mode+svd.signal.dat'
openw,7,'field3auto.crossday.10mode+svd.err.dat'
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode'
printf,5,'10-0mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='1'
printf,5,'10-1mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='2'
printf,5,'10-2mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='3'
printf,5,'10-3mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='4'
printf,5,'10-4mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,'''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='5'
printf,5,'10-5mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='6'
printf,5,'10-6mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='8'
printf,5,'10-8mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
;
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='10'
printf,5,'10-10mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f3',mode='10mode',/svdfg,svdmode='12'
printf,5,'10-12mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;correlate_auto_crossday,auto,meanauto,field='f3',mode='10mode',/dosigcut,/svdfg,svdmode='8'
;printf,5,'10-8mode'
;printf,5,'corr:',auto
;printf,5,'mean:',meanauto
;printf,5,'mean Kelvin:',sqrt(abs(meanauto))
;printf,5,''
;correlate_auto_crossday,auto,meanauto,field='f3',mode='10mode',/dosigcut,/svdfg,svdmode='10'
;printf,5,'10-10mode'
;printf,5,'corr:',auto
;printf,5,'mean:',meanauto
;printf,5,'mean Kelvin:',sqrt(abs(meanauto))
;printf,5,''
;correlate_auto_crossday,auto,meanauto,field='f3',mode='10mode',/dosigcut,/svdfg,svdmode='12'
;printf,5,'10-12mode'
;printf,5,'corr:',auto
;printf,5,'mean:',meanauto
;printf,5,'mean Kelvin:',sqrt(abs(meanauto))
;printf,5,''
;
close,5
close,6
close,7
;
;
;
;if field eq 'f4' then begin

openw,5,'field4auto.crossday.10mode+svd.dat'
openw,6,'field4auto.crossday.10mode+svd.signal.dat'
openw,7,'field4auto.crossday.10mode+svd.err.dat'
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode'
printf,5,'10-0mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='1'
printf,5,'10-1mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='2'
printf,5,'10-2mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='3'
printf,5,'10-3mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='4'
printf,5,'10-4mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='5'
printf,5,'10-5mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='6'
printf,5,'10-6mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='8'
printf,5,'10-8mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='10'
printf,5,'10-10mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
correlate_auto_crossday,auto,meanauto,err,field='f4',mode='10mode',/svdfg,svdmode='12'
printf,5,'10-12mode'
printf,5,'corr:',auto
printf,5,'mean:',meanauto
printf,5,'mean Kelvin:',sqrt(abs(meanauto))
printf,5,''
printf,6,sqrt(abs(meanauto))
printf,7,err
;
close,5
close,6
close,7
;
;



end
