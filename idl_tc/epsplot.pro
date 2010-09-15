; davor window-Funktion ausschalten
PRO EpsPlot, File, XSize = SizeX, YSize = SizeY, font_size=chsz
 print, "Preparing to plot to file '"+File+"'..."
 on_error,1

;  setcolors
; if n_elements(chsz) le 0 then chsz = 12
 set_plot, "PS"
 device, filename=File, /encapsulated, /color, bits_per_pixel=8
; device,/times,font_index=20,/italic
; device, font_size=chsz
;  device, yoffset=2, xoffset=0.5

; !p.font  = 0  ;; use fonts supplied by the graphics device
; !p.thick = 2.5
; !p.multi=[1,1,1]

; !X.THICK = 2.0
; !Y.THICK = 2.0

; if ( keyword_set( SizeX ) eq 0 ) then SizeX = 16.0  ; in cm
; if ( keyword_set( SizeY ) eq 0 ) then SizeY = 12.0
; A4: xs=29.7cm, ys=21cm

; device, xsize = SizeX, ysize = SizeY

END
