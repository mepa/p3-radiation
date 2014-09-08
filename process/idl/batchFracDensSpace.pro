 pro batchFracDensSpace, sstart, send, step

; File Options
;directory = "/work/01707/mepa/Rad"
;directory = "/work/01707/mepa/Rad_1.0sigma8"
;directory = "/work/01707/mepa/Rad_1.0sigma8/newChem"
;directory = "/scratch/01707/mepa/Rad_1.0sigma8"
directory = "/scratch/01707/mepa/Rad_res512"
file      = "/rad_hdf5_chk_"

;directory = "/data1/r900-3/mepa/SUBFIND/120725/Flash2Gadget_Shell_Script"
;file      = "/final_hdf5_chk_"

print, file

; Output Options
plotps = 0
plotgif = 1
charsize = 1.4
;if plotps eq 1 then begin
;    device,xsize=16.0
;    device,ysize=12.0
;end

;if plotgif eq 1 then begin
 ;   window, 10, xsize=800,ysize=600
;end

;!P.Multi = [0,2,2,0,0]
;loadct,39

for number = sstart,send,step do begin

    if (number ge 1)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''

    filename = directory + file + prefix + String(strcompress(number,/remove))
    print, filename

    ;;;;;
    ; get redshift of file
    file_identifier = H5F_OPEN(filename)
    dataset = H5D_OPEN(file_identifier, "real scalars")
    real_scalars = H5D_READ(dataset)
    for i=1, (size(real_scalars))[3] do begin
        if (stregex(real_scalars[i-1].name, '^scale', /BOOLEAN)) then scale = real_scalars[i-1].value
    endfor
    ;;;;;
    
    redshift = (1.0 / scale) - 1.0
    oneplusred = 1.0 + redshift
    
    print, "redshift = ", redshift

    outfileH = 'h_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    
    outfileHplu = 'hplu_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    
    outfileHtwo = 'htwo_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    

    print, "reading density..."
    dens = loaddata_nomerge(filename,'dens')
    hfrac =  loaddata_nomerge(filename,'h')
    hplufrac =  loaddata_nomerge(filename,'hplu')
    htwofrac =  loaddata_nomerge(filename,'htwo')
    
    print, "***"
    print, "size of dens struct"
    print, size(dens)
    print, "***"


    dens = dens * oneplusred^3.0
    
    numdens = dens/1.67e-24

    print, dens(10), hfrac(10), hplufrac(10), htwofrac(10)


    ;;; HYDROGEN FRACTION ;;;
    xr = [1e-4, 1e4]
    yr = [0.5e-5, 1]
    ;yr = [0.755, 0.76]

    plot, numdens, hfrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='H abundance'

    ;plot, numdens, hfrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr,  xstyle=1, ystyle=1

    void = cgSnapshot(FILENAME=outfileH, /NoDialog)


    ;;; HPLUS FRACTION ;;;
    
    xr = [1e-4, 1e4]
    yr = [1e-6, 1]
    ;yr = [1e-5, 1e-3]

    plot, numdens, hplufrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='Hplus abundance'

    ;plot, numdens, hplufrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr,  xstyle=1, ystyle=1

    void = cgSnapshot(FILENAME=outfileHplu, /NoDialog)


    ;;; HTWO FRACTION ;;;

    xr = [1e-4, 1e4]
    ;yr = [1e-14, 1e-2] ; used for Fry
    yr = [1e-18, 1e-1] ; used for RadTrans
    ;yr = [1e-6, 1e-3] ; used for Rad

    plot, numdens, htwofrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='Htwo abundance'

    ;plot, numdens, htwofrac, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr,  xstyle=1, ystyle=1

    void = cgSnapshot(FILENAME=outfileHtwo, /NoDialog)
endfor

end
