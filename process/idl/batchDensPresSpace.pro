 pro batchDensPresSpace, sstart, send, step

; File Options
;directory = "/work/01707/mepa/Rad"
;directory = "/work/01707/mepa/Rad_1.0sigma8"
;directory = "/work/01707/mepa/Rad_1.0sigma8/zInitial_149.4"
;directory = "/work/01707/mepa/Rad_1.0sigma8/newChem"
;directory = "/scratch/01707/mepa/zInitial_74.2/eintSwitch_1.0"
;directory = "/scratch/01707/mepa/Rad_1.0sigma8"
directory = "/scratch/01707/mepa/Rad_res512"
file      = "/rad_hdf5_chk_"
;file      = "/final_hdf5_chk_"
;num       = 0022
;num       = 0007
;sstart    = 1
;send      = 10
;step      = 1

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


;if (sstart gt send and step gt 0) then send = sstart
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
    outfile = 'pres_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    

    print, "reading density..."
    dens =  loaddata_nomerge(filename,'dens', XCOORDS=x, YCOORDS=y, ZCOORDS=z)
    ;dens = loaddata_nomerge(filename,'dens')
    pres = loaddata_nomerge(filename,'pres')
    
    print, "***"
    print, "size of dens struct"
    print, size(dens)
    print, "***"
    ;print, size(pres)
    print, "***"
    print, "size of lrefine struct"
    print, size(lrefine)
    print, "***"
    print, "size of x coord struct"
    print, size(x)
    print, "size of y coord struct"
    print, size(y)
    print, "size of z coord struct"
    print, size(z)
    print, "***"
    
    densSize = size(dens)
    
    totalBlocks = densSize[1]
    
    nxb = densSize[2]
    nyb = densSize[3]
    nzb = densSize[4]

    dens = dens * oneplusred^3.0
    
    pres = pres * oneplusred

    numdens = dens/1.67e-24

    xr = [1e-4, 2e5]
    ;yr = [1e0, 1e5]

    ;xr = [1e-4, 1e4]
    ;yr = [1e0, 1e5]

    ;xr = [1e-4, 1e5]
    ;yr = [5e-1, 1e5]

    ;xr = [1e-4, 1e5]
    ;yr = [3e0, 5e6]

    ;xr = [1e-2, 1e2]
    ;yr = [1e2, 1e3]

    ;plot, numdens, preserature, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='preserature (K)'
    plot, numdens, pres, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='pressure (g s^-2 cm^-1)'
    
     void = cgSnapshot(FILENAME=outfile, /NoDialog)
endfor

end
