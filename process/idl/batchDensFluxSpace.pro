 pro batchDensFluxSpace, sstart, send, step

; File Options
;directory = "/data1/r900-1/ctss/FLASH_3.2/TestCosmo/"
;directory = "/scratch/copland/d1/ctss/FLASH_storage/TestCosmo2/"

;directory = "/ranger/scratch/01707/mepa/Final"
;directory = "/ranger/scratch/01707/mepa/Fry"
;directory = "/ranger/scratch/01707/mepa/Rad"
;directory = "/ranger/scratch/01707/mepa/FryRad"
directory = "/ranger/scratch/01707/mepa/FryRad062"
;directory = "/data1/r900-3/mepa/Work/FryRad062"
;file      = "/1star_hdf5_chk_"
file      = "/final_hdf5_chk_"
;num       = 0022
;num       = 0007
;sstart    = 1
;send      = 10
;step      = 1
;box_size  = 5.0*3.08e24
box_size  = 2.0*3.08e24
print, file

centerx = box_size / 2.0
centery = centerx
centerz = centerx

;radius = 0.5*3.08e24
radius = 3.08e24

xmin = centerx - radius
xmax = centerx + radius
ymin = centery - radius
ymax = centery + radius
zmin = centerz - radius
zmax = centerz + radius

;bhmass    = 1.35D34

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



;number = 42
;if (sstart gt send and step gt 0) then send = sstart
for number = sstart,send,step do begin

    if (number ge 1)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''

    filename = directory + file + prefix + String(strcompress(number,/remove))
    print, filename
    ;outfile = 'phase_' + String(strcompress(number, /remove)) + '.png'
    ;outfile = 'phase.png'
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
    outfile = 'flux_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    

    print, "reading density..."
    ;dens =  loaddata_nomerge(filename,'dens', XCOORDS=x, YCOORDS=y, ZCOORDS=z)
    dens = loaddata_nomerge(filename,'dens')
    flux = loaddata_nomerge(filename,'rad1')
    ;c011 =  loaddata_nomerge(filename,'c011')
    ;c012 =  loaddata_nomerge(filename,'c012')
    ;c013 =  loaddata_nomerge(filename,'c013')
    ;c014 =  loaddata_nomerge(filename,'c014')
    ;c015 =  loaddata_nomerge(filename,'c015')
    ;c016 =  loaddata_nomerge(filename,'c016')
    ;temperature = loaddata_nomerge(filename,'temp')

    print, "***"
    print, "size of dens struct"
    print, size(dens)
    print, "***"

    ;lrefine = getRefineLevel(filename, 'dens')
    ;ndim = determine_file_dimensionality(filename)
    ;print, 'ndim is: ', ndim
    ;get_lrefine_max_min, filename, ndim, MAX=maxLevels, MIN=minLevels
    ;lrefine = get_lrefine_max_min(filename, ndim);, MIN=lrefinemin, MAX=lrefinemax)
    ;print, 'lrefinemin is: ', minLevels
    ;print, 'lrefinemax is: ', maxLevels
    ;print, 'lrefine is: ', lrefine

    ;numCells = 2.0^(maxLevels + 3)
    ;print, 'numcells is: ', numCells
    ;dx = box_size / oneplusred / numCells
    ;print, 'dx is: ', dx

    ;print, size(temperature)
    ;print, "***"
    ;print, "size of lrefine struct"
    ;print, size(lrefine)
    
    ;densSize = size(dens)
    
    ;totalBlocks = densSize[1]
    
    ;nxb = densSize[2]
    ;nyb = densSize[3]
    ;nzb = densSize[4]

    ;temperature = temperature/oneplusred^2.0

    dens = dens * oneplusred^3.0
    numDens = dens / 1.67e-24

    ;flux = flux * oneplusred^2.0

    ;tauCellH = numDens * sigmaH * dx
    ;tauCellHe = numDens * sigmaHe * dx
    ;tauCellHep = numDens * sigmaHep * dx
    ;tauCell = tauCellH + tauCellHe + tauCellHep

    ;tau = 

    absorbedFlux = 6.0e6 - flux

    xr = [1e-2, 1e0]
    yr = [1e3, 1e8]

    plot, numDens, flux, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='number density (cm^-3)', ytitle='flux (cm^-2 s^-1)'

    oplot, numDens, absorbedFlux, color=240, psym=3


    ;plot, numDens, flux, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xtitle='number density (cm^-3)', ytitle='flux (cm^-2 s^-1)'
    
     void = cgSnapshot(FILENAME=outfile, /NoDialog)
endfor

end
