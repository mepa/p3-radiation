 pro batchHpluDensTempSpace, sstart, send, step

; File Options
;directory = "/data1/r900-1/ctss/FLASH_3.2/TestCosmo/"
;directory = "/scratch/copland/d1/ctss/FLASH_storage/TestCosmo2/"

;directory = "/ranger/scratch/01707/mepa/Final"
;directory = "/ranger/scratch/01707/mepa/Fry"
directory = "/ranger/scratch/01707/mepa/Rad"
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
    outfile = 'hplu_phase_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    

    print, "reading density..."
    dens =  loaddata_nomerge(filename,'dens', XCOORDS=x, YCOORDS=y, ZCOORDS=z)
    hplufrac =  loaddata_nomerge(filename,'hplu', XCOORDS=x, YCOORDS=y, ZCOORDS=z)
    ;dens = loaddata_nomerge(filename,'dens')
    temperature = loaddata_nomerge(filename,'temp')
    
    ; variable name does not matter for xyz load coords
    ;print, "reading x..."
    ;x = load_coords(filename, 'dens', 1)
    ;print, "reading y..."
    ;y = load_coords(filename, 'dens', 2)
    ;print, "reading z..."
    ;z = load_coords(filename, 'dens', 3)
    
    
    ;list = where((x gt xmin) and (x lt xmax) and (y gt ymin) and (y lt ymax) and (z gt zmin) and (z lt zmax))
    
    ;dens = dens(list)
    ;temperature = temperature(list)

    ;lrefine = getRefineLevel(filename, 'dens')
    
         ; Let's compute the enclosed baryonic
         ; mass within radius of x_center, y_center, and z_center

    massEnclosed = 0.0

    print, "***"
    print, "size of dens struct"
    print, size(dens)
    print, "***"
    ;print, size(temperature)
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
    
    ;get to physical coordinates!
    ;;x = x / oneplusred
    ;y = y / oneplusred
    ;z = z / oneplusred
   

    ; Let's find out what the enclosed mass radially around a given 
    ; point it. Note: does not include mass on other side of grid
    ; given periodic BCs

    
  
 
   
    temperature = temperature/oneplusred^2.0

    numdens = hplufrac*dens/1.67e-24

    xr = [1e-4, 1e4]
    yr = [1e0, 1e5]

    plot, numdens, temperature, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1
    
    n = findgen(100)
    
    n = n/100
    n= n/2
    
    ; line of adiabatic collapse
    oplot, n, 10.0 * (n/0.1)^(2.0/3.0), color=0
    
    n = n * 1000   
    ;oplot, n, 3.0 * n, color=2
    oplot, n, 10.0 * n, color=2
  
    
    mytemp = fltarr(100)

    ;dens1 = findgen(100) * 3.0e-25
    ;temperature1 = dens1^(2.0/3.0) * 1.0e14

    ;oplot, dens1, temperature1, color=0

    ;print, dens1
    ;print, temperature1
    
     void = cgSnapshot(FILENAME=outfile, /NoDialog)
endfor

end
