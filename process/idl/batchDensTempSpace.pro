 pro batchDensTempSpace, sstart, send, step

; File Options
directory = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128"
file      = "/radCosmoLW_hdf5_chk_"

box_size  = 1.0*3.08e24
print, file

centerx = box_size / 2.0
centery = centerx
centerz = centerx

radius = 0.5*3.08e24
;radius = 3.08e24

xmin = centerx - radius
xmax = centerx + radius
ymin = centery - radius
ymax = centery + radius
zmin = centerz - radius
zmax = centerz + radius

; Output Options
plotps = 0
plotgif = 1
charsize = 1.4

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
    outfile = 'plots/phase/phase_' + String(strcompress(number, /remove)) + '_z' + String(strcompress(redshift, /remove)) + '.png'    

    print, "reading density..."
    dens =  loaddata_nomerge(filename,'dens', XCOORDS=x, YCOORDS=y, ZCOORDS=z)
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
   
    temperature = temperature/oneplusred^2.0

    numdens = dens/1.67e-24

    ;xr = [1e-4, 2e5]
    ;yr = [1e0, 1e5]

    ;xr = [1e-4, 1e4]
    ;yr = [1e0, 1e5]

    xr = [1e-5, 1e5]
    yr = [1e-1, 1e5]

    ;xr = [1e-4, 1e5]
    ;yr = [3e0, 5e6]

    ;xr = [1e-2, 1e2]
    ;yr = [1e2, 1e3]

    plot, numdens, temperature, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1, xtitle='physical number density (cm^-3)', ytitle='temperature (K)'
    
    ; line of adiabatic collapse
    n1 = findgen(100)
    ;n1 = n1/100
    ;n1 = n1/2
    ;oplot, n1, 10.0 * (n1/0.1)^(2.0/3.0), color=0
    ;oplot, n1/2000, 10.0 * (n1*100)^(2.0/3.0), color=0
    oplot, n1/1000, 10.0 * (n1/10)^(2.0/3.0), color=0

    ;; Jeans floor    
    ;n2 = findgen(50)
    ;oplot, n2*10, n2*10, color=2  
    ;;oplot, n2*2000, n2*2, color=2 
 
    ; virial temperature
    mu = 1.2
    ;n3 = findgen(100)
    ;T3 = identity(100)
    ;;Mvir = 4.0e11 * exp(-0.5 * redshift) ; 1.4*sigma8 run
    ;Mvir = 1.0e10 * exp(-0.5 * redshift) ; 1.0*sigma8 run **CHECK**
    Mvir = 6.0e10 * exp(-0.5 * redshift) ; 1 Mpc box (maybe not correct scaling)
    Tvir = 2.0e4 * (mu / 1.2) * (Mvir / 1.0e8)^(2.0/3.0) * (1.0 + redshift) / 10.0
    oplot, xr, [Tvir, Tvir], color=3

    print, Mvir, Tvir

    ; CMB temperature
    Tcmb = 2.2725 * oneplusred
    oplot, xr, [Tcmb, Tcmb], color=4

    ;mytemp = fltarr(100)

    ;dens1 = findgen(100) * 3.0e-25
    ;temperature1 = dens1^(2.0/3.0) * 1.0e14

    ;oplot, dens1, temperature1, color=0

    ;print, dens1
    ;print, temperature1
    
     void = cgSnapshot(FILENAME=outfile, /NoDialog)
endfor

end
