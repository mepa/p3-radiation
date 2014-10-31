pro maxDensRedshift, sstart, send, step

; File Options
directory = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128"
file      = "/radCosmoLW_hdf5_chk_"

pi = 3.1415927e0
grav_const = 6.6726e-8

hubble_const = 2.27502706e-18 ; 70.2 km/s/Mpc
omega_m = 0.275e0
omega_v = 0.725e0
omega_b = 0.0458

pfname = 'plots/maxDens_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
openw,lun,pfname,/get_lun

outfile = 'plots/maxdens_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'

num_pts = floor((send-sstart)/step)+1
maxdensities = dblarr(num_pts) 
mean_baryon_densities = dblarr(num_pts) 
redshifts = dblarr(num_pts) 

count = 0L
for number = sstart,send,step do begin

    if (number ge 1)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''
    
    filename = directory + file + prefix + String(strcompress(number,/remove))
    print, filename

    if (number eq 228) then continue ;TACC lost this file

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
    redshifts[count] = redshift

    ;;;;;
    ; get tree
    ;print, "start getting tree..."
    read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params

    mass_enclosed_baryon = 0.0

    max_x = 0.0
    max_y = 0.0
    max_z = 0.0
    
    max_value = 0.0

    ; ENCLOSED BARYONIC MASS
    ; loop over blocks
    for block_id = 0L, params.totBlocks - 1L, 1L do begin
    
        if (tree[block_id].nodetype EQ 1) then begin
            
            dx = tree[block_id].size[0] / params.nxb
            dy = tree[block_id].size[1] / params.nyb
            dz = tree[block_id].size[2] / params.nzb
            
            
            block_x_min = tree[block_id].bndBox[0,0]
            block_y_min = tree[block_id].bndBox[0,1]
            block_z_min = tree[block_id].bndBox[0,2]
            
            block_x_max = tree[block_id].bndBox[1,0]
            block_y_max = tree[block_id].bndBox[1,1]
            block_z_max = tree[block_id].bndBox[1,2]
            
                        
            ; loop over cells
            for ii=0, params.nxb-1 do begin
                for jj=0, params.nyb-1 do begin
                    for kk=0, params.nzb-1 do begin
                        
                        x_cell = block_x_min + (ii + 0.5) * dx
                        y_cell = block_y_min + (jj + 0.5) * dy
                        z_cell = block_z_min + (kk + 0.5) * dz
                        
                        if dens(0,block_id,ii,jj,kk) gt max_value then begin
                            
                            max_value = dens(0,block_id,ii,jj,kk)
                            max_x = x_cell
                            max_y = y_cell
                            max_z = z_cell
                            
                        endif
                        
                    endfor
                endfor    
            endfor
        endif      
    endfor

    max_dens_physical = max_value * oneplusred^3.0
    maxdensities[count] = max_dens_physical

    hubble_param_sq = hubble_const^2.0 * (omega_m *oneplusred^3.0 + omega_v)
    
    rho_crit = 3.0 * hubble_param_sq / (8.0 * pi * grav_const) / oneplusred^3.0
    rho_background = rho_crit * omega_m

    rho_crit_approx = 3.0 * hubble_const^2.0 / (8.0 * pi * grav_const)
    rho_background_approx = rho_crit_approx * omega_b
    mean_baryon_densities[count] = rho_background_approx * oneplusred^3.0
    
    print, "max density value and location"
    print, max_value, max_x, max_y, max_z
    print, "redshift"
    print, redshift
    print, "physical max density"
    print, max_dens_physical
    print, "mean baryon density"
    print, mean_baryon_densities[count]
    print, "background density"
    print, rho_background

    printf, lun, FORMAT='(I,F,E,E,E,E)', number, redshift, max_dens_physical, max_x, max_y, max_z

    count = count + 1
endfor

close, lun
free_lun, lun

;xr = [1e-4, 1e4]
;yr = [1e-3, 1e5]

;plot, redshifts, maxdensities, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1

;plot, redshifts, maxdensities, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xtitle='redshift', ytitle='phys max dens (g cm^-3)'

plot, redshifts, maxdensities/1.67e-24, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xtitle='redshift', ytitle='phys max num dens (cm^-3)'

oplot, redshifts, mean_baryon_densities/1.67e-24, color=1
    
void = cgSnapshot(FILENAME=outfile, /NoDialog)

end
