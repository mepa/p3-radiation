pro m200Redshift, sstart, send, step

; File Options
;directory = "/data1/r900-4/mepa/Radiation/Work/Fry"
directory = "/scratch/cerberus/d4/mepa/stampede/Rad"
file      = "/rad_hdf5_chk_"

;num_bins = 10000

;;comoving r_min/r_max
;r_min = 1.0d18
;r_max = 1.0d21

num_bins = 100

;comoving r_min/r_max
r_min = 1.0d21
r_max = 1.0d24
    
pi = 3.1415927
h0 = 2.27502706e-18 ; 70.2 km/s/Mpc
G =  6.6726e-8
omega_m = 0.275
omega_v = 0.725

; Output Options
plotps = 0
plotgif = 1
charsize = 1.4

pfname1 = 'halo_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
openw,lun1,pfname1,/get_lun

for number = sstart,send,step do begin

    if (number ge 1)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''

    filename = directory + file + prefix + String(strcompress(number,/remove))
    

    print, 'FILE: ', filename

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Read input data

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ;;;; scalefactor and redshift

    file_identifier = H5F_OPEN(filename)
    dataset = H5D_OPEN(file_identifier, "real scalars")
    real_scalars = H5D_READ(dataset)
    for i=1, (size(real_scalars))[3] do begin
        if (stregex(real_scalars[i-1].name, '^scalefactor', /BOOLEAN)) then scale = real_scalars[i-1].value
        if (stregex(real_scalars[i-1].name, '^time', /BOOLEAN)) then time = real_scalars[i-1].value
    endfor
    
    redshift = 1.0/scale - 1.0
    oneplusred = 1.0 + redshift

    print, redshift

    hubble_param_sq = h0^2 * (omega_m * oneplusred^3.0 + omega_v)
    rho_crit = 3.0 * hubble_param_sq / 8.0 / pi / G / oneplusred^3.0
    density_thresh = 200.0 * rho_crit
    
    rho_crit_approx = 3.0 * h0^2.0 / 8.0 / pi / G 
    rho_background_approx = rho_crit_approx * omega_m
    density_thresh_approx = 200.0 * rho_background_approx
    
    print, density_thresh, density_thresh_approx


    
    ;;;;;;;;;;;;;;;;;;;;;
    ;;;; Grid data

    read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='gpot', TREE=tree, DATA=gpot, PARAMETERS=params
    
;    read_amr_no_part, filename, VAR_NAME='pde', TREE=tree, DATA=pde, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=pden, PARAMETERS=params


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Halo Center


    min_spacing = 1.0e55

    max_value = 0.0
    min_value = 0.0

       for block_id = 0L, params.totBlocks - 1L, 1 do begin
        
        if (tree[block_id].nodetype EQ 1) then begin
            
            dx = tree[block_id].size[0] / params.nxb
            dy = tree[block_id].size[1] / params.nyb
            dz = tree[block_id].size[2] / params.nzb

            if (dx lt min_spacing) then begin
                min_spacing = dx
            endif
            
            
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
                        
                        
                        ; let's center it on minimum of grav potential!
                        if gpot(0,block_id,ii,jj,kk) lt min_value then begin
                            
                            min_value = gpot(0,block_id,ii,jj,kk)
                            max_x = x_cell
                            max_y = y_cell
                            max_z = z_cell
                            
                        endif
                        
                    endfor
                endfor    
            endfor
        endif      
     endfor

     x_center = max_x
     y_center = max_y
     z_center = max_z

     ;x_center = 3.39e24
     ;y_center = 3.11e24
     ;z_center = 3.12e24

     print, "GPOT Center: ", max_x, max_y, max_z


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Sum quantities into bins


    baryon_mass = dblarr(num_bins+1) 
    pden_mass = dblarr(num_bins+1)
    volume = dblarr(num_bins+1)
    total_mass = dblarr(num_bins+1)

    radius = dblarr(num_bins+1)

    sum_mass = 0.0
    sum_gas_mass = 0.0

    for block_id = 0L, params.totBlocks - 1L, 1L do begin
        
        if (tree[block_id].nodetype EQ 1) then begin
            
            dx = tree[block_id].size[0] / params.nxb ;/ oneplusred
            dy = tree[block_id].size[1] / params.nyb; / oneplusred
            dz = tree[block_id].size[2] / params.nzb; / oneplusred
            
            block_x_min = tree[block_id].bndBox[0,0]; / oneplusred
            block_y_min = tree[block_id].bndBox[0,1]; / oneplusred
            block_z_min = tree[block_id].bndBox[0,2]; / oneplusred
            
            ; loop over cells
            for ii=0, params.nxb-1 do begin
                for jj=0, params.nyb-1 do begin
                    for kk=0, params.nzb-1 do begin
                        
                        x_cell = block_x_min + (ii + 0.5) * dx
                        y_cell = block_y_min + (jj + 0.5) * dy
                        z_cell = block_z_min + (kk + 0.5) * dz

                        vol = dx*dy*dz
                        mass_cell = dens[0,block_id,ii,jj,kk] * vol
                        pden_mass_cell = pden[0,block_id,ii,jj,kk] * vol

                        sum_mass = sum_mass + mass_cell + pden_mass_cell
                        sum_gas_mass = sum_gas_mass + mass_cell


                        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                        ;;;; GPOT Centered

                        deltax = x_cell - max_x
                        deltay = y_cell - max_y
                        deltaz = z_cell - max_z
                        
    
                        distanceFromCenter = sqrt((deltax)^2.0  + $
                                                  (deltay)^2.0  + $ 
                                                  (deltaz)^2.0 )
                        

                        if(distanceFromCenter lt r_max) then begin
                            
                            if (distanceFromCenter lt r_min) then begin
                                bin = 0
                            endif else begin
                                bin = floor(double(num_bins) * alog(distanceFromCenter/r_min) / alog(r_max / r_min)) + 1
                            endelse
    
                            volume[bin] = volume[bin] + vol
     

                            baryon_mass[bin] = baryon_mass[bin] + mass_cell
                            pden_mass[bin] = pden_mass[bin] + pden_mass_cell
                            total_mass[bin] = total_mass[bin] + pden_mass_cell + mass_cell

                        endif

                    endfor
                endfor    
            endfor
        endif      
    endfor




;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Output radial profile



    pfname = 'profile_' + prefix + String(strcompress(number,/remove)) + '.txt'
    openw,lun,pfname,/get_lun

    scalecube = scale^3.0

    printf, lun, time, redshift

    for i=0, num_bins do begin
        radius[i] = r_min * (r_max/r_min)^(double(i) / double(num_bins))

        if (volume[i] EQ 0.0 or baryon_mass[i] EQ 0.0) then begin
          avg_density = 0.0
          bar_density = 0.0
          pden_density = 0.0
        endif else begin
          avg_density = total_mass[i]/volume[i]
          bar_density = baryon_mass[i]/volume[i]
          pden_density = pden_mass[i]/volume[i]
        endelse
        
        printf, lun, FORMAT='(E,E,E,E)', radius[i]*scale, volume[i]*scalecube, avg_density/scalecube, bar_density/scalecube, pden_density/scalecube
    endfor

    close, lun
    free_lun, lun


    
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Find the virial mass/radius

    
    
    ;; Sum interior bins for average enclosed properties
    for i=1, num_bins do begin
        total_mass[i] = total_mass[i-1] + total_mass[i]
        volume[i] = volume[i-1] + volume[i]
    endfor
    
    R200 = 0.0
    M200 = 0.0
    for i=0, num_bins do begin
        avg_density = total_mass[i]  / volume[i]
        
        if(avg_density lt density_thresh) then begin
            print, "R200: ", radius[i]*scale, ' cm   (', radius[i]*scale / 3.0856e18, ' parsec)'
            print, "M200: ", total_mass[i], ' g   (', total_mass[i] / 2e33, ' M_sun)'
            print, "Masses: ", sum_mass, sum_gas_mass;, sum_metal_mass

            R200 = radius[i]*scale
            M200 = total_mass[i]

            break
        endif
    endfor

    R200_approx = 0.0
    M200_approx = 0.0    
    for i=0, num_bins do begin
        avg_density = total_mass[i]  / volume[i]
        
        if(avg_density lt density_thresh_approx) then begin
            print, "R200_approx: ", radius[i]*scale, ' cm   (', radius[i]*scale / 3.0856e18, ' parsec)'
            print, "M200_approx: ", total_mass[i], ' g   (', total_mass[i] / 2e33, ' M_sun)'
            print, "Masses: ", sum_mass, sum_gas_mass;, sum_metal_mass
            
            R200_approx = radius[i]*scale
            M200_approx = total_mass[i]

            break
        endif
    endfor

    printf, lun1, FORMAT='(I,F,E,E,E,E)', number, redshift, R200, M200, R200_approx, M200_approx
            

    
endfor ; multiple checkpoints

close, lun1
free_lun, lun1
 
end
