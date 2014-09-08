pro make200, sstart, send, step

; File Options
;directory = "/ranger/scratch/01468/grapenut/current"
;file      = "/1star_hdf5_chk_"
directory = "/ranger/scratch/01707/mepa/Fry"
file      = "/final_hdf5_chk_"


;sstart    = 90
;send      = 90
;step      = 1
box_size  = 2.0*3.08D24

num_bins = 100

;comoving r_min/r_max
r_min = 1.0d21
r_max = 1.0d24
    
; distance within which to sample halo particle velocities, 100pc
halo_sample_distance = 3.0e24  
       
pi = 3.141592
h0 = 2.280e-18
G =  6.67e-8
omega_m = 0.275
density_thresh = (18.0 * pi^2.0) * 3.0 * h0^2.0 / 8.0 / pi / G * omega_m

print, density_thresh


; Output Options
plotps = 0
plotgif = 1
charsize = 1.4

;number = 50
;if (sstart gt send and step gt 0) then send = sstart
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


    
    ;;;;;;;;;;;;;;;;;;;;;
    ;;;; Grid data

    read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='gpot', TREE=tree, DATA=gpot, PARAMETERS=params
    
    read_amr_no_part, filename, VAR_NAME='pde', TREE=tree, DATA=pde, PARAMETERS=params

    read_amr_no_part, filename, VAR_NAME='velx', TREE=tree, DATA=velx, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='vely', TREE=tree, DATA=vely, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='velz', TREE=tree, DATA=velz, PARAMETERS=params
    
    max_lrefine = max(tree[*].lrefine)
    



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


star_x = max_x
star_y = max_y
star_z = max_z

print, "Star Center: ", star_x, star_y, star_z
print, "GPOT Center: ", max_x, max_y, max_z

dist_offset = sqrt((star_x - max_x)^2.0 + (star_y - max_y)^2.0 + (star_z - max_z)^2.0)

print, 'Distance   : ', dist_offset

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Average velocity offset of halo DM particles

       
       vx_offset = 0.0
       vy_offset = 0.0
       vz_offset = 0.0

print, "VELOCITY: ", vx_offset, vy_offset, vz_offset


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Sum quantities into bins


    baryon_mass = dblarr(num_bins+1) 
    pde_mass = dblarr(num_bins+1)
    volume = dblarr(num_bins+1)
    total_mass = dblarr(num_bins+1)
    metal_mass = dblarr(num_bins+1)
    total_velocity = dblarr(num_bins+1)
    metal_velocity = dblarr(num_bins+1)

    star_baryon_mass = dblarr(num_bins+1) 
    star_pde_mass = dblarr(num_bins+1)
    star_volume = dblarr(num_bins+1)
    star_total_mass = dblarr(num_bins+1)
    star_metal_mass = dblarr(num_bins+1)
    star_total_velocity = dblarr(num_bins+1)
    star_metal_velocity = dblarr(num_bins+1)

    radius = dblarr(num_bins+1)

    star_volume = 0.0
    sum_mass = 0.0
    sum_gas_mass = 0.0
    sum_metal_mass = 0.0

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
                        pde_mass_cell = pde[0,block_id,ii,jj,kk] * vol
                      
                        ;metal = c12[0,block_id,ii,jj,kk] + n14[0,block_id,ii,jj,kk] + o16[0,block_id,ii,jj,kk] + ne20[0,block_id,ii,jj,kk] + mg24[0,block_id,ii,jj,kk] + si28[0,block_id,ii,jj,kk] + s32[0,block_id,ii,jj,kk] + fe56[0,block_id,ii,jj,kk]
                        metal = 0.0
                        metal_mass_cell = metal * mass_cell
                        sum_mass = sum_mass + mass_cell + pde_mass_cell
                        sum_gas_mass = sum_gas_mass + mass_cell
                        sum_metal_mass = sum_metal_mass + metal_mass_cell
                        

                        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
                        ;;;; GPOT Centered

                        deltax = x_cell - max_x
                        deltay = y_cell - max_y
                        deltaz = z_cell - max_z
                        
    
                        distanceFromCenter = sqrt((deltax)^2.0  + $
                                                  (deltay)^2.0  + $ 
                                                  (deltaz)^2.0 )
                        

                        if(distanceFromCenter lt r_max) then begin
                            
                            vx = double(velx(0,block_id,ii,jj,kk)) - vx_offset
                            vy = double(vely(0,block_id,ii,jj,kk)) - vy_offset
                            vz = double(velz(0,block_id,ii,jj,kk)) - vz_offset
                            
                            if (distanceFromCenter lt r_min) then begin
                              bin = 0
                              radial_velocity = 0.0
                            endif else begin
                              bin = floor(double(num_bins) * alog(distanceFromCenter/r_min) / alog(r_max / r_min)) + 1
                              radial_velocity = (vx*deltax + vy*deltay + vz*deltaz) / distanceFromCenter
                            endelse
    
                            volume[bin] = volume[bin] + vol
     

                            total_velocity[bin] = total_velocity[bin] + mass_cell * radial_velocity
                            metal_velocity[bin] = metal_velocity[bin] + metal_mass_cell * radial_velocity
                            
                            baryon_mass[bin] = baryon_mass[bin] + mass_cell
                            pde_mass[bin] = pde_mass[bin] + pde_mass_cell
                            total_mass[bin] = total_mass[bin] + pde_mass_cell + mass_cell
                            metal_mass[bin] = metal_mass[bin] + metal_mass_cell
                            

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
          pde_density = 0.0
          metal_density = 0.0
          avg_velocity = 0.0
          avg_metal_velocity = 0.0
        endif else begin
          avg_density = total_mass[i]/volume[i]
          bar_density = baryon_mass[i]/volume[i]
          pde_density = pde_mass[i]/volume[i]
          metal_density = metal_mass[i]/volume[i]
          avg_velocity = total_velocity[i] / baryon_mass[i]
          avg_metal_velocity = metal_velocity[i] / metal_mass[i]
        endelse
        
        
        ;printf, lun, FORMAT='(E,E,E,E,E,E,E,E,E,E,E,E,E)', radius[i]*scale, volume[i]*scalecube, avg_density/scalecube, bar_density/scalecube, metal_density/scalecube, avg_velocity*scale, avg_metal_velocity*scale;, star_volume[i]*scalecube, star_avg_density/scalecube, star_bar_density/scalecube, star_metal_density/scalecube, star_avg_velocity*scale, star_avg_metal_velocity*scale

        printf, lun, FORMAT='(E,E,E,E)', radius[i]*scale, volume[i]*scalecube, avg_density/scalecube, bar_density/scalecube
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

    for i=0, num_bins do begin
        avg_density = total_mass[i]  / volume[i]
        
        if(avg_density lt density_thresh) then begin
            print, "R200: ", radius[i]*scale, ' cm   (', radius[i]*scale / 3.0856e18, ' parsec)'
            print, "M200: ", total_mass[i], ' g   (', total_mass[i] / 2e33, ' M_sun)'
            print, "Masses: ", sum_mass, sum_gas_mass, sum_metal_mass
            
            break
        endif
    endfor

endfor ; multiple checkpoints
 
end
