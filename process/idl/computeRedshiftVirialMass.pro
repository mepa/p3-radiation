pro computeRedshiftVirialMass, sstart, send, step

; File Options
;directory = "/data1/r900-1/ctss/FLASH_3.2/TestCosmo/"
directory = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128"
file      = "/radCosmoLW_hdf5_chk_"

;center of halo/max(DENS) and z=12.87
;x_center = 1.5275718995488837E+24
;y_center = 1.4393370395294772E+24
;z_center = ?

;bhmass    = 1.35D34

; Output Options
plotps = 0
plotgif = 1
charsize = 1.4

numOfFiles = send - sstart + 1

fileRedshift = dblarr(numOfFiles)
mass     = dblarr(numOfFiles)
Msun = 1.989e33

print, "numOfFiles = ", numOfFiles

pfname = 'plots/minGpot_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
openw,lun,pfname,/get_lun

for number = sstart, send, step do begin

   index = number - sstart

   if (number ge 1)   then prefix = '000'
   if (number ge 10)   then prefix = '00'
   if (number ge 100)  then prefix = '0'
   if (number ge 1000) then prefix = ''
   
   filename = directory + file + prefix + String(strcompress(number,/remove))
   print, filename
   
   if (number eq 228) then continue ;TACC lost this file
   ;;;;;
   ; get redshift
   file_identifier = H5F_OPEN(filename)
   dataset = H5D_OPEN(file_identifier, "real scalars")
   real_scalars = H5D_READ(dataset)
   for i=1, (size(real_scalars))[3] do begin
      if (stregex(real_scalars[i-1].name, '^scale', /BOOLEAN)) then scale = real_scalars[i-1].value
   endfor
   ;;;;;
    
   redshift = (1.0 / scale) - 1.0
   oneplusred = 1.0 + redshift

   print, "index=", index
   print, "redshift=", redshift
   fileRedshift[index] = redshift
    
   ; get tree and data
   read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
   read_amr_no_part, filename, VAR_NAME='gpot', TREE=tree, DATA=gpot, PARAMETERS=params
    
   read_amr_no_part, filename, VAR_NAME='pde', TREE=tree, DATA=pde, PARAMETERS=params
   read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=pden, PARAMETERS=params
    
    max_lrefine = max(tree[*].lrefine)
    
    ; get particle info - 
    
   ; read_write_particle_positions, filename, POSITIONS=positions, P_MASS=p_mass, $
    ;  NUMBER_PARTICLES=number_particles

    print, "max lrefine = ", max_lrefine

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; Now the meat of the program
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    
    ;;;;;;; Where to we want these radial bins to be centered?
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


          ;using_gpot = .FALSE.            
          ; loop over cells
          for ii=0, params.nxb-1 do begin
             for jj=0, params.nyb-1 do begin
                for kk=0, params.nzb-1 do begin
                   
                   x_cell = block_x_min + (ii + 0.5) * dx
                   y_cell = block_y_min + (jj + 0.5) * dy
                   z_cell = block_z_min + (kk + 0.5) * dz
                                           
                   ; let's center it on minimum of grav potential!
                   if gpot(0,block_id,ii,jj,kk) lt min_value then begin
                      ;using_gpot = .TRUE.
                      
                      min_value = gpot(0,block_id,ii,jj,kk)
                      max_x = x_cell
                      max_y = y_cell
                      max_z = z_cell
                      
                   endif
                        
                   ; alternatively, max density
                   ;if dens(0,block_id,ii,jj,kk) gt max_value then begin
                   ;    
                   ;    max_value = dens(0,block_id,ii,jj,kk)
                   ;    max_x = x_cell
                   ;    max_y = y_cell
                   ;    max_z = z_cell
                   ;    
                   ;endif
                            
                endfor
             endfor    
          endfor
       endif      
    endfor
    
    x_center = max_x
    y_center = max_y
    z_center = max_z

    printf, lun, FORMAT='(I,F,E,E,E,E)', number, redshift, min_value, max_x, max_y, max_z

    ;;;;;;;;;;;;;;;

    print, "Centering radial profile at"
    print, "x=", max_x
    print, "y=", max_y
    print, "z=", max_z
    print, "smallest cell, comoving=", min_spacing

    ;;;;;
    num_bins = 25
    baryon_density = dblarr(num_bins)
    pde_density = dblarr(num_bins)
    pden_density = dblarr(num_bins)
    DM_density = dblarr(num_bins)
    radial_coordinate = dblarr(num_bins)
    
    ;comoving r_min/r_max
    r_min = min_spacing / 2.0 
    r_max = 2.0e23
  
    total_mass_enclosed_baryon = 0.0
    total_mass_enclosed_pde = 0.0
    total_mass_enclosed_pden = 0.0
    
    total_mass_enclosed = 0.0
    
    
    omega_m = 0.25
    ;density_thresh = (18.0 * 3.141592^2)*3.0 * (2.280e-18)^2 / 8.0 / 3.141592 / 6.67e-8
    density_thresh =  (18.0 * 3.141592^2) * 3.0 * (2.280e-18)^2 / 8.0 / 3.141592 / 6.67e-8
       
    print, "density threshold = ", density_thresh
     
    for i=0, num_bins - 1 do begin

       mass_enclosed_baryon = 0.0
       mass_enclosed_pde = 0.0
       mass_enclosed_pden = 0.0
               
       print, "stepping through radius in baryons...", i
           
       min_r = r_min * 10^((i-1) * alog10(r_max / r_min) / (num_bins -1))
       max_r =  r_min * 10^(i * alog10(r_max / r_min) / (num_bins-1))
       if (i eq 0) then begin
          min_r = 0.0
          max_r = r_min
       endif
       
       print, "limits between...", min_r, max_r
       
       ;print, i, min_r, max_r
        
       vol_shell = (4.0*3.141592/3.0) * ((max_r^3.0) - (min_r^3.0))
       
       total_volume =  (4.0*3.141592/3.0) * (max_r^3.0)
       
    
       ; Radial profile of baryons starting from x_center, y_center, z_center  - ----
    
       ; loop over blocks

       ; trying to find point where average total density of halo
       ; is equal to 18 * pi^2 * rho_tot_background
    
       for block_id = 0L, params.totBlocks - 1L, 1 do begin
          
          if (tree[block_id].nodetype EQ 1) then begin
             
             dx = tree[block_id].size[0] / params.nxb ;/ oneplusred
             dy = tree[block_id].size[1] / params.nyb ; / oneplusred
             dz = tree[block_id].size[2] / params.nzb ; / oneplusred
             
             if (dx lt min_spacing) then begin
                min_spacing = dx
             endif
                         
             block_x_min = tree[block_id].bndBox[0,0] ; / oneplusred
             block_y_min = tree[block_id].bndBox[0,1] ; / oneplusred
             block_z_min = tree[block_id].bndBox[0,2] ; / oneplusred
             
             block_x_max = tree[block_id].bndBox[1,0] ; / oneplusred
             block_y_max = tree[block_id].bndBox[1,1] ; / oneplusred
             block_z_max = tree[block_id].bndBox[1,2] ; / oneplusred
                         
             ; loop over cells
             for ii=0, params.nxb-1 do begin
                for jj=0, params.nyb-1 do begin
                   for kk=0, params.nzb-1 do begin
                      
                      x_cell = block_x_min + (ii + 0.5) * dx
                      y_cell = block_y_min + (jj + 0.5) * dy
                      z_cell = block_z_min + (kk + 0.5) * dz
                      
                      deltax = x_cell - x_center
                      deltay = y_cell - y_center
                      deltaz = z_cell - z_center
                          
                      distanceFromCenter = sqrt((deltax)^2.0  + $
                                                (deltay)^2.0  + $ 
                                                (deltaz)^2.0 )

                      ; Check to see if particular cell lies within shell
                      if (distanceFromCenter LT max_r) then begin
                            
                         mass_cell = dens[0,block_id,ii,jj,kk] * dx * dy * dz
                         pde_mass_cell = pde[0,block_id,ii,jj,kk] *dx*dy*dz
                         pden_mass_cell =  pden[0,block_id,ii,jj,kk] *dx*dy*dz
                            
                         ;print, mass_cell
                            
                         mass_enclosed_baryon = mass_enclosed_baryon + mass_cell
                         mass_enclosed_pde = mass_enclosed_pde +pde_mass_cell
                         mass_enclosed_pden = mass_enclosed_pden +pden_mass_cell
                            
                         total_mass_enclosed_baryon = total_mass_enclosed_baryon + mass_cell
                         total_mass_enclosed_pde = total_mass_enclosed_pde +pde_mass_cell
                         total_mass_enclosed_pden = total_mass_enclosed_pden +pden_mass_cell
                            
                         ;total_mass_enclosed =total_mass_enclosed +  mass_cell + pden_mass_cell
                         total_mass_enclosed =total_mass_enclosed +  mass_cell + pde_mass_cell

                      endif
                      
                   endfor
                endfor    
             endfor
          endif      
       endfor
       
       avg_tot_density = total_mass_enclosed / total_volume
       
       if (avg_tot_density le density_thresh) then begin
        
          print, "virial radius=", (min_r + max_r)/2.0
          print, "virial mass (total) = ", total_mass_enclosed
        
          mass[index] = total_mass_enclosed
        
          BREAK
        
       endif
         
    endfor
    
 endfor
close, lun
free_lun, lun

print, fileRedshift
print, "****"
print, mass

openw, 33, "massHist.dat"
for i=0, numOfFiles-1, step do begin
   printf, 33, fileRedshift[i], mass[i]
endfor
close, 33

plot, fileRedshift, mass/Msun, /ylog, color=0, background='FFFFFF'xl, xlabel='redshift', ylabel='virial mass (M_sun)'

print, "total mass enclosed in baryons = ", total_mass_enclosed_baryon
print, "total mass enclosed pde = ", total_mass_enclosed_pde
print, "total mass enclosed pden = ", total_mass_enclosed_pden
min_spacing = min_spacing / oneplusred
print, "smallest cell size, physical=", min_spacing
print, "maximum refinement level on grid", max_lrefine
print, "redshift=", redshift

end
