pro mVirRedshift, sstart, send, step

;directory = "/work/01707/mepa/Rad_1.0sigma8/newChem"
directory = "/scratch/01707/mepa/Rad_1.0sigma8"
file      = "/rad_hdf5_chk_"

mu = 0.6

pi = 3.1415927e0
grav_const = 6.6726e-8
k = 1.3806503e-16
mass_H = 1.67e-24

pc = 3.0856775807e18
kpc = 3.0856775807e21
Mpc = 3.0856775807e24
Msun = 1.9889225e33

box_size = 2.0e0 * Mpc

hubble_const = 2.27502706e-18 ; 70.2 km/s/Mpc
omega_m = 0.275e0
omega_v = 0.725e0

pfname1 = 'halo_dm_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '_maxpden' + '.txt'
openw,lun1,pfname1,/get_lun

for number = sstart,send,step do begin

   if (number ge 1)   then prefix = '000'
   if (number ge 10)   then prefix = '00'
   if (number ge 100)  then prefix = '0'
   if (number ge 1000) then prefix = ''

   filename = directory + file + prefix + String(strcompress(number,/remove))
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Read input data
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;; Redshift
   ;;;; scale_factor (redshift)
   
   file_identifier = H5F_OPEN(filename)
   dataset = H5D_OPEN(file_identifier, "real scalars")
   real_scalars = H5D_READ(dataset)
   for i=1, (size(real_scalars))[3] do begin
      if (stregex(real_scalars[i-1].name, '^scalefactor', /BOOLEAN)) then scale_factor = real_scalars[i-1].value
      if (stregex(real_scalars[i-1].name, '^time', /BOOLEAN)) then time = real_scalars[i-1].value
   endfor
   
   redshift = 1.0 / scale_factor - 1.0
   one_plus_redshift = (1.0 + redshift)

   print, "********************"
   print, "redshift:", redshift
   
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;; Dark Matter Particles
   ;;;; n_particles, posx, posy, posz, mass, tag
   
   file_identifier = H5F_OPEN(filename)
   file_contents = h5_parse(filename)
   
   tag_contents = tag_names(file_contents)
   tag_where = where(tag_contents EQ "PARTICLE_NAMES")
   
   dataset = H5D_OPEN(file_identifier, "tracer particles")
   particles_read = H5D_READ(dataset)
   H5D_CLOSE, dataset

   numDims = (size(particles_read))[0]
   if (numDims EQ 1) then begin $
      particles_read = reform(particles_read, (size(particles_read, /DIMENSIONS))[0], 1)
   endif
   ;print, "********************"
   ;print, "size(particles_read):", size(particles_read)
   
   n_particles = (size(particles_read, /DIMENSIONS))[1]
   
   print, "n_particles:", n_particles
   
   min_mass = 1.0e6 * Msun
   max_mass = 0.0
   x = dblarr(n_particles)
   y = dblarr(n_particles)
   z = dblarr(n_particles)
   M_particle = dblarr(n_particles)
   for i = 0L, n_particles-1L do begin
      
      posx = particles_read(9,i) 
      posy = particles_read(10,i)
      posz = particles_read(11,i)
      mass = particles_read(5,i) 
      ;tag  = particles_read(13,i)
 
      x[i] = posx
      y[i] = posy
      z[i] = posz
      M_particle[i] = mass
      if (mass lt min_mass) then min_mass = mass
      if (mass gt max_mass) then max_mass = mass

   endfor 
   print, "min_mass, max_mass (Msun):", min_mass / Msun, max_mass / Msun


   ;;;;;;;;;;;;;;;;;;;;;;;;;;; Grid Data
   ;;;; gpot, pden, dens

   read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
   read_amr_no_part, filename, VAR_NAME='gpot', TREE=tree, DATA=gpot, PARAMETERS=params
   read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=pden, PARAMETERS=params


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Find Halo Center

   min_gpot = 0.0
   max_pden = 0.0
   max_dens = 0.0
   
   min_x_gpot = 0.0
   min_y_gpot = 0.0
   min_z_gpot = 0.0

   max_x_pden = 0.0
   max_y_pden = 0.0
   max_z_pden = 0.0

   max_x_dens = 0.0
   max_y_dens = 0.0
   max_z_dens = 0.0
   
   for block_id = 0L, params.totBlocks - 1L, 1 do begin
      
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
                          
                  ; find min gpot
                  if gpot(0,block_id,ii,jj,kk) lt min_gpot then begin
                            
                     min_gpot = gpot(0,block_id,ii,jj,kk)
                     min_x_gpot = x_cell
                     min_y_gpot = y_cell
                     min_z_gpot = z_cell
                     
                  endif

                  ; find max pden
                  if pden(0,block_id,ii,jj,kk) gt max_pden then begin
                     
                     max_pden = pden(0,block_id,ii,jj,kk)
                     max_x_pden = x_cell
                     max_y_pden = y_cell
                     max_z_pden = z_cell
                     
                  endif

                  ; find max dens
                  if dens(0,block_id,ii,jj,kk) gt max_dens then begin
                     
                     max_dens = dens(0,block_id,ii,jj,kk)
                     max_x_dens = x_cell
                     max_y_dens = y_cell
                     max_z_dens = z_cell
                     
                  endif
                  
               endfor
            endfor    
         endfor
      endif      
   endfor

   ;x_center = min_x_gpot
   ;y_center = min_y_gpot
   ;z_center = min_z_gpot

   x_center = max_x_pden
   y_center = max_y_pden
   z_center = max_z_pden

   ;x_center = 3.39e24
   ;y_center = 3.11e24
   ;z_center = 3.12e24

   print, "min_gpot: ", min_gpot, min_x_gpot, min_y_gpot, min_z_gpot
   print, "max_pden: ", max_pden, max_x_pden, max_y_pden, max_z_pden
   print, "max_dens: ", max_dens, max_x_dens, max_y_dens, max_z_dens
   print, "x_center, y_center, z_center: ", x_center, y_center, z_center


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Calculate Threshold Density

   hubble_param_sq = hubble_const^2.0 * (omega_m / scale_factor^3.0 + omega_v)

   rho_crit = 3.0 * hubble_param_sq / (8.0 * pi * grav_const) / one_plus_redshift^3.0 ; exact
   rho_background = rho_crit * omega_m
   rho_crit_200 = 200.0 * rho_crit
   rho_background_200 = 200.0 * rho_background
   
   rho_crit_approx = 3.0 * hubble_const^2.0 / (8.0 * pi * grav_const) ; approximate 
   rho_background_approx = rho_crit_approx * omega_m
   rho_crit_200_approx = 200.0 * rho_crit_approx
   rho_background_200_approx = 200.0 * rho_background_approx

   print, "rho_crit_200: ", rho_crit_200
   print, "rho_background_200: ", rho_background_200
   print, "rho_crit_200_approx: ", rho_crit_200
   print, "rho_background_200_approx: ", rho_crit_200

 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Sort Particles By Radial Distance from Halo Center

   ;print, "********************"
   ;print, "Begin sorting..."
   r = dblarr(n_particles)
   r_sorted = dblarr(n_particles)
   for i = 0L, n_particles - 1L do begin
      r[i] = sqrt( (x[i] - x_center)^2 + (y[i] - y_center)^2 + (z[i] - z_center)^2 )
   endfor
   r_sorted_indices = sort(r)
   r_sorted = r[r_sorted_indices]
   ;print, r_sorted
   ;print, "...finished sorting"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Sum Up Enclosed Mass

   ;print, "********************"
   ;print, "Begin calculating enclosed mass..."
   M_enclosed = dblarr(n_particles)
   ;M_enclosed[0] = 0.0                             ; don't include particle at r[1]
   M_enclosed[0] = M_enclosed[r_sorted_indices[0]] ; include particle at r[1]
   for i = 1L, n_particles - 1L do begin     
      M_enclosed[i] = M_enclosed[i - 1] + M_particle[r_sorted_indices[i]]
   endfor
   ;print, "...finished calculating enclosed mass"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Find Virial Radius and Mass

   ;print, "********************"
   ;print, "Begin calculating virial radius..."
   r_200 = 0.0
   M_200 = 0.0
   j = 0
   rho = dblarr(n_particles)
   for i = 0L, n_particles - 1L do begin
      rho[i] = M_enclosed[i] / (4.0 * pi * r_sorted[i]^3 / 3.0)
      if(rho[i] ge rho_crit_200) then begin
         r_200 = r_sorted[i]
         M_200 = M_enclosed[i]
         j = j + 1
      endif
      if(rho[i] ge rho_crit_200_approx) then begin
         r_200_approx = r_sorted[i]
         M_200_approx = M_enclosed[i]
      endif
      if(rho[i] ge rho_background_200) then begin
         r_background_200 = r_sorted[i]
         M_background_200 = M_enclosed[i]
      endif
      if(rho[i] ge rho_background_200_approx) then begin
         r_background_200_approx = r_sorted[i]
         M_background_200_approx = M_enclosed[i]
      endif
;   if(M_enclosed[i] le 1.0e9*Msun) then begin
;      r_9 = r_sorted[i]
;      M_9 = M_enclosed[i]
;      rho_9 = rho[i]
;   endif
;   if(r_sorted[i] le 3.0*kpc/scale_factor) then begin
;      r_3kpc = r_sorted[i]
;      M_3kpc = M_enclosed[i]
;      rho_3kpc = rho[i]
;   endif

   endfor
   ;print, "...finished calculating virial radius"

   print, "********************"
   print, "Num Particles Within r_200:", j
   print, "********************"
   print, "r_200 (physical pc):", r_200 * scale_factor / pc
   print, "r_200 (comoving pc):", r_200 / pc
   print, "M_200 (Msun):", M_200 / Msun
   print, "********************"
   print, "r_200_approx (physical pc):", r_200_approx * scale_factor / pc
   print, "r_200_approx (comoving pc):", r_200_approx / pc
   print, "M_200_approx (Msun):", M_200_approx / Msun
   print, "********************"
   print, "r_background_200 (physical pc):", r_background_200 * scale_factor / pc
   print, "r_background_200 (comoving pc):", r_background_200 / pc
   print, "M_background_200 (Msun):", M_background_200 / Msun
   print, "********************"
   print, "r_background_200_approx (physical pc):", r_background_200_approx * scale_factor / pc
   print, "r_background_200_approx (comoving pc):", r_background_200_approx / pc
   print, "M_background_200_approx (Msun):", M_background_200_approx / Msun
   print, "********************"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Calculate Virial Temperature and Circular
;;;;;;;;;;;;;;;;;;;;;;;; Velocity

   v_circ = sqrt(grav_const * M_200 / r_200)
   ;v_circ = sqrt(grav_const * 1.1 * M_200 / r_200)
   ;v_vir = v_circ

   ;T_vir = mu * mass_H * v_vir^2 / 3 / k

   Mvir = M_200 / Msun
   ;Mvir = 1.1 * M_200 / Msun
   Tvir = 2.0e4 * (mu / 1.2) * (Mvir / 1.0e8)^(2/3) * (1.0 + redshift) / 10.0

   print, "v_circ (km/s):", v_circ / 1.0e4
   ;print, "T_vir (K):", T_vir
   print, "Tvir (K):", Tvir
   print, "********************"


                                ;write(11,format), redshift, M_200,
                                ;r_200, j, x_center, y_center,
                                ;z_center, mass, mcrit200,
                                ;M_background_200_approx

   printf, lun1, FORMAT='(I,F,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E)', number, redshift, r_200, M_200, r_background_200_approx, M_background_200_approx, min_gpot, min_x_gpot, min_y_gpot, min_z_gpot, max_pden, max_x_pden, max_y_pden, max_z_pden, max_dens, max_x_dens, max_y_dens, max_z_dens, v_circ, Tvir

endfor

close, lun1
free_lun, lun1

end
