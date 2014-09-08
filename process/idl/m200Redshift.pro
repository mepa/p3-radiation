pro m200Redshift, sstart, send, step

; File Options
;directory = "/scratch/01707/mepa/Rad_1.0sigma8"
;directory = "/scratch/01707/mepa/Rad_gridRes128_partRes512"
directory = "/scratch/01707/mepa/Rad_res512"
file      = "/rad_hdf5_chk_"

num_bins = 20000L

;comoving r_min/r_max
r_min = 5.0d21
r_max = 5.0d22

;num_bins = 100

;;comoving r_min/r_max
;r_min = 1.0d21 ;325 pc
;r_max = 1.0d24 ;325 kpc
    
pi = 3.1415927
h0 = 2.27502706e-18 ; 70.2 km/s/Mpc
G =  6.6726e-8
omega_m = 0.275
omega_v = 0.725

;outfile = 'halofrac_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'  
outfile = 'halomass_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'  

pfname1 = 'halo_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '_maxpden' + '.txt'
openw,lun1,pfname1,/get_lun

num_pts = floor((send-sstart)/step)+1

M200s = dblarr(num_pts)
bM200s = dblarr(num_pts)
dM200s = dblarr(num_pts)
bM200_over_dM200 = dblarr(num_pts)
bM200_over_dM200_approx = dblarr(num_pts)
redshifts = dblarr(num_pts) 

count = 0L
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
    redshifts[count] = redshift

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
    
    read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=pden, PARAMETERS=params


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;; Halo Center


    min_spacing = 1.0e55
    
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

    ;x_center = max_x_dens
    ;y_center = max_y_dens
    ;z_center = max_z_dens
    
    ;x_center = 3.39e24
    ;y_center = 3.11e24
    ;z_center = 3.12e24

    print, "min_gpot: ", min_gpot, min_x_gpot, min_y_gpot, min_z_gpot
    print, "max_pden: ", max_pden, max_x_pden, max_y_pden, max_z_pden
    print, "max_dens: ", max_dens, max_x_dens, max_y_dens, max_z_dens
    print, "x_center, y_center, z_center: ", x_center, y_center, z_center
 
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

                        deltax = x_cell - x_center
                        deltay = y_cell - y_center
                        deltaz = z_cell - z_center
                        
    
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



    ;pfname = 'profile_' + prefix + String(strcompress(number,/remove)) + '_' + String(strcompress(num_bins,/remove)) + '.txt'
    ;openw,lun,pfname,/get_lun

    scalecube = scale^3.0

    ;printf, lun, time, redshift

    for i=0L, num_bins do begin
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
        
        ;printf, lun, FORMAT='(E,E,E,E,E)', radius[i]*scale, volume[i]*scalecube, avg_density/scalecube, bar_density/scalecube, pden_density/scalecube
    endfor

    ;close, lun
    ;free_lun, lun


    
   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Find the virial mass/radius

    
    
    ;; Sum interior bins for average enclosed properties
    for i=1, num_bins do begin
       baryon_mass[i] = baryon_mass[i-1] + baryon_mass[i]
       pden_mass[i] = pden_mass[i-1] + pden_mass[i]
       total_mass[i] = total_mass[i-1] + total_mass[i]
       volume[i] = volume[i-1] + volume[i]
     endfor

    print, 'r_min: ', r_min*scale/3.0856e18, ' phy pc'
    print, 'r_max: ', r_max*scale/3.0856e18, ' phy pc'
        
    R200 = 0.0
    M200 = 0.0
    for i=0, num_bins do begin
        avg_density = total_mass[i]  / volume[i]
        
        if(avg_density lt density_thresh) then begin
            print, "R200: ", radius[i]*scale, ' cm   (', radius[i]*scale / 3.0856e18, ' phy pc)'
            print, "M200: ", total_mass[i], ' g   (', total_mass[i] / 2e33, ' M_sun)'
            print, "bM200: ", baryon_mass[i], ' g   (', baryon_mass[i] / 2e33, ' M_sun)'
            print, "dM200: ", pden_mass[i], ' g   (', pden_mass[i] / 2e33, ' M_sun)'
            print, "Masses: ", sum_mass, sum_gas_mass;, sum_metal_mass

            R200 = radius[i]*scale
            M200 = total_mass[i]
            bM200 = baryon_mass[i]
            dM200 = pden_mass[i]

            M200s[count] = M200
            bM200s[count] = bM200
            dM200s[count] = dM200

            bM200_over_dM200[count] = bM200 / dM200

            print, "bM200_over_dM200: ", bM200_over_dM200[count]

            break
         endif
    endfor

    if(R200 le r_min*scale) then begin
       print, 'WARNING: R200 = ',R200/scale, ' is le r_min = ', r_min, ' com cm'
    endif

    R200_approx = 0.0
    M200_approx = 0.0    
    for i=0, num_bins do begin
        avg_density = total_mass[i]  / volume[i]
        
        if(avg_density lt density_thresh_approx) then begin
            print, "R200_approx: ", radius[i]*scale, ' cm   (', radius[i]*scale / 3.0856e18, ' phy pc)'
            print, "M200_approx: ", total_mass[i], ' g   (', total_mass[i] / 2e33, ' M_sun)'
            print, "bM200_approx: ", baryon_mass[i], ' g   (', baryon_mass[i] / 2e33, ' M_sun)'
            print, "dM200_approx: ", pden_mass[i], ' g   (', pden_mass[i] / 2e33, ' M_sun)'
            print, "Masses: ", sum_mass, sum_gas_mass;, sum_metal_mass
            
            R200_approx = radius[i]*scale
            M200_approx = total_mass[i]
            bM200_approx = baryon_mass[i]
            dM200_approx = pden_mass[i]

            bM200_over_dM200_approx[count] = bM200_approx / dM200_approx

            print, "bM200_over_dM200_approx: ", bM200_over_dM200_approx[count]

            break
        endif
     endfor

    if(R200_approx le r_min*scale) then begin
       print, 'WARNING: R200_approx = ',R200_approx/scale, ' is le r_min = ', r_min, ' com cm'
    endif

    printf, lun1, FORMAT='(I,F,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E,E)', number, redshift, R200, M200, bM200, dM200, R200_approx, M200_approx, bM200_approx, dM200_approx, min_gpot, min_x_gpot, min_y_gpot, min_z_gpot, max_pden, max_x_pden, max_y_pden, max_z_pden, max_dens, max_x_dens, max_y_dens, max_z_dens

    count = count + 1
     
endfor ; multiple checkpoints

close, lun1
free_lun, lun1

plot, redshifts, M200s, /ylog, background='FFFFFF'xl, color=0, psym=-0, linestyle=0, xtitle='redshift', ytitle='mass within R200'

oplot, redshifts, bM200s, linestyle=1

oplot, redshifts, dM200s, linestyle=2
;xyouts, 

;plot, redshifts, bM200_over_dM200, /ylog, background='FFFFFF'xl, color=0, psym=-0, linestyle=0, xtitle='redshift', ytitle='baryonic to dark matter fraction within R200'

;oplot, redshifts, bM200_over_dM200_approx, color=1, linestyle=1


void = cgSnapshot(FILENAME=outfile, /NoDialog)
 
end
