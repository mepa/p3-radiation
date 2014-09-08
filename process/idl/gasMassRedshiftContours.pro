pro gasMassRedshiftContours, sstart, send, step

; File Options
directory = "/scratch/01707/mepa/Rad_1.0sigma8"
;directory = "/data1/r900-3/mepa/SUBFIND/120725/Flash2Gadget_Shell_Script"
;directory = "/data1/r900-4/mepa/Radiation/Work/Fry"
;directory = "/data1/r900-3/mepa/Work/FryRad062"
;file      = "/final_hdf5_chk_"

file      = "/rad_hdf5_chk_"

linecolors

mass_sun = 2.0D33

count = 0L
num_pts = floor((send-sstart)/step)+1
print, "num_pts =", num_pts
ionized_gas = dblarr(num_pts)
neutral_gas = dblarr(num_pts)
molecular_gas = dblarr(num_pts)
h_gas = dblarr(num_pts)
hplu_gas = dblarr(num_pts)
hel_gas = dblarr(num_pts)
hep_gas = dblarr(num_pts)
hepp_gas = dblarr(num_pts)
hmin_gas = dblarr(num_pts)
htwo_gas = dblarr(num_pts)
htwp_gas = dblarr(num_pts)
hd_gas = dblarr(num_pts)
deut_gas = dblarr(num_pts)
dplu_gas = dblarr(num_pts)
elec_gas = dblarr(num_pts)
h_gas_core = dblarr(num_pts)
gas_10neg2 = dblarr(num_pts)
gas_10neg1 = dblarr(num_pts)
gas_10zero = dblarr(num_pts)
gas_10pos1 = dblarr(num_pts)
gas_10pos2 = dblarr(num_pts)
gas_10pos3 = dblarr(num_pts)
gas_10pos4 = dblarr(num_pts)

redshifts = dblarr(num_pts) 

pfname = 'gasContours_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
openw,lun,pfname,/get_lun

outfile = 'gasContours_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'

for number = sstart,send,step do begin

    num_cells = 0
    
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

    redshifts[count] = redshift

    ;;;;;
    ; get tree
    ;print, "start getting tree..."
    read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='h', TREE=tree, DATA=h, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='hplu', TREE=tree, DATA=hplu, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='hel', TREE=tree, DATA=hel, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='hep', TREE=tree, DATA=hep, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='hepp', TREE=tree, DATA=hepp, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='hmin', TREE=tree, DATA=hmin, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='htwo', TREE=tree, DATA=htwo, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='htwp', TREE=tree, DATA=htwp, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='deut', TREE=tree, DATA=deut, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='dplu', TREE=tree, DATA=dplu, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='hd', TREE=tree, DATA=hd, PARAMETERS=params
    ;read_amr_no_part, filename, VAR_NAME='elec', TREE=tree, DATA=elec, PARAMETERS=params
    ;print, "...done"

    mass_10neg2 = 0.0D
    mass_10neg1 = 0.0D
    mass_10zero = 0.0D
    mass_10pos1 = 0.0D
    mass_10pos2 = 0.0D
    mass_10pos3 = 0.0D
    mass_10pos4 = 0.0D
    ;print, "start looping over blocks..."
    ; loop over blocks
    for block_id = 0, params.totBlocks - 1, 1 do begin
        
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
                        
                        ;x_cell = block_x_min + (ii + 0.5) * dx
                        ;y_cell = block_y_min + (jj + 0.5) * dy
                        ;z_cell = block_z_min + (kk + 0.5) * dz
                        
                        mass_cell = dens[0,block_id,ii,jj,kk] * dx * dy * dz
                        
                        mass_h_cell = h[0,block_id,ii,jj,kk] * mass_cell
                        mass_hplu_cell = hplu[0,block_id,ii,jj,kk] * mass_cell
                        mass_htwo_cell = htwo[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_hel_cell = hel[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_hep_cell = hep[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_hepp_cell = hepp[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_hmin_cell = hmin[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_htwp_cell = htwo[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_deut_cell = deut[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_dplu_cell = dplu[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_hd_cell = hd[0,block_id,ii,jj,kk] * mass_cell
                        ;mass_elec_cell = elec[0,block_id,ii,jj,kk] * mass_cell

                        ;mass_total_cell = mass_h_cell + mass_hplu_cell + mass_htwo_cell + mass_hel_cell + mass_hep_cell + mass_hepp_cell + mass_hmin_cell + mass_htwp_cell + mass_deut_cell + mass_dplu_cell + mass_hd_cell + mass_elec_cell

                        mass_total_cell_approx = mass_h_cell + mass_hplu_cell + mass_htwo_cell

                        ;total_mass_frac = h[0,block_id,ii,jj,kk] + $
                        ;  hplu[0,block_id,ii,jj,kk] + $
                        ;  hel[0,block_id,ii,jj,kk] + $
                        ;  hep[0,block_id,ii,jj,kk] + $
                        ;  hepp[0,block_id,ii,jj,kk] + $
                        ;  hmin[0,block_id,ii,jj,kk] + $
                        ;  htwo[0,block_id,ii,jj,kk] + $
                        ;  htwp[0,block_id,ii,jj,kk] + $
                        ;  deut[0,block_id,ii,jj,kk] + $
                        ;  dplu[0,block_id,ii,jj,kk] + $
                        ;  hd[0,block_id,ii,jj,kk] + $
                        ;  elec[0,block_id,ii,jj,kk]
                        ;if (total_mass_frac le 1.0D-5) then begin
                        ;    print, "WARNING: total mass fraction = ", total_mass_frac 
                        ;endif
                        

                        pmass = 1.67e-24
                        numdens = dens[0,block_id,ii,jj,kk] * oneplusred^3.0 / pmass
                        if (numdens GT 1.0e-2) then begin
                            mass_10neg2 = mass_10neg2 + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e-1) then begin
                            mass_10neg1 = mass_10neg1 + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e0) then begin
                            mass_10zero = mass_10zero + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e1) then begin
                            mass_10pos1 = mass_10pos1 + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e2) then begin
                            mass_10pos2 = mass_10pos2 + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e3) then begin
                            mass_10pos3 = mass_10pos3 + mass_total_cell_approx
                        endif
                        if (numdens GT 1.0e4) then begin
                            mass_10pos4 = mass_10pos4 + mass_total_cell_approx
                        endif
   
                    endfor
                endfor    
            endfor
        endif       ; leaf blocks
    endfor     ;blocks

    gas_10neg2[count] = mass_10neg1 / mass_sun
    gas_10neg1[count] = mass_10neg2 / mass_sun
    gas_10zero[count] = mass_10zero / mass_sun
    gas_10pos1[count] = mass_10pos1 / mass_sun
    gas_10pos2[count] = mass_10pos2 / mass_sun
    gas_10pos3[count] = mass_10pos3 / mass_sun
    gas_10pos4[count] = mass_10pos4 / mass_sun
        
    
    printf, lun, FORMAT='(I,F,E,E,E,E,E,E,E)', number, redshift, gas_10neg2[count], gas_10neg1[count], gas_10zero[count], gas_10pos1[count], gas_10pos2[count], gas_10pos3[count], gas_10pos4[count]
    
    ;print, h_gas[count], hel_gas[count], hplu_gas[count], hep_gas[count], hepp_gas[count], htwo_gas[count], hmin_gas[count], htwp_gas[count], hd_gas[count], deut_gas[count], dplu_gas[count], elec_gas[count]
    
    count = count + 1

endfor

close, lun
free_lun, lun



xr = [16.6,17.1]
yr = [1.0D0, 1.0D10]

plot, redshifts, gas_10neg2, psym=-3, /ylog, background='FFFFFF'xl, color=1, linestyle=0, yrange=yr, ystyle=1, xtitle='redshift', ytitle='total gas mass in a 2 Mpc comoving box (Msun)'
;plot, redshifts, h_gas, psym=-3, /ylog, background='FFFFFF'xl, color=0, linestyle=0, xtitle='redshift', ytitle='total gas mass in H, H+, H2, H(>10^{-4}) [Msun]'
oplot, redshifts, gas_10neg2, psym=-3, color=1, linestyle=1
oplot, redshifts, gas_10neg1, psym=-3, color=2, linestyle=1
oplot, redshifts, gas_10zero, psym=-3, color=3, linestyle=2
oplot, redshifts, gas_10pos1, psym=-3, color=4, linestyle=3
oplot, redshifts, gas_10pos2, psym=-3, color=5, linestyle=4
oplot, redshifts, gas_10pos3, psym=-3, color=6, linestyle=5
oplot, redshifts, gas_10pos4, psym=-3, color=7, linestyle=0

print, gas_10neg2[0], gas_10neg1[0], gas_10zero[0], gas_10pos1[0], gas_10pos2[0], gas_10pos3[0], gas_10pos4[0]
print, gas_10neg2[count-1], gas_10neg1[count-1], gas_10zero[count-1], gas_10pos1[count-1], gas_10pos2[count-1], gas_10pos3[count-1], gas_10pos4[count-1]


    
;plot, ionizedmass, redshift, /ylog, background='FFFFFF'xl, Color=0, psym=-4, yrange=yr, ystyle=1
;plot, ionizedmass, redshift, /ylog, background='FFFFFF'xl, Color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1
;plot, redshifts, ionized_gas, /ylog, background='FFFFFF'xl, color=0, psym=0
;plot, redshifts, hplu_gas, psym=-4, /ylog, background='FFFFFF'xl, color=0, linestyle=0
;oplot, redshifts, hel_gas, linstyle=1
;oplot, redshifts, hep_gas, linestyle=3
;oplot, redshifts, hepp_gas, linestyle=4
;oplot, redshifts, hmin_gas, linestyle=0
;oplot, redshifts, htwp_gas
;oplot, redshifts, hd_gas
;oplot, redshifts, deut_gas
;oplot, redshifts, dplu_gas
;Oplot, elec_gas
    
void = cgSnapshot(FILENAME=outfile, /NoDialog)
end
