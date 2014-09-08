pro calcTotalMass, sstart, send, step

; File Options
;directory = "/data1/r900-4/mepa/Radiation/Work/Fry"
;directory = "/work/01707/mepa/Rad"
;directory = "/scratch/01707/mepa/zInitial_74.2/eintSwitch_1.0"
;directory = "/work/01707/mepa/Rad_1.0sigma8/zInitial_149.4"
directory = "/work/01707/mepa/Rad_1.0sigma8/newChem"
file      = "/rad_hdf5_chk_"
    
pi = 3.1415927
h0 = 2.27502706e-18 ; 70.2 km/s/Mpc
G =  6.6726e-8
omega_m = 0.275
omega_v = 0.725

;pfname1 = 'halo_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
;openw,lun1,pfname1,/get_lun

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

    scalecube = scale^3.0

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
    
    read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=pden, PARAMETERS=params


    total_mass = 0.0
    total_gas_mass = 0.0
    total_dm_mass = 0.0
    total_volume = 0.0

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

                        vol_cell = dx*dy*dz
                        dens_cell = dens[0,block_id,ii,jj,kk]
                        pden_cell = pden[0,block_id,ii,jj,kk]
                        gas_mass_cell = dens_cell * vol_cell
                        dm_mass_cell = pden_cell * vol_cell
                        
                        total_volume = total_volume + vol_cell
                        total_gas_mass = total_gas_mass + gas_mass_cell
                        total_dm_mass = total_dm_mass + dm_mass_cell
                        total_mass = total_mass + gas_mass_cell + dm_mass_cell
      
                     endfor
                 endfor    
             endfor
         endif      
     endfor
    
    ave_baryon_dens = total_gas_mass / total_volume
    ave_darkmatter_dens = total_dm_mass / total_volume
    ave_total_dens = total_mass / total_volume

    print, "(physical g/cm^3)"
    print, "ave_baryon_dens :", ave_baryon_dens / scalecube
    print, "ave_darkmatter_dens:", ave_darkmatter_dens / scalecube
    print, "ave_total_dens:", ave_total_dens / scalecube
                                ;print, "ave_baryon_dens /
                                ;ave_darkmatter_dens:" ave_baryon_dens
                                ;/ ave_darkmatter_dens
    print, "total_dm_mass / total_gas_mass:", total_dm_mass / total_gas_mass
  



    
endfor ; multiple checkpoints

;close, lun1
;free_lun, lun1
 
end
