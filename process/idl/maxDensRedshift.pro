pro maxDensRedshift, sstart, send, step

; File Options
;directory = "/ranger/scratch/01707/mepa/Final"
;directory = "/ranger/scratch/01707/mepa/Fry"
;directory = "/ranger/scratch/01707/mepa/Rad"
;directory = "/ranger/scratch/01707/mepa/FryRad"
directory = "/ranger/scratch/01707/mepa/FryRad062"
;directory = "/data1/r900-3/mepa/SUBFIND/120725/Flash2Gadget_Shell_Script"
;directory = "/data1/r900-4/mepa/Radiation/Work/Fry"
file      = "/final_hdf5_chk_"

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

pfname = 'maxDens_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.txt'
openw,lun,pfname,/get_lun

outfile = 'maxdens_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'

num_pts = floor((send-sstart)/step)+1
maxdensities = dblarr(num_pts) 
redshifts = dblarr(num_pts) 

count = 0L
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

    max_dens_physical = max_value*(1.0 + redshift)^3.0
    maxdensities[count] = max_dens_physical
    
    print, "max density value and location"
    print, max_value, max_x, max_y, max_z
    print, "redshift
    print, redshift
    print, "physical max density"
    print, max_dens_physical

    printf, lun, FORMAT='(I,F,E)', number, redshift, max_dens_physical

    count = count + 1
endfor

close, lun
free_lun, lun

xr = [1e-4, 1e4]
yr = [1e0, 1e5]

;plot, redshifts, maxdensities, /xlog, /ylog, background='FFFFFF'xl, color=0, psym=3, xrange = xr, yrange=yr, xstyle=1, ystyle=1

plot, redshifts, maxdensities, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xtitle='redshift', ytitle='phys max dens (g cm^-3)'

plot, redshifts, maxdensities/1.67e-24, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xtitle='redshift', ytitle='phys max num dens (cm^-3)'
    
void = cgSnapshot(FILENAME=outfile, /NoDialog)

end
