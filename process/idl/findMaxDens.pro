pro findMaxDens

; File Options
directory = "/data1/r900-4/mepa/Radiation/Work/cosmoDMonly_256/"
file      = "cosmoDMonly_256_hdf5_chk_"
;num       = 0018
;sstart    = 206 ;?
;send      = 235 ;?
;step      = 1
;box_size  = 13.5769814D24 ;3.08D24

;bhmass    = 1.35D34

; Output Options
;plotps = 0
;plotgif = 1
;charsize = 1.4

number = 89 

;if (sstart gt send and step gt 0) then send = ssta27
;for number = sstart,send,step do begin

if (number ge 0)   then prefix = '000'
if (number ge 10)   then prefix = '00'
if (number ge 100)  then prefix = '0'
if (number ge 1000) then prefix = ''

filename = directory + file + prefix + String(strcompress(number,/remove))

;;;;;
; get redshift of file
file_identifier = H5F_OPEN(filename)
dataset = H5D_OPEN(file_identifier, "real scalars")
real_scalars = H5D_READ(dataset)
for i=1, (size(real_scalars))[3] do begin
    if (stregex(real_scalars[i-1].name, '^redshift', /BOOLEAN)) then redshift = real_scalars[i-1].value
endfor
;;;;;

;;;;;
; get scalefactor of file
for i=1, (size(real_scalars))[3] do begin
    if (stregex(real_scalars[i-1].name, '^scalefactor', /BOOLEAN)) then scaleFactor = real_scalars[i-1].value
endfor

scale = 1.0/(1.0 + redshift)
redshiftConverted = 1.0/scaleFactor - 1.0    



;oneplusred = 1.0 + redshift

; get tree - 
read_amr_no_part, filename, VAR_NAME='pden', TREE=tree, DATA=dens, PARAMETERS=params

; get particle info - 

;read_write_particle_positions, filename, POSITIONS=positions, P_MASS=p_mass, $
;NUMBER_PARTICLES=number_particles



;      print, "redshift = ", redshift
;      lrefine = getRefineLevel(filename, 'dens')
;      print, "reading density..."
;      dens = loaddata_nomerge(filename,'dens')
;      ;temperature = loaddata_nomerge(filename,'temp')
    
;      ; variable name does not matter for xyz load coords
;      print, "reading x..."
;      x = load_coords(filename, 'dens', 1)
;      print, "reading y..."
;      y = load_coords(filename, 'dens', 2)
;      print, "reading z..."
;      z = load_coords(filename, 'dens', 3)
    
mass_enclosed_baryon = 0.0

max_x = 0.0
max_y = 0.0
max_z = 0.0

max_value = 0.0


; ENCLOSED BARYONIC MASS (not including sink particles)
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

print, "max density value and location"
print, max_value, max_x, max_y, max_z
print, "redshift"
print, redshiftConverted
print, "physical max density"
print, max_value*(1.0 + redshiftConverted)^3.0




end
