pro jeans_criteria, chk=number

; File Options

directory = "/scratch/01707/mepa/Rad_res512/"
file = "rad_hdf5_chk_"

num       = 0043
sstart    = 206
send      = 235
step      = 1
box_size  = 3.08D24

kb = 1.38e-16
newton = 6.67e-8
mh = 1.67e-24
mu = 1.2
gamma = 5.0 / 3.0
yhe = 0.08
pi = 3.141592



; Output Options
plotps = 0
plotgif = 1
charsize = 2.0
charsize1 = 1.0
thick = 2
exthick = 4


    if (number ge 0)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''

    filename = directory + file + prefix + String(strcompress(number,/remove))
    print, filename

    redshift = get_redshift(filename)
    oneplusred = 1.0 + redshift
    
       
    ; get tree and data
    read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
    read_amr_no_part, filename, VAR_NAME='temp', TREE=tree, DATA=temp, PARAMETERS=params
    dens = reform(dens)
    temp = reform(temp)
    
    dens = dens * (1.0 + redshift)^3.0
    temp = temp * (1.0 + redshift)^(-2.0)

    numdens = dens / (mh * (1.0 + 4.0 * yhe))
    sound_speed = sqrt(gamma * kb * temp / (mh * mu))
    jeans_length = sqrt(pi * sound_speed^2 / (newton * dens) )
    
   
    
    max_lrefine = max(tree[*].lrefine)

    print, "max lrefine = ", max_lrefine
    
    coords = locate_dens_max(dens,tree,params)
    
    centerx = coords[0]
    centery = coords[1]
    centerz = coords[2]

    r_max = 1.0e23
    
    range_xl = centerx - r_max
    range_xu = centerx + r_max
    range_yl = centery - r_max
    range_yu = centery + r_max
    range_zl = centerz - r_max
    range_zu = centerz + r_max
    
    good_blocks = where(tree.nodetype EQ 1 AND $
                        tree.bndbox[0,0] LE range_xu AND tree.bndbox[1,0] GT range_xl AND $
                        tree.bndbox[0,1] LE range_yu AND tree.bndbox[1,1] GT range_yl AND $
                        tree.bndbox[0,2] LE range_zu AND tree.bndbox[1,2] GT range_zl, cblocks)

    good_blocks = where(tree.nodetype EQ 1, cblocks)
    
    
    plot, [0], /nodata, /xlog, /ylog, xrange=[1.0e-4, 1.0e9], yrange=[1.0, 1.0e4], color=0, background='FFFFFF'xl, thick=thick, charthick=thick, xthick=thick, ythick=thick, charsize=charsize
    
    
    for block_id1 = 0.0, cblocks - 1, 1 do begin
        
        block_id = good_blocks[block_id1]
        
        ; physical x cell size
        dx = tree[block_id].size[0] / params.nxb / oneplusred
        block_lj = jeans_length[block_id,*,*,*]/dx

        oplot, numdens[block_id,*,*,*], block_lj, color=0, psym=3

    endfor

    oplot, [1.0e-5, 1.0e20], [24.0, 24.0], color=0, linestyle=2
    oplot, [1.0e-5, 1.0e20], [48.0, 48.0], color=0, linestyle=2
    
end


