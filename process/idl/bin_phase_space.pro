pro bin_phase_space, xbin, ybin, dens, xmin, xmax, ymin, ymax, nx, ny, good_blocks, cblocks, tree, params, CELLDATA=cell_data

print, "Starting phase binning.."

ycoeff = alog10(ymax/ymin)
xcoeff = alog10(xmax/xmin)
m=ceil(cblocks/10)
for blk = 0.0, cblocks - 1, 1 do begin
    block_id = good_blocks[blk]
    dx = tree[block_id].size[0] / params.nxb 
    dy = tree[block_id].size[1] / params.nyb 
    dz = tree[block_id].size[2] / params.nzb 
    for ii=0, params.nxb-1 do begin
        for jj=0, params.nyb-1 do begin
            for kk=0, params.nzb-1 do begin
                x = xbin[block_id,ii,jj,kk]
                y = ybin[block_id,ii,jj,kk]
                if ((x gt xmin) and (x lt xmax)) then begin
                    if ((y gt ymin) and (y lt ymax)) then begin
                        y_bin = ceil(ny * alog10(y / ymin) / ycoeff) - 1
                        x_bin = ceil(nx * alog10(x / xmin) / xcoeff) - 1
                        weight = dens[block_id,ii,jj,kk]*dx*dy*dz
                        cell_data[x_bin, y_bin] += weight
                    endif
                endif
            endfor
        endfor
    endfor
    if ((blk mod m) EQ 0) then print, "binning phase space...", blk, cblocks
endfor

cell_data = alog10(cell_data)


end
