pro write_phase_data, cell_data, nx, ny, data_filename

openw,33,data_filename
for i = 0, nx-1 do begin
    for j= 0 , ny-1 do begin
        s = cell_data[i,j]
        printf, 33, s
    endfor
endfor
close, 33

end
