pro read_phase_data, cell_data, nx, ny, data_filename

openr,33,data_filename
for i = 0, nx-1 do begin
    for j= 0 , ny-1 do begin
        readf, 33, s
        cell_data[i,j] = s
        print, cell_data[i,j]
    endfor
endfor
close, 33

end
