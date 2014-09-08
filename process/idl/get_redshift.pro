function get_redshift, filename


found_scalefactor = 0
file_identifier = H5F_OPEN(filename)
dataset = H5D_OPEN(file_identifier, "real scalars")
real_scalars = H5D_READ(dataset)
for i=1, (size(real_scalars))[3] do begin
    if (stregex(real_scalars[i-1].name, '^scalefactor', /BOOLEAN)) then begin
        scalefactor = real_scalars[i-1].value
        found_scalefactor = 1
    endif
endfor

if(found_scalefactor  eq 1) then begin
    
    redshift = 1.0 / scalefactor - 1.0

endif else begin
    
    redshift = 0.0
    
endelse

return, redshift


end
