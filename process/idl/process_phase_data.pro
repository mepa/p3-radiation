pro process_phase_data, cell_data, cmax, cmin

  zero_cells = where(cell_data le 0.0, num_zero)
  if(num_zero ge 1) then begin
      cell_data[zero_cells] = 0.0
  endif
  
  good_cells = where(cell_data gt 0.0, num_good)
  if(num_good eq 0) then begin
      print, "[process_phase_data]: no good cells found!!!!"
  endif
  minval = min(cell_data[good_cells])
  maxval = max(cell_data[good_cells])
  print, "Min and max of cell data=", minval, maxval
  ;; Enforces that color range only extends
  ;; 5 dex - purely stylistic:
  minval = max(maxval-7.0, minval)
  cell_data = cmin + bytscl(temporary(cell_data), min=minval, max=maxval,$
                                   top=cmax-cmin)

end
