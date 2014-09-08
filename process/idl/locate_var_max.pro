FUNCTION locate_var_max, var, tree, params

; Goal: Find the comoving position of the maximum value of var
;
; var - reformed variable variable from read_amr_no_part
; tree - tree from read_amr_no_part
; params - params from read_amr_no_part

  print, "Locating maximum of variable... "

  var_max = max(var, max_index)
  var_max_indices = array_indices(var, max_index)
  lb = var_max_indices[0]
  i  = var_max_indices[1]  
  j  = var_max_indices[2]
  k  = var_max_indices[3]
  ;print, lb, i, j, k
  bxl = tree[lb].bndBox[0,0]
  byl = tree[lb].bndBox[0,1]
  bzl = tree[lb].bndBox[0,2]
  dx = tree[lb].size[0] / params.nxb
  dy = tree[lb].size[1] / params.nyb
  dz = tree[lb].size[2] / params.nzb
  xc = dindgen(params.nxb) * dx + bxl + dx/2.0
  yc = dindgen(params.nyb) * dy + byl + dy/2.0
  zc = dindgen(params.nzb) * dz + bzl + dz/2.0
  x_center = xc[i]
  y_center = yc[j]
  z_center = zc[k]
  
  max_coords = [x_center,y_center,z_center, dx]
  print, "var max coords=", x_center, y_center, z_center
  print, "var max value=", max(var)
  
  return, max_coords



END
