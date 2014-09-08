FUNCTION locate_var_min, var, tree, params

; Goal: Find the comoving position of the maximum value of var
;
; var - reformed variable variable from read_amr_no_part
; tree - tree from read_amr_no_part
; params - params from read_amr_no_part

  print, "Locating maximum of variable... "

  var_min = min(var, min_index)
  var_min_indices = array_indices(var, min_index)
  lb = var_min_indices[0]
  i  = var_min_indices[1]  
  j  = var_min_indices[2]
  k  = var_min_indices[3]
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
  
  min_coords = [x_center,y_center,z_center, dx]
  print, "var min coords=", x_center, y_center, z_center
  print, "var min value=", min(var)
  
  return, min_coords



END
