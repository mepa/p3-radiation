function loaddata_nomerge, filename, var, $
                   SAMPLE=sample, $
                   DOUBLE=double, $
                   XCOORDS=x, YCOORDS=y, ZCOORDS=z, $
                   XLCOORD=xl, XRCOORD=xr, $
                   XRANGE=xrange, YRANGE=yrange, ZRANGE=zrange, $
                   TIME=time, UNIFORM_1D=uniform_1d

;
; simple proceedure to read a single variable (or derived variable), var,
; from the file, filename.  The data is read in, and put onto a
; uniformly gridded mesh at the resolution of the finest AMR mesh (or 
; subsampled by sample levels if desired).  The coordinates of the
; data are returned through the optional keywords, XCOORDS, YCOORDS,
; and ZCOORDS.
;
; USAGE
; 
; Here is the basic way to read a file in and grab the energy.
;
; energy = loaddata('myCheckPointFile','ener')
; 
; optionally
; 
; energy = loaddata('myCheckPointFile','ener', XCOORDS=x, YCOORDS=y,
; TIME=t)
;



if (n_elements(filename) EQ 0) then begin
    print, 'ERROR: no filename specified to get_var_list'
    return, -1
endif

if n_elements(var) EQ 0 then begin
    print, 'ERROR: no variable specified'
    return, -1
endif

if n_elements(sample) EQ 0 then sample = 0

if n_elements(uniform_1d) EQ 0 then uniform_1d = 0

if n_elements(double) EQ 0 then double = 0

;------------------------------------------------------------------------------
; read in the data
;------------------------------------------------------------------------------
itype = determine_file_type(filename)

if (double) then begin
    read_amr_no_part, filename, VAR_NAME=var, $
      TREE=tree, DATA=unk, PARAMETERS=params
endif else begin
    read_amr_no_part, filename, VAR_NAME=var, $
      TREE=tree, DATA=unk, PARAMETERS=params
endelse

time = params.time

; set the ranges
if n_elements(xrange) EQ 0 then begin
    xrange = fltarr(2)
    xrange[0] = min(tree[*].bndBox[0,0])     ; set to minimum of x coord
    xrange[1] = max(tree[*].bndBox[1,0])     ; set to maximum of x coord
endif

if (params.ndim GE 2) then begin
    if n_elements(yrange) EQ 0 then begin
        yrange = fltarr(2)
        yrange[0] = min(tree[*].bndBox[0,1]) ; set to minimum of x coord
        yrange[1] = max(tree[*].bndBox[1,1]) ; set to maximum of x coord
    endif
endif

if (params.ndim EQ 3) then begin
    if n_elements(zrange) EQ 0 then begin
        zrange = fltarr(2)
        zrange[0] = min(tree[*].bndBox[0,2]) ; set to minimum of x coord
        zrange[1] = max(tree[*].bndBox[1,2]) ; set to maximum of x coord
    endif
endif

sData = reform(unk, params.totBlocks, $
                                          params.nxb, $
                                          params.nxb, $
                                          params.nzb)
   

;------------------------------------------------------------------------------
; merge it onto a uniform grid

 


return, sData
end


