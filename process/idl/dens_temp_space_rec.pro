 pro dens_temp_space_rec, chk=number

; File Options

directory = "/scratch/01707/mepa/Rad_res512/"
file = "rad_hdf5_chk_"

kb = 1.38d-16
mh = 1.67d-24
me = 9.11d-28
mu = 1.4 
yhe = 0.08
gravity = 6.67e-8
pi = 3.141592
H0 = 2.279e-18
omega_b = 0.0458
colormap_min = 0
colormap_max = 254
start_time = systime(/seconds)


;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;
; User Parameters

; y = temp
; x = yn
    ;x_min = 1.0d-24
    ;x_max = 1.0e-15
    x_min = 1.0d-4*mh
    x_max = 1.0d+4*mh
    y_min = 10.0
    y_max =20000.0
    nx = 200
    ny = 200
    trim = 0.00
    cell_max_mass = 1.0d40

; Make the x-axis linear scaling? (may not work...)
    linear_x = 0

    write_location = ''
    file_start = 'dts_cutout_' 
    xtitle = 'Density [g cm!U-3!N]'
    ytitle = 'Temperature [K]'

    ; option 1 = do from scratch
    ; option 2 = work w/ already made data file
    option = 1

 
;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;

    if (number ge 0)   then prefix = '000'
    if (number ge 10)   then prefix = '00'
    if (number ge 100)  then prefix = '0'
    if (number ge 1000) then prefix = ''
    strnum =  prefix + String(strcompress(number,/remove))
    filename = directory + file + strnum
    fileout = write_location + file_start + strnum + '.eps'
    data_filename = write_location + file_start + strnum + '.dat'
   
    redshift = get_redshift(filename)
    oneplusred = 1.0 + redshift
    print, "processing file ", filename
    print, "redshift = ", redshift

    T_cmb = 2.725 * oneplusred
    
    cell_data  = dblarr(nx, ny)
    
    if(option eq 1) then begin
        
        ;;; Read in necessary data
        read_amr_no_part, filename, VAR_NAME='dens', TREE=tree, DATA=dens, PARAMETERS=params
        read_amr_no_part, filename, VAR_NAME='temp', TREE=tree, DATA=temp, PARAMETERS=params
        read_amr_no_part, filename, VAR_NAME='gpot', TREE=tree, DATA=gpot
        dens = reform(dens)
        temp = reform(temp)
        gpot = reform(gpot)
        
        ;;; Create desired variables
        dens = dens * (1.0 + redshift)^3.0
        temp = temp / (1.0 + redshift)^(2.0)
    
        yn = dens / (mh * (1.0 + 4.0 * yhe))
     
        ;;; Get blocks in range
        coords = locate_var_max(dens,tree,params)
        print, "centering at=", coords
        coords1 = locate_var_min(gpot,tree,params)
        print, "gpot min=", coords1
        
        r_max = 50.0 * 3.08d18
        r_max = r_max * (1.0 + redshift)
        get_good_blocks, coords, r_max, tree, GOODBLOCKS=good_blocks, CBLOCKS=cblocks
        
        
        ;;; Do the binning:
        bin_phase_space, dens, temp, yn, x_min, x_max, y_min, y_max, nx, ny, good_blocks, $
          cblocks, tree, params, CELLDATA=cell_data
        
        ;;; Read out cell_data
        write_phase_data, cell_data, nx, ny, data_filename
     
  endif

  
  ; If reading in cell_data
  if(option eq 2) then begin
      read_phase_data, cell_data, nx, ny, data_filename
  endif
  
 
  ;; Make the image
  
  process_phase_data, cell_data, colormap_max, colormap_min
  

  set_plot, 'PS'
  device, /color, bits_per_pixel=8, filename=fileout, encapsulated=1
  loadct, 12
  device, xsize=20

  plot_phase_space, cell_data, fileout, x_min, x_max, y_min, y_max, nx, ny, trim, xtitle, ytitle, LINX=linear_x
  
  
;;; Supplementary stuff:
 

      ;; lines of constant jeans mass:
      jeans_masses = [1.0e3,1.0e2,1.0e1,1.0, 0.1]
      njeans = size(jeans_masses)
      numdens_jeans = findgen(50) - 20.0
      numdens_jeans = 10.0^(numdens_jeans)
      print, "number of jeans masses=", njeans[1]
      for i=0,njeans[1]-1 do begin
          ;print, "making mj curve for mj =", jeans_masses[i]
          ;oplot, numdens_jeans, numdens_jeans^(1.0/3.0) * (jeans_masses[i] / (20.0 * 1.9))^(2.0/3.0), color=0
      endfor

      ;; CMB temperature line

      oplot, [1.0e-30, 1.0e30], [T_cmb, T_cmb], color=0, linestyle=2


; sink threshold
      jeans_dens_coeff = pi / (4.0^2.0 * gravity)
      temp = [10.0,1.0e2, 1.0e3, 1.0e4, 1.0e5]
      dx2min = (3.08568d24 / 8.0 / 2.0^(19.0) / (1.0 + redshift))^2
      cs2 =  kb * temp / mh
      oplot, jeans_dens_coeff * cs2 / dx2min, temp, color=0
      print, jeans_dens_coeff * cs2 / dx2min
      print, temp
  
      

      
      print, "wrote file ", fileout
      print, 'Time to run = ', systime(/seconds) - start_time
      device, /close
      set_plot, 'x'

  end
  
