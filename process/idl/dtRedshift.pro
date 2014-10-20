pro dtRedshift, sstart, send, step

directory = "/scratch/01707/mepa/Rad_1Mpc/RadCosmoLW_res128"
file      = "/radCosmoLW_hdf5_chk_"

outfile = 'plots/timestep_' + String(strcompress(sstart, /remove)) + '_' + String(strcompress(send, /remove)) + '_' + String(strcompress(step, /remove)) + '.png'

yr = 3.1536e7
Myr = 1.0e6 * yr
kyr = 1.0e3 * yr

num_pts = floor((send-sstart)/step)+1
print, "num_pts:", num_pts
redshifts = dblarr(num_pts) 
timesteps = dblarr(num_pts) 
times = dblarr(num_pts) 

count = 0L
for number = sstart,send,step do begin
   
   if (number ge 1)   then prefix = '000'
   if (number ge 10)   then prefix = '00'
   if (number ge 100)  then prefix = '0'
   if (number ge 1000) then prefix = ''
   
   filename = directory + file + prefix + String(strcompress(number,/remove))
   print, filename
   
   if (number eq 228) then continue ;TACC lost this file
   ;;;;;
   ; get redshift, time and dt of file
   file_identifier = H5F_OPEN(filename)
   dataset = H5D_OPEN(file_identifier, "real scalars")
   real_scalars = H5D_READ(dataset)
   for i=1, (size(real_scalars))[3] do begin
      if (stregex(real_scalars[i-1].name, '^scale', /BOOLEAN)) then scale = real_scalars[i-1].value
      if (stregex(real_scalars[i-1].name, '^time', /BOOLEAN)) then time = real_scalars[i-1].value
      if (stregex(real_scalars[i-1].name, '^dt', /BOOLEAN)) then dt = real_scalars[i-1].value
   endfor
   ;;;;;
    
   redshift = (1.0 / scale) - 1.0
   
   print, "redshift = ", redshift
   redshifts[count] = redshift
   
   print, "time = ", time, "sec (", time/kyr, "kyr)"
   times[count] = time
   
   print, "timestep = ", dt, "sec (", dt/kyr, "kyr)"
   timesteps[count] = dt

   count = count + 1
   
endfor

xr = [10, 25]
yr = [5, 2000]

plot, redshifts, timesteps/kyr, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xrange=xr, yrange=yr, xtitle='redshift', ytitle='timestep (kyr)'

;plot, redshifts, times/Myr, /ylog, background='FFFFFF'xl, color=0, psym=-3, linestyle=0, xtitle='redshift', ytitle='age of universe (Myr)'

;oplot, redshifts, mean_baryon_densities/1.67e-24, color=1
    
void = cgSnapshot(FILENAME=outfile, /NoDialog)

end
