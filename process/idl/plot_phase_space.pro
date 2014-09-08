pro plot_phase_space, cell_data, fileout, x_min, x_max, y_min, y_max, nx, ny, trim, xtitle, ytitle, LINX=linear_x


; Plotting Options
charsize = 1.5
thick = 3



       xrange = [x_min, x_max]
       yrange = [y_min, y_max]


        if(linear_x ne 1) then begin
            plot, [0], /nodata, /noerase,  $
              xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
              _extra=extra_keywords, color=255, /xlog, /ylog, $
              xtitle=xtitle, ytitle=ytitle, $
              charsize=charsize, thick=thick, xthick=thick, ythick=thick, $
              charthick=thick, xtickformat='exponent', ytickformat='exponent'
        endif 
      
        if(linear_x eq 1) then begin
             plot, [0], /nodata, /noerase,  $
              xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
              _extra=extra_keywords, color=255, /ylog, $
              xtitle=xtitle, ytitle=ytitle, $
              charsize=charsize, thick=thick, xthick=thick, ythick=thick, $
              charthick=thick, ytickformat='exponent'
        endif


       ; make rectangles
       
       if(linear_x eq 1) then begin
           deltax = (x_max - x_min)/nx
       endif else begin
           deltax = (alog10(x_max) - alog10(x_min))/nx
       endelse
       deltay = (alog10(y_max) - alog10(y_min))/ny
       
       loadct, 54, file='mycolors.tbl'
       for i = 0, nx-1 do begin
           for j= 0 , ny-1 do begin
             
               xl = 10.0^(alog10((x_max / x_min)^(double(i) / nx) * x_min) + trim * deltax)
               xu = 10.0^(alog10((x_max / x_min)^(double(i+1) / nx) * x_min) - trim * deltax)
               
               if(linear_x eq 1) then begin
                   xl = double(i) * deltax + x_min + deltax * trim
                   xu = double(i+1) * deltax + x_min - deltax * trim
               endif
         
               yl = 10.0^(alog10((y_max / y_min)^(double(j) / nx) * y_min) + trim * deltay)
               yu = 10.0^(alog10((y_max / y_min)^(double(j+1) / nx) * y_min) - trim * deltay)
               
               xcell = [xl,xu,xu,xl]
               ycell = [yl,yl,yu,yu]
               
               if(cell_data[i,j] ge 1) then begin
                   polyfill, xcell, ycell, color=cell_data[i,j]
               endif
               
           endfor
       endfor    
  
       loadct, 12
       
       if(linear_x ne 1) then begin
           plot, [0], /nodata, /noerase,  $
             xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
             _extra=extra_keywords, color=0, /xlog, /ylog, $
             xtitle=xtitle, ytitle=ytitle, $
             charsize=charsize, thick=thick, xthick=thick, ythick=thick, $
             charthick=thick, xtickformat='exponent', ytickformat='exponent'
       endif 
       
       if(linear_x eq 1) then begin
           plot, [0], /nodata, /noerase,  $
             xrange=xrange, xstyle=1, yrange=yrange, ystyle=1, $
             _extra=extra_keywords, color=0, /ylog, $
             xtitle=xtitle, ytitle=ytitle, $
             charsize=charsize, thick=thick, xthick=thick, ythick=thick, $
             charthick=thick, ytickformat='exponent'
       endif
       


end
