function charis_get_phot_bandwidth,photfiltnames

bandwidth=fltarr(n_elements(photfiltnames))

for i=0L,n_elements(photfiltnames)-1 do begin

case photfiltnames[i] of 

   'Y': begin
         bandwidth[i]=0.0996
       end

   'J': begin
         bandwidth[i]=0.163
        end 
    
   'H': begin
         bandwidth[i]=0.296
        end

   'Ks': begin
         bandwidth[i]=0.311
         end


   '[3.1]':begin
          bandwidth[i]=0.1549
          end

   '[3.3]':begin
           bandwidth[i]=0.0555
           end


   'Lp': begin
          bandwidth[i]=0.7
         end

   '[4.05]': begin
           bandwidth[i]=0.068
            end


   'Ms': begin
           bandwidth[i]=0.241
         end
    else: begin
           print,'cannot find filter'
           stop
          end
endcase

endfor

return,bandwidth
end
