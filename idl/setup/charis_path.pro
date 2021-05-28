function charis_path,pathname=pathname

;Basic structure is...
;directory_name = [path to CHARIS DPP]/[subdirectory]

;e.g. my path to CHARIS DPP is '~/idl_tools/ADI_dl/charis_dpp_new/'
;so for modeldirectory it is modeldirectory=[path to CHARIS DPP]/models/'
; ... or ... modeldirectory='~/idl_tools/ADI_dl/charis_dpp_new/models/'

;*****MODIFY THESE****

;1.
charispath='~/idl_tools/ADI_dl/charis_dpp_new/'

;these models have are from the online substellar atmosphere model server: 
;http://svo2.cab.inta-csic.es/theory/newov2/
;format lte[temp]-[logg]-[metallicity].BT-Sett.7.dat.txt
;e.g lte020-4.0-0.0.BT-Settl.7.dat.txt  = 2000K, log(g)=4, solar metallicity

;2.
btsettldirectory='/Users/thaynecurrie/Research/Planets/DI/ames/btsettl/'


;********


;****do NOT modify anything below this line

modeldirectory=charispath+'models/'
wavecaldirectory=charispath+'cals/'
filterresponsedirectory=charispath+'tools/filters/filter_response/'

planetdirectory=charispath+'models/planetmodels/'
modelspectrumdirectory=charispath+'models/starmodels/'

bonnefoydirectory=charispath+'models/planetmodels/bonn_speclib/'
montrealdirectory=charispath+'models/planetmodels/Montreal/'

;*not implemented yet, obsolete due to Cruz gravity standards?
;spexdirectory='~/Research/Planets/DI/spex/'

case pathname of

'modeldir': begin
 pathdirectory=modeldirectory
           end

'wavecalpath': begin
 pathdirectory=wavecaldirectory
            end

'planetdir': begin
 pathdirectory=planetdirectory
            end 

 'filtresponsedir': begin
pathdirectory=filterresponsedirectory
            end

'modelspecdir': begin
pathdirectory=modelspectrumdirectory
            end

'montrealdir': begin
pathdirectory=montrealdirectory
              end
'bonnefoydir': begin
pathdirectory=bonnefoydirectory
              end
'btsettldir': begin
pathdirectory=btsettldirectory
              end

;'spexdir': begin
;pathdirectory=spexdirectory
;           end

endcase

return,pathdirectory

end
