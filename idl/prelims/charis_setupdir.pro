pro charis_setupdir

;04/27/2019 - changed to simplify.
;03/05/2018 - copied from setupdir to package with CHARIS post-processing pipeline.
;does directory structure setup

file_mkdir,'data'
file_mkdir,'reduc'
file_mkdir,'reduc/prep'
file_mkdir,'reduc/reg'
file_mkdir,'reduc/rsub'
file_mkdir,'reduc/proc'

file_mkdir,'data/raw'

file_mkdir,'psfmodel'

;should do this at some point: make a log directory just like GPI DRP

file_mkdir,'log'
;for future implementation

end
