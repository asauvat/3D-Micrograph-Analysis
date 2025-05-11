#########################################################
# ~3D spheroid segmentation's tutorial~                 #
#                                                       #
# DIU "modélisation multi-échelles complexité tumorale "#
# Author : Allan Sauvat                                 #
# Date : 16/05/25                                       #
#########################################################

#=================================================================================================================================================================#
# import Python modules

library(reticulate)
use_condaenv('r-reticulate')
#
np = import('numpy', as = 'np', convert = F)
sk = import('skimage', as ='sk', convert = F)
scp = import('scipy', as = 'scp', convert = F)
tfl = import('tifffile', as = 'tfl', convert = F)
nap = import('napari', as = 'nap', convert = F)

#=================================================================================================================================================================#
# Load required R libraries

library(EBImage)
library(RBioFormats)
#
library(magrittr)
library(xml2)

#=================================================================================================================================================================#
# Get images information

ips = list.files('MICROGRAPHS', full.names = T) # list file names inf image folder and store in array
#
met =  lapply(ips,read.metadata);names(met) = basename(ips)
ome = lapply(ips,function(ipi)read.omexml(ipi));names(ome) = basename(ips) # read OME-metadata as xml
ome = lapply(ome,function(x)read_xml(x) %>% as_list()) # parse OME-XML
#
str(ome[[1]]) # overview of OME-XML structure, after conversion to list
str(met[[1]]$coreMetadata) # overview of base image metadata
#
ps = sapply(c('PhysicalSizeZ','PhysicalSizeY','PhysicalSizeX'),function(s)attr(ome[[1]]$OME$Image$Pixels,s)) %>% as.numeric # get pixel size attributes, and convert to numeric
nc = met[[1]]$coreMetadata$sizeC # number of channels in images
cn = sapply(1:nc,function(i)attr(ome[[1]]$OME$Image$Pixels[[i]],'Name'))# Get channel names

#=================================================================================================================================================================#
# Segment images

batch_mode = TRUE # control intermediate display
if(!dir.exists('SEGS')){dir.create('SEGS')} # Where segmentations will be saved
if(!dir.exists('CELLDATA')){dir.create('CELLDATA')} # Where cell data will be saved
#
dat = list() # initialize list that will contain cell data
count = 1 # counter for loop
#
for(ipi in ips){
  print(paste0('Analyzing image:',basename(ipi),',',count,'/',length(ips)))
  
  # Import image stack in environment------------------#
  #----------------------------------------------------#
  
  img = sk$io$imread(ipi)
  print(img$shape) # Z,C,Y,X
  
  # Vizualize interactively in 3D----------------------#
  #----------------------------------------------------#
  
  if(!batch_mode){
    vw = nap$Viewer()
    vw$add_image(img, scale = ps, channel_axis = 1L)
    nap$run()
    rm(vw);gc()
  }
  
  # 3D binarization ---------------------------------#
  #-------------------------------------------------#
  
  kball = sk$morphology$ball(5L);kball = kball[c(3L,5L,9L),,] # Create Structuring Element
  #
  mk = img[,1L,,] > sk$filters$threshold_local(img[,1L,,], block_size = c(11L,25L,25L), method = 'mean', offset=-2000) #3D adaptive threshold
  mk = sk$morphology$binary_closing(mk, footprint = kball) # Morphology : dilation followed by erosion
  mk = sk$morphology$binary_opening(mk, footprint = kball) # Morphology : erosion followed by dilation
  #
  if(!batch_mode){
    vw = nap$Viewer()
    vw$add_image(img, scale = ps, channel_axis = 1L)
    vw$add_labels(mk, scale = ps)
    nap$run()
    rm(vw);gc()
  }
  
  # Nuclei labelization-----------------------------#
  #-------------------------------------------------#
  
  # Compute distance transform for all slices at once
  dst = scp$ndimage$distance_transform_edt(mk, sampling=ps) 
  
  # Find local maxima [internal marker]
  lmax = sk$feature$peak_local_max(dst, min_distance=15L, labels=mk, footprint=np$ones(c(3L,3L,3L)), exclude_border = TRUE, p_norm=2L)   
  if(!batch_mode){
    View(py_to_r(lmax))
  }
  
  # Create markers for watershed
  markers = np$zeros_like(mk, dtype='bool') # stack filled with 0s
  markers[lmax[,0L], lmax[,1L], lmax[,2L]] = TRUE  # set TRUE where local maxima are located
  
  # Apply dilation and label
  markers = scp$ndimage$label(sk$morphology$binary_dilation(markers, np$ones(c(3L,5L,5L))),structure=np$ones(c(3L,3L,3L)))[0] 
  
  # Apply watershed segmentation
  labs = sk$segmentation$watershed(-dst, markers=markers, mask=mk, watershed_line=TRUE)
  
  # Remove small objects and relabel
  labs = sk$segmentation$relabel_sequential(sk$morphology$remove_small_objects(labs, 2000), offset=1)[0]   
  
  if(!batch_mode){
    vw = nap$Viewer()
    vw$add_image(img, scale = ps, channel_axis = 1L)
    vw$add_labels(labs, scale = ps)
    nap$run()
    rm(vw);gc()
  }
  
  # Write segmentation file 
  clb = sk$color$label2rgb(labs) # convert labeled matrix to RGB, ZYXC
  clb = np$transpose(clb, c(0L,3L,1L,2L)) # ZYXC to ZCYX
  #
  tif = tfl$TiffWriter(gsub('MICROGRAPHS','SEGS',ipi),bigtiff=FALSE,ome = TRUE);# create writer
  targs = list(photometric='minisblack', compression='lzw',metadata=list(Pixels = list(PhysicalSizeX = ps['x']*10**-3, PhysicalSizeXUnit = 'mm', #here we could replace "." with "," for local system compatibility
                                                                                       PhysicalSizeY = ps['y']*10**-3, PhysicalSizeYUnit = 'mm',
                                                                                       PhysicalSizeZ = ps['z']*10**-3, PhysicalSizeZUnit = 'mm'),
                                                                         axes='ZCYX'))
  do.call(tif$write,c(list(data = clb), targs))
  tif$close()
  
  # Features computing------------------------------#
  #-------------------------------------------------#
  
  poi = c("centroid","num_pixels","solidity","intensity_mean","label") #new features
  NFi = sk$measure$regionprops_table(label_image = labs, intensity_image=np$transpose(img,c(0L,2L,3L,1L)), properties = poi,spacing = rev(ps)) %>% py_to_r() %>% do.call(what='cbind') %>% data.frame()
  #
  if(!batch_mode){
    View(NFi)
  }
  write.csv(NFi, paste0('CELLDATA/',gsub('.tiff','.csv',basename(ipi))), row.names = FALSE)
  #
  dat = c(dat,list(NFi)) # append computing feature to list
  count = count+1 # increment counter
}
names(dat) = basename(ips) # name data according to images
rm(list = setdiff(ls(),c('ips','cn','ps','dat'))) # clean environment
gc()

#=================================================================================================================================================================#
# Basic data processing

# Retrieve unique treatment groups------------------#
#---------------------------------------------------#
groups = sapply(strsplit(names(dat),'_'),function(x)x[[1]]) # give group labels from file names
print(groups)
uni_gp = unique(groups) # get unique labels

# Quick view on parameters of interest--------------#
#---------------------------------------------------#
poi = c('num_pixels','intensity_mean.0') # variables of interest : nuclear volume and GFP intensity
#
par(mfrow=c(1,2)) # create grid for plotting multiple graphs
for(p in poi){
  boxplot(lapply(dat,function(x)x[,p]), outline = FALSE, main = p, las = 1)
}
# Aggregate data for graphical interpreation---------#
#---------------------------------------------------#
sdat = lapply(uni_gp,function(gpi){
  y = dat[grep(gpi,names(dat))] %>% do.call(what = rbind) # extract element from dat list, that contain the groupID pattern
  y = y[,poi] # select columns of interest
  return(y)
})
names(sdat) = uni_gp

# Graphical interpretation-------------------------#
#---------------------------------------------------#
plot(log(sdat$CTR), pch = 19, col = scales::alpha('blue',0.5), xlab = 'Nuclear volume [log]', ylab = 'GFP intensity [log]', main = 'Cell clustering',bty = 'n')
points(log(sdat$OXA), pch = 19, col = scales::alpha('red',0.5))
legend(x = 'topright', bty = 'n', legend = names(sdat), fill = c('blue','red'))
th = locator(n=1) # prompt user to interactively click on plot to choose threshold
abline(v=th$x, h = th$y, col = 'black', lty = 2, lwd = 2, xpd=F)

# Cell data labeling---------------------------------#
#---------------------------------------------------#

dat = lapply(dat,function(x){
  x$isSmall = log(x[,poi[1]]) < th$x # could help detecting pyknotic cells
  x$isGFP = log(x[,poi[2]]) > th$y # nuclei which are positive for the GFP marker
  return(x)
})

# Data reducing--------------------------------------#
#---------------------------------------------------#

doi = c('isSmall','isGFP')
#
red = lapply(dat,function(x){
  y = lapply(doi,function(d)table(x[,d])) %>% do.call(what = rbind)
  rownames(y) = doi
  return(y)
})






