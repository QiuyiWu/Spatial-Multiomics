# Rscript Code/assemble_manuscript_components.R COMPONENT_ROOT OUTPUT_DIRECTORY
# Figure 3: components 1/2 of A26, A27, A27 projected onto A26.
# S3: all 10 A26 components. S4: all 10 A27 components.
# Component PNGs are embedded unchanged; 600-dpi TIFF export does not add source detail.
args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
root <- normalizePath(args[1]);out <- normalizePath(args[2],mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
component_file <- function(group,i) file.path(root,paste0(group,'_components'),paste0(group,'_comp',i,'.png'))
files <- unlist(lapply(c('A26','A27','A27_to_A26'),function(g)vapply(1:10,function(i)component_file(g,i),character(1))))
stopifnot(all(file.exists(files)))
# Bound all nonwhite pixels, retaining labels, ticks, legends and a small margin.
bounds <- function(file) {
 im <- png::readPNG(file)
 ink <- pmin(im[,,1],im[,,2],im[,,3]) < .98
 if(dim(im)[3]==4L) ink <- ink & im[,,4]>.01
 rr <- which(rowSums(ink)>0);cc <- which(colSums(ink)>0)
 stopifnot(length(rr)>0,length(cc)>0)
 c(max(1,min(rr)-8),min(dim(im)[1],max(rr)+8),
   max(1,min(cc)-8),min(dim(im)[2],max(cc)+8))
}
all_bounds <- lapply(files,bounds);names(all_bounds)<-files
union_bounds <- function(fs) {
 b<-do.call(rbind,all_bounds[fs]);c(min(b[,1]),max(b[,2]),min(b[,3]),max(b[,4]))
}
group_bounds <- lapply(c('A26','A27','A27_to_A26'),function(g)union_bounds(vapply(1:10,function(i)component_file(g,i),character(1))))
names(group_bounds)<-c('A26','A27','A27_to_A26')
write.csv(data.frame(file=files,do.call(rbind,all_bounds)),file.path(out,'component_crop_bounds.csv'),row.names=FALSE)
width<-183;outer<-3;gap<-4;letter<-4;rowgap<-3
render <- function(name,groups,components,heading=NULL,rowlabels=NULL) {
 nr<-nrow(components);nc<-ncol(components)
 label_width<-if(is.null(rowlabels))0 else 12
 heading_height<-if(is.null(heading))0 else 9
 cell_width<-(width-2*outer-label_width-(nc-1)*gap)/nc
 image_width<-cell_width-letter
 ratios<-vapply(groups,function(g){b<-group_bounds[[g]];(b[2]-b[1]+1)/(b[4]-b[3]+1)},numeric(1))
 image_height<-image_width*max(ratios)
 height<-2*outer+heading_height+nr*image_height+(nr-1)*rowgap
 draw<-function(){
  grid::grid.newpage();grid::pushViewport(grid::viewport(xscale=c(0,width),yscale=c(0,height)))
  if(!is.null(heading))grid::grid.text(heading,x=.5,y=grid::unit(height-outer-heading_height/2,'native'),gp=grid::gpar(fontsize=12,fontface='bold'))
  for(r in seq_len(nr)) {
   top<-height-outer-heading_height-(r-1)*(image_height+rowgap)
   if(!is.null(rowlabels))grid::grid.text(rowlabels[r],x=grid::unit(outer+label_width/2,'native'),y=grid::unit(top-image_height/2,'native'),gp=grid::gpar(fontsize=11,fontface='bold'))
   for(c in seq_len(nc)){
    group<-if(length(groups)==1L)groups else groups[r]
    f<-component_file(group,components[r,c]);im<-png::readPNG(f);b<-group_bounds[[group]]
    im<-im[b[1]:b[2],b[3]:b[4],,drop=FALSE]
    x<-outer+label_width+(c-1)*(cell_width+gap)
    h<-image_width*(b[2]-b[1]+1)/(b[4]-b[3]+1)
    grid::grid.raster(im,x=grid::unit(x+letter,'native'),y=grid::unit(top,'native'),just=c('left','top'),width=grid::unit(image_width,'mm'),height=grid::unit(h,'mm'),interpolate=FALSE)
    grid::grid.text(letters[(r-1)*nc+c],x=grid::unit(x,'native'),y=grid::unit(top-.5,'native'),just=c('left','top'),gp=grid::gpar(fontsize=11,fontface='bold'))
   }
  };grid::popViewport()
 }
 tiff(file.path(out,paste0(name,'.tiff')),width=width,height=height,units='mm',res=600,compression='lzw',type='cairo');draw();dev.off()
 png(file.path(out,paste0(name,'_preview.png')),width=width,height=height,units='mm',res=150,type='cairo');draw();dev.off()
 cat(name,width,height,'mm\n')
}
render('Figure3',c('A26','A27','A27_to_A26'),matrix(rep(1:2,3),nrow=3,byrow=TRUE),rowlabels=c('A26','A27','A27\n=>\nA26'))
render('Figure S3','A26',matrix(1:10,nrow=5,byrow=TRUE),heading='A26')
render('Figure S4','A27',matrix(1:10,nrow=5,byrow=TRUE),heading='A27')
writeLines(c('Figure3: row 1 A26 components 1/2; row 2 A27 components 1/2; row 3 A27-to-A26 components 1/2.',
             'Figure S3: A26 components 1 through 10 in row-major order.',
             'Figure S4: A27 components 1 through 10 in row-major order.',
             'Paired maps retain their source PNG signs, labels, orientation, and legends.'),file.path(out,'component_assembly.txt'))
