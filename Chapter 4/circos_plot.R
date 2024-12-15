modules1UP <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/negative_weights_factor3/modules1.rds")
colnames(modules1UP[["corrplot"]][["corr"]])=paste0('F3_UP_',colnames(modules1UP[["corrplot"]][["corr"]]))

modules1DOWN <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/positive_weight_factor3/modules1.rds")
colnames(modules1DOWN[["corrplot"]][["corr"]])=paste0('F3_DOWN_',colnames(modules1DOWN[["corrplot"]][["corr"]]))

cell_type1UP <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/negative_weights_factor3/cell_type1.rds")
cell_type1DOWN <- readRDS("/rds/general/project/ukdrmultiomicsproject/live/Proteomics/synapticproteomics/synaptic_cytosolic/2_synaptic_cytosolic_integration/positive_weight_factor3/cell_type1.rds")


## mat corr
library(circlize )


mat1=t(cbind(modules1UP[["corrplot"]][["corr"]],
             modules1DOWN[["corrplot"]][["corr"]]))

col_fun1 = colorRamp2(c(-0.7, 0, 0.7), c("#2166AC", "white", "#B2182B"))
#mat 0

# mat modules
mat0=data.frame(rownames(mat1))
rownames(mat0)=rownames(mat1)
#mat0$rownames.mat1.=seq(1:length(mat0$rownames.mat1.))
#rownames(mat0)=str_replace_all(rownames(mat0), "ME", "")
mat0$module=rownames(mat0)

modNames=rownames(mat1)

col_fun0=list()
for(i in modNames){
  col_fun0[i]<-i
}
names(col_fun0)<-substring(rownames(mat1), 3)


### matrix for cell types p values
list_cell=list()
list_cell_types=list(cell_type1UP,cell_type1DOWN)
names(list_cell_types)=c("cell_type1UP","cell_type1DOWN")

library(readr)
list_cell=list()
i="cell_type1UP" 
for(i in names(list_cell_types)){
  df <- do.call("rbind", list_cell_types[[paste(i)]][['plots']])
  df=df[,1]
  
  df <- do.call("rbind", df)
  
  df$color=sub("\\..*", "", rownames(df))
  
  
  df2=df[,c('q','CellType','color')]
  df3=reshape(df2, idvar = "color", timevar = "CellType", direction = "wide")
  rownames(df3)=df3$color
  df3$color=NULL
  df3=abs(df3)
  rownames(df3)=paste0('F',parse_number(i),'_','ME',rownames(df3))
  colnames(df3)=str_remove(colnames(df3), 'q.')
  list_cell[[paste(i)]]=df3
}



pvalues_list= do.call("rbind", list_cell)

# Define a function to reformat the strings
reformat_string <- function(input_string) {
  # Extract the relevant parts using regex
  result <- gsub("cell_type(\\d+)(UP|DOWN)\\.F(\\d+)\\_ME(\\d+)", "F\\3_\\2_ME\\4", input_string)
  return(result)
}

# Apply the function to each element in the vector
reformatted_strings <- sapply(rownames(pvalues_list), reformat_string)


rownames(pvalues_list)=as.character(reformatted_strings )     
rownames(pvalues_list)=c("F3_UP_ME1"   ,"F3_UP_ME2" ,  "F3_DOWN_ME1" ,"F3_DOWN_ME3" ,"F3_DOWN_ME4")
pvalues_list["F3_DOWN_ME2",]=c(NA,NA,NA,NA,NA,NA,NA,NA,NA)
pvalues_list=pvalues_list[rownames(mat1),]


col_fun2 = colorRamp2(c(0,0.05 ,1), c("black",'yellow',"white"))

## mat1 and pval_mat 

#mat 0
#circos.par(start.degree = 30, gap.degree = 1, points.overflow.warning = FALSE)

mat00=data.frame(rownames(mat1))
rownames(mat00)=rownames(mat1)


add_significance_stars <- function(p_values, row_order, col_order) {
  for (i in 1:length(row_order)) {
    for (j in 1:length(col_order)) {
      p_val = p_values[i, j]
      if (p_val < 0.001) {
        stars_text <- "***"
      } else if (p_val < 0.01) {
        stars_text <- "**"
      } else if (p_val < 0.05) {
        stars_text <- "*"
      } else {
        stars_text <- ""
      }
      
      if (stars_text != "") {
        cell_center_x <- (i - 0.5)  # Calculate center x-coordinate of the cell
        cell_center_y <- (((length(col_order)-j)+1) - 0.6)  # Calculate center y-coordinate of the cell
        
        circos.text(cell_center_x, cell_center_y, stars_text, col = "black", cex = 0.8)
      }
    }
  }
}


pval_mat_label=rbind(modules1UP$moduleTraitPvalue,
                     modules1DOWN$moduleTraitPvalue)
rownames(pval_mat_label)=rownames(mat1)


### PLOT
circos.clear()
circos.par( gap.after = c(70))


circos.heatmap(mat00, 
               col = col_fun0, 
               track.height = 0.05, 
               bg.border = "gray50",
               rownames.side = "outside",
               rownames.cex = 0.8,
               colnames.cex=0.8)
circos.heatmap(mat1[,c(1,2,3,5,6)], col = col_fun1,track.height = 0.3,cell.border='black',cluster = TRUE)

# Add significance stars to the mat1 heatmap
add_significance_stars(p_values = pval_mat_label[,c(1,2,3,5,6)],
                       row_order = rownames(mat1[,c(1,2,3,5,6)]),
                       col_order = colnames(mat1[,c(1,2,3,5,6)]))

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = rev(colnames(mat1[,c(1,2,3,5,6)]))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), #x coordinate
                (1:n),#y coordinate
                cn, #label
                cex = 1, adj = c(0, 1), facing = "inside")
  }
}, bg.border = NA)


colnames(pvalues_list)[9]='Vascular'
circos.heatmap(pvalues_list, col = col_fun2,track.height = 0.3,cell.border='grey')
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  if(CELL_META$sector.numeric.index == 1) { # the last sector
    cn = rev(colnames(pvalues_list))
    n = length(cn)
    circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), #x coordinate
                (1:n),#y coordinate
                cn, #label
                cex = 0.8, adj = c(0, 1), facing = "inside")
  }
}, bg.border = NA)


lgd = Legend(title = "Pearson's correlation", col_fun = col_fun1,at=c(-0.7,0,0.7))
lgd2 = Legend(title = "BH-corrected q-values", col_fun = col_fun2,at=c(0,0.05,1))
lgd_list_vertical2 = packLegend(lgd,lgd2)

draw(lgd_list_vertical2, x = unit(50, "mm"), y = unit(1, "mm"), just = c("right", "bottom"))



