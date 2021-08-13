library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
options(stringsAsFactors = FALSE)

source("heatmap.3.smallLabels.R")
source("misc_rnaseq_funcs.R")

matrix.file <- "Pacuta_Pos_Met.csv"

## load matrix file
print('Reading matrix file.')
data1 = read.table(matrix.file, header=T, com='', row.names=1, check.names=F, sep=',')
data = as.matrix(data1)

myheatcol = colorpanel(75, 'purple','black','yellow')
samples_data = read.table("Pacuta_Sample_Info_Treatment.csv", h=T, sep=",")
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name', 'TP', 'Weight')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% rep_names, drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
  samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
  sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}


## Normalize data
#data = data[rowSums(data)>=10,]
#cs = colSums(data)
#data = t( t(data)/cs) * 1e6; ## CPM
mode(data)<-'numeric'
data = sweep(data,2,as.numeric(samples_data$Weight),FUN='/')
data = log2(data+1) ## log2


## Get colored side bars
sample_factoring = colnames(data)
for (i in 1:nsamples) {
  sample_type = sample_types[i]
  replicates_want = sample_type_list[[sample_type]]
  sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix
write.table(data, file=paste(matrix.file,"_Treatment_norm.dat",sep=''), quote=F, sep='\t');



## Check we have at least 2 rows and columns
if (nrow(data) < 2) { stop("

**** Sorry, at least two rows are required for this matrix.

");}
if (ncol(data) < 2) { stop("

**** Sorry, at least two columns are required for this matrix.

");}

sample_cor = cor(data, method='pearson', use='pairwise.complete.obs')
sample_dist = dist(t(data), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
write.table(sample_cor, file=paste(matrix.file,"_Treatment_cor.dat",sep=''), quote=F, sep='\t')



## Plot heatmap
pdf(paste(matrix.file,"_Treatment_cor_matrix.pdf",sep=''))
if (is.null(hc_samples)) { RowV=NULL; ColV=NULL} else { RowV=as.dendrogram(hc_samples); ColV=RowV }
heatmap.3(sample_cor, dendrogram='both', Rowv=RowV, Colv=ColV, col = myheatcol, 
          scale='none', symm=TRUE, key=TRUE,density.info='none', trace='none', symkey=FALSE, symbreaks=F, 
          margins=c(10,10), cexCol=0.1, cexRow=0.1, cex.main=0.75, 
          main=paste("sample correlation matrix\n", matrix.file, ".norm.data") , 
          ColSideColors=sampleAnnotations, RowSideColors=t(sampleAnnotations)
)
dev.off()