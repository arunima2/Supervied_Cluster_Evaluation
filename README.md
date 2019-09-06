# Supervised_Cluster_Evaluation
Allows to view the underlying structure of the data with respect to the ground truth labels

## Usage
Rscript Assess_Clusters_with_Labels.R #File_with_data# #File_with_labels# #T/F- whether file_with_data has headers# #delimiter of file_with_data# #number of clusters expected in the data# #output_folder-previously created#
### Example
Rscript Assess_Clusters_with_Labels.R test.csv test_labels.txt F , 4 output/

## Output
### Unsupervised clustering of data, colored by known labels
Unrooted_Dendogram.pdf
Fan_Dendogram.pdf
Cluster_Histograms_Kmeans.pdf
Cluster_Histograms_Hclust.pdf

### Dimensionality reduction, colored by known labels
PCA.pdf
TSNE.pdf

### PAM/Nearest Shrunken Centroids analysis
(http://statweb.stanford.edu/~tibs/PAM/)
PAM_CV.pdf
PAM_Confusion_Matrices.pdf (Most important to see misclassification, row is True Label and column is Predicted Label)
PAM_CV_probs.pdf
PAM_Centroid_Threshold_X.pdf (50 percentile threshold)
PAM_Example_Features_Threshold_X.pdf (75 percentile threshold)

## Requirements
install.packages(stringr)
install.packages(cluster)
install.packages(fpc)
install.packages(NMF)
install.packages(ggfortify)
install.packages(Rtsne)
install.packages(ape)
install.packages(pamr)
install.packages(corrplot)

### Note!
If the input file is very large (Num of features>>>>Num of data points) consider doing some dimensionality reduction or feature selection prior to analyzing
output_folder should have a "/" at the end for successful file creation

### Example data from
Rakhlin, A., Shvets, A., Iglovikov, V., Kalinin, A.: Deep Convolutional Neural Networks for Breast Cancer Histology Image Analysis. arXiv:1802.00752 [cs.CV]
