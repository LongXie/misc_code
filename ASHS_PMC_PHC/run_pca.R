# arguments
options(echo=TRUE)
args <- commandArgs(TRUE)
require(ggplot2);

# read demon info
data.demong_fn=args[1]
#data.demong_fn='/home/longxie/ASHS_PHC/thickness_newlabel/exp/exp800/gshoot/shape_avg/iter_4/pca/IDs.csv'
#data.demong_fn='/data/jet/longxie/ASHS_PHC/thickness_newlabel/exp/exp506/MultiTemps//group1/VTGeoShoot/shavg/pca/left/IDs_left_1.csv'
X<-read.csv(data.demong_fn,header = T, na.strings = "");

# read momenta
data.momenta_fn=args[2];
#data.momenta_fn='/home/longxie/ASHS_PHC/thickness_newlabel/exp/exp800/gshoot/shape_avg/iter_4/pca/initial_momenta.txt'
#data.momenta_fn='/data/jet/longxie/ASHS_PHC/thickness_newlabel/exp/exp506/MultiTemps//group1/VTGeoShoot/shavg/pca/left/initial_momenta_left_1.txt'
X.mom.raw<-read.table(data.momenta_fn, header = FALSE);
#names(X.mom.raw)<-c("ID","SIDE",paste('M',1:(dim(X.mom.raw)[2]-2),sep=''));

# Read the mesh points
data.points_fn=args[3];
#data.points_fn='/home/longxie/ASHS_PHC/thickness_newlabel/exp/exp800/gshoot/shape_avg/iter_4/pca/template_points.txt'
#data.points_fn='/data/jet/longxie/ASHS_PHC/thickness_newlabel/exp/exp506/MultiTemps//group1/VTGeoShoot/shavg/pca/left/template_points_left_1.txt'
X.q <- read.table(data.points_fn, header = FALSE);

# read label probabilities
data.labelprobs_fn=args[4];
#data.labelprobs_fn='/home/longxie/ASHS_PHC/thickness_newlabel/exp/exp800/gshoot/shape_avg/iter_4/pca/label_probs.txt'
#data.labelprobs_fn='/data/jet/longxie/ASHS_PHC/thickness_newlabel/exp/exp506/MultiTemps//group1/VTGeoShoot/shavg/pca/left/label_probs_left_1.txt'
X.labelprobs<-read.table(data.labelprobs_fn, header = FALSE);

# Compute the kernel matrix
dmat = as.matrix(dist(X.q));
sigma.gs = 2.0;
G.dmat = exp(-dmat^2 / (2 * sigma.gs^2))

# Function to map velocities to momenta
momentum_to_velocity<-function (p)
{
  # Split into columns for x, y, z momentum
  pmat = matrix(as.double(p), ncol=3, byrow=TRUE);
  
  # Multiply by the G matrix
  vmat = cbind(G.dmat %*% as.matrix(pmat[,1]), 
               G.dmat %*% as.matrix(pmat[,2]), 
               G.dmat %*% as.matrix(pmat[,3]));
  
  # Flatten into a vector
  matrix(t(vmat), nrow=1)
}

# Function to map momenta to velocities
velocity_to_momentum<-function (v)
{
  # Split into columns for x, y, z momentum
  vmat = matrix(as.double(v), ncol=3, byrow=TRUE);
  
  # Multiply by the G matrix
  pmat = cbind(solve(G.dmat,as.matrix(vmat[,1])), 
               solve(G.dmat,as.matrix(vmat[,2])), 
               solve(G.dmat,as.matrix(vmat[,3])));
  
  # Flatten into a vector
  matrix(t(pmat), nrow=1)
}

# This maps momentum to velocity
X.vel.raw = cbind(X.mom.raw[,1],
                  t(apply(X.mom.raw[,2:dim(X.mom.raw)[2]],1,momentum_to_velocity)));

# Add the other important columns
#X.vel = merge(X[,c('ID')],X.vel.raw)
#i.vel = 2:dim(X.vel)[2]

pca.vel = prcomp(X.vel.raw[,2:dim(X.mom.raw)[2]], retx = T, center = F, scale. = F);
X.pca.vel = cbind(X, pca.vel$x);

pca.var = data.frame(Mode=1:length(pca.vel$sdev),
                     CumVar = cumsum(pca.vel$sdev^2) / sum(pca.vel$sdev^2));
#ggplot(data=pca.var, aes(x=Mode,y=CumVar)) +
#  geom_point() + theme_classic() + 
#  geom_hline(aes(yintercept=0.95)) +
#  labs(title="Variance explained by PCA modes")

#ggplot(data=X.pca.vel, aes(x=PC1/pca.vel$sdev[1], y=PC2/pca.vel$sdev[2])) +
#  geom_point() + theme_classic() +
#  labs(x="PCA Mode 1", y = "PCA Mode 2")


write.vtk<-function(p, file)
{
  sink(file)
  
  # Write VTK header
  cat("# vtk DataFile Version 4.0\n")
  cat("vtk output\n")
  cat("ASCII\n")
  cat("DATASET POLYDATA\n")
  cat(paste("POINTS", dim(p)[1], "float\n"))
  write.table(as.matrix(X.q),row.names=F,col.names=F);
  cat(paste("\nPOINT_DATA",dim(p)[1],"\n"));
  cat("FIELD FieldData 1\n");
  cat(paste("InitialMomentum 3",dim(p)[1],"double\n"));
  write.table(p,row.names=F,col.names=F);
  sink(NULL)
}

write.vtkmultiple<-function(p, G, LP, pca, file)
{
  sink(file)
  
  # Write VTK header
  cat("# vtk DataFile Version 4.0\n")
  cat("vtk output\n")
  cat("ASCII\n")
  cat("DATASET POLYDATA\n")
  cat(paste("POINTS", dim(p)[1], "float\n"))
  write.table(as.matrix(X.q),row.names=F,col.names=F);
  cat(paste("\nPOINT_DATA",dim(p)[1],"\n"));
  
  
  cat("FIELD FieldData 7\n");
  cat(paste("InitialMomentum 3",dim(p)[1],"double\n"));
  write.table(p,row.names=F,col.names=F);
  
  # write label probabilities
  cat(paste("LabelProbs ", dim(LP)[2], dim(LP)[1],"double\n"));
  write.table(LP,row.names=F,col.names=F);
  
  # write G matrix (for now)
  G.matrix = matrix(G, ncol=dim(p)[1], byrow = TRUE);
  cat(paste("G ", dim(G.matrix)[2], dim(G.matrix)[1],"double\n"));
  write.table(G.matrix,row.names=F,col.names=F);
  
  # write loading
  Erep = t(matrix(rep(pca$sdev, dim(p)[1]), ncol = dim(p)[1]));
  cat(paste("PCALoading ", dim(Erep)[2], dim(Erep)[1],"double\n"));
  write.table(Erep,row.names=F,col.names=F);
  
  # write rotation
  TransX = matrix(, nrow = dim(p)[1], ncol = dim(Erep)[2]);
  TransY = matrix(, nrow = dim(p)[1], ncol = dim(Erep)[2]);
  TransZ = matrix(, nrow = dim(p)[1], ncol = dim(Erep)[2]);
  for (i in 1:dim(pca$rotation)[2])
  {
    Trans.data = matrix(as.double(pca$rotation[,i]), ncol=3, byrow=TRUE);
    TransX[,i] = Trans.data[,1];
    TransY[,i] = Trans.data[,2];
    TransZ[,i] = Trans.data[,3];
  }
  cat(paste("PCATransX ", dim(TransX)[2], dim(TransX)[1],"double\n"));
  write.table(TransX,row.names=F,col.names=F);
  cat(paste("PCATransY ", dim(TransY)[2], dim(TransY)[1],"double\n"));
  write.table(TransY,row.names=F,col.names=F);
  cat(paste("PCATransZ ", dim(TransZ)[2], dim(TransZ)[1],"double\n"));
  write.table(TransZ,row.names=F,col.names=F);
  
  sink(NULL)
}


mom.mode1 = matrix(velocity_to_momentum(pca.vel$rotation[,1] * pca.vel$sdev[1]), 
                   ncol=3, byrow = TRUE);

mom.mode2 = matrix(velocity_to_momentum(pca.vel$rotation[,2] * pca.vel$sdev[2]), 
                   ncol=3, byrow = TRUE);

mom.mode3 = matrix(velocity_to_momentum(pca.vel$rotation[,3] * pca.vel$sdev[3]), 
                   ncol=3, byrow = TRUE);

mom.mode4 = matrix(velocity_to_momentum(pca.vel$rotation[,4] * pca.vel$sdev[4]), 
                   ncol=3, byrow = TRUE);

mom.mode5 = matrix(velocity_to_momentum(pca.vel$rotation[,5] * pca.vel$sdev[5]), 
                   ncol=3, byrow = TRUE);

data.OUTROOT = args[5]
#data.OUTROOT = '/data/jet/longxie/ASHS_PHC/thickness_newlabel/exp/exp506/MultiTemps//group1/VTGeoShoot/shavg/pca/left'
data.mode1_pos_fn = file.path(data.OUTROOT, 'mode1_vector_pos.vtk')
data.mode1_neg_fn = file.path(data.OUTROOT, 'mode1_vector_neg.vtk')
write.vtk(mom.mode1, data.mode1_pos_fn)
write.vtk(-mom.mode1, data.mode1_neg_fn)

data.mode2_pos_fn = file.path(data.OUTROOT, 'mode2_vector_pos.vtk')
data.mode2_neg_fn = file.path(data.OUTROOT, 'mode2_vector_neg.vtk')
write.vtk(mom.mode2, data.mode2_pos_fn)
write.vtk(-mom.mode2, data.mode2_neg_fn)

data.mode3_pos_fn = file.path(data.OUTROOT, 'mode3_vector_pos.vtk')
data.mode3_neg_fn = file.path(data.OUTROOT, 'mode3_vector_neg.vtk')
write.vtk(mom.mode3, data.mode3_pos_fn)
write.vtk(-mom.mode3, data.mode3_neg_fn)

data.mode4_pos_fn = file.path(data.OUTROOT, 'mode4_vector_pos.vtk')
data.mode4_neg_fn = file.path(data.OUTROOT, 'mode4_vector_neg.vtk')
write.vtk(mom.mode4, data.mode4_pos_fn)
write.vtk(-mom.mode4, data.mode4_neg_fn)

data.mode5_pos_fn = file.path(data.OUTROOT, 'mode5_vector_pos.vtk')
data.mode5_neg_fn = file.path(data.OUTROOT, 'mode5_vector_neg.vtk')
write.vtk(mom.mode5, data.mode5_pos_fn)
write.vtk(-mom.mode5, data.mode5_neg_fn)

# Save other information for template fitting
# compute mean initial momentum
X.mean = matrix(apply(X.mom.raw[,2:dim(X.mom.raw)[2]], 2, mean), ncol = 3, byrow = TRUE);
LabelProb = t(matrix(t(X.labelprobs[,2:dim(X.labelprobs)[2]]), ncol = dim(X.labelprobs)[2]-1, byrow = TRUE));
data.template_fn = file.path(data.OUTROOT, 'template_allinfo.vtk')
write.vtkmultiple(X.mean, G.dmat, LabelProb, pca.vel, data.template_fn)
