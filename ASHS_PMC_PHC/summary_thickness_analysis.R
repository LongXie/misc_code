# arguments
options(echo=TRUE)
args <- commandArgs(TRUE)
library(pROC)

# Globals
data.demog_fn=args[1];
#data.root='~/ASHS_PHC/thickness_newlabel';
#data.demog_dir=sprintf('%s/analysis_input', data.root);
#data.demog_fn=sprintf('%s/demog.csv', data.demog_dir);
#data.demog_fn='/home/longxie/ASHS_PHC/thickness_newlabel/analysis_input/restatisticanalysis/demog.csv';


data.summary_thick_fn=args[2];
#data.summary_thick_fn=sprintf('%s/exp/exp506/MultiTemps/analysis/summary_thick.csv', data.root);
#data.summary_thick_fn=sprintf('%s/exp/exp180/analysis/summary_thick.csv', data.root);
#data.summary_thick_fn='/home/longxie/ASHS_PHC/thickness_newlabel/analysis_input/restatisticanalysis/SingleGroup.csv';
#data.summary_thick_fn='/data/picsl/longxie/ASHS/thickness_newlabel/exp/summary_thick.csv';
#data.summary_thick_fn=sprintf('%s/all_summary_measures/all_summary_measures_withFS.csv', data.root);

data.out_dir=args[3];
#data.out_dir='/home/longxie/ASHS_PHC/thickness_newlabel/analysis_input/tmp.csv'

# read data
demog=read.csv(data.demog_fn);
thk_data=read.csv(data.summary_thick_fn);

# Residualize for age&icv
num_thick=ncol(thk_data);
print(num_thick)
z.raw=thk_data[,2:num_thick]
z.res=lm(as.matrix(z.raw) ~ demog$Age + demog$ICV)$residuals;
#z.res=lm(as.matrix(z.raw) ~ demog$Age )$residuals;

# Compute t-stat
stats=matrix(ncol=num_thick-1, nrow=4)
colnames(stats) <- colnames(z.raw)
rownames(stats) <- c('T-statistic', 'P-value', 'AUC', 'AUC95%R')

for(i in 1:ncol(z.raw))
{
  stats[1,i] = coef(summary(lm(z.raw[,i]~Age+ICV+Group,data=demog)))['Groupmci',3]
  stats[2,i] = coef(summary(lm(z.raw[,i]~Age+ICV+Group,data=demog)))['Groupmci',4]
  roc_tmp = roc(demog$Group ~ z.res[,i], direction=">");
  stats[3,i] = roc_tmp$auc;
  tmp = ci.auc(roc_tmp);
  stats[4,i] = tmp[3] - tmp[2];
  
  
  
}

# save output
stats
write.csv(stats, data.out_dir)
