m=z.res
# anterior ba35
leftaba35_multi_roc=roc(demog$Group ~ m[,5], direction=">")
leftaba35_unified_roc=roc(demog$Group ~ m[,17], direction=">")
leftaba35_single_roc=roc(demog$Group ~ m[,29], direction=">")
leftaba35_fs_roc=roc(demog$Group ~ m[,67], direction=">")

# ba35
leftba35_multi_roc=roc(demog$Group ~ m[,2], direction=">")
leftba35_unified_roc=roc(demog$Group ~ m[,14], direction=">")
leftba35_single_roc=roc(demog$Group ~ m[,26], direction=">")

# erc
lefterc_multi_roc=roc(demog$Group ~ m[,1], direction=">")
lefterc_unified_roc=roc(demog$Group ~ m[,13], direction=">")
lefterc_single_roc=roc(demog$Group ~ m[,25], direction=">")
lefterc_fs_roc=roc(demog$Group ~ m[,65], direction=">")

# tests
roc.test(leftaba35_fs_roc,leftaba35_single_roc, alternative = c("greater"))
roc.test(leftaba35_fs_roc,leftaba35_multi_roc, alternative = c("greater"))
roc.test(leftaba35_fs_roc,leftaba35_unified_roc, alternative = c("greater"))


roc.test(leftaba35_fs_roc,leftba35_single_roc, alternative = c("greater"))
roc.test(leftaba35_fs_roc,leftba35_multi_roc, alternative = c("greater"))
roc.test(leftaba35_fs_roc,leftba35_unified_roc, alternative = c("greater"))

roc.test(lefterc_fs_roc,lefterc_single_roc, alternative = c("greater"))
roc.test(lefterc_fs_roc,lefterc_multi_roc, alternative = c("greater"))
roc.test(lefterc_fs_roc,lefterc_unified_roc, alternative = c("greater"))

roc.test(lefterc_fs_roc,leftaba35_single_roc, alternative = c("greater"))
roc.test(lefterc_fs_roc,leftaba35_multi_roc, alternative = c("greater"))
roc.test(lefterc_fs_roc,leftaba35_unified_roc, alternative = c("greater"))

## confidence interval
ci.auc(leftaba35_single_roc)
