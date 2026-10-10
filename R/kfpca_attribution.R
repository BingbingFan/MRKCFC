# KFPCA-derived internal runtime
#
# The following internal source files are derived from KFPCA version 2.0:
# basicfun.R, kernfun.R, GCVbw1D.R, GCVbw2D.R, GetGCVbw1D.R,
# GetGCVbw2D.R, GetMatList.R, MeanEst.R, MeanObsGrid.R, FPCscoreLSE.R,
# KFPCA.R, and predict.KFPCA.R.
#
# Upstream KFPCA authors: Rou Zhong and Jingxiao Zhang.
# Upstream source: https://github.com/BingbingFan/KFPCA
# Upstream license: GPL (>= 3).
#
# MRKCFC integration and modifications by Bingbing Fan:
# - incorporated the runtime closure into the MRKCFC namespace;
# - replaced direct generalized-inverse score calculations in FPCscoreLSE()
#   and predict.KFPCA() with MRKCFC's solve()-first covariance solver;
# - retained these functions as internal implementation details of MRKCFC.
