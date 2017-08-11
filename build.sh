#!/bin/bash
cd ..
rm -rf PeakSegJoint-release
cp -r PeakSegJoint PeakSegJoint-release
##grep -v Remotes PeakSegJoint/DESCRIPTION > PeakSegJoint-release/DESCRIPTION
#                                f seconds
# 1:      test-functions-trivial.R   0.077
# 2:    test-functions-converted.R   0.423
# 3:     test-functions-features.R   0.526
# 4:   test-functions-likelihood.R   2.928
# 5:       test-functions-binSum.R  29.041
# 6:     test-functions-no-zeros.R  38.127
# 7: test-functions-PeakSegJoint.R  56.224
# 8:    test-functions-real-data.R 674.881
rm PeakSegJoint-release/test-functions-real-data.R
PKG_TGZ=$(R CMD build PeakSegJoint-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
