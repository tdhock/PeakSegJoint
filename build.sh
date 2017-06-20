#!/bin/bash
cd ..
rm -rf PeakSegJoint-release
cp -r PeakSegJoint PeakSegJoint-release
grep -v Remotes PeakSegJoint/DESCRIPTION > PeakSegJoint-release/DESCRIPTION
rm -rf PeakSegJoint-release/inst
rm PeakSegJoint-release/data/H3K4me3.PGP.immune.4608.RData
rm PeakSegJoint-release/man/H3K4me3.PGP.immune.4608.Rd
rm PeakSegJoint-release/tests/testthat/test-bash-*.R
PKG_TGZ=$(R CMD build PeakSegJoint-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
