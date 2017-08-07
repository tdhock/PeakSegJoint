#!/bin/bash
cd ..
rm -rf PeakSegJoint-release
cp -r PeakSegJoint PeakSegJoint-release
grep -v Remotes PeakSegJoint/DESCRIPTION > PeakSegJoint-release/DESCRIPTION
rm -rf PeakSegJoint-release/inst
PKG_TGZ=$(R CMD build PeakSegJoint-release|grep building|sed 's/.*‘//'|sed 's/’.*//')
R CMD INSTALL $PKG_TGZ
R CMD check --as-cran $PKG_TGZ
