#
#
# MASKS
#
#
maskFull = "$H$He$Li$Be$B$C$N$O$F$Ne$Na$Mg$Al$Si$P$S$Cl$Ar$K$Ca$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$Ga$Ge$As$Se$Br$Kr$Rb$Sr$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$In$Sn$Sb$Te$I$Xe$Cs$Ba$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$Tl$Pb$Bi$Po$At$Rn$Fr$Ra$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$Rf$Db$Sg$Bh$Hs$Mt$Ds$Rg$Cn$Nh$Fl$Mc$Lv$Ts$Og$"

mask2p = "$B$C$N$O$F$"
mask3p = "$Al$Si$P$S$Cl$"
mask4p = "$Ga$Ge$As$Se$Br$"
mask5p = "$In$Sn$Sb$Te$I$"
mask6p = "$Tl$Pb$Bi$Po$At$"
maskP  = mask2p+mask3p+mask4p+mask5p+mask6p

mask3d = "$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$"
mask4d = "$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$"
mask5d = "$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$"
maskD  = mask3d+mask4d+mask5d

mask4f = "$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$"
mask5f = "$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$"
maskF  = mask4f+mask5f

maskG1   = "$H$Li$Na$K$Rb$Cs$Fr$"
maskG2   = "$Be$Mg$Ca$Sr$Ba$Ra$"
maskG3   = "$Sc$Y$Lu$Lr$La$Ac$"
maskG4   = "$Ti$Zr$Hf$Rf$Ce$Th$"
maskG5   = "$V$Nb$Ta$Db$Pr$Pa$"
maskG6   = "$Cr$Mo$W$Sg$Nd$U$"
maskG7   = "$Mn$Tc$Re$Bh$Pm$Np$"
maskG8   = "$Fe$Ru$Os$Hs$Sm$Pu$"
maskG9   = "$Co$Rh$Ir$Mt$Eu$Am$"
maskG10 = "$Ni$Pd$Pt$Ds$Gd$Cm$"
maskG11 = "$Cu$Ag$Au$Rg$Tb$Bk$"
maskG12 = "$Zn$Cd$Hg$Cn$Dy$Cf$"
maskG13 = "$B$Al$Ga$In$Tl$Nh$Ho$Es$"
maskG14 = "$C$Si$Ge$Sn$Pb$Fl$Er$Fm$"
maskG15 = "$N$P$As$Sb$Bi$Mc$Tm$Md$"
maskG16 = "$O$S$Se$Te$Po$Lv$Yb$No$"
maskG17 = "$F$Cl$Br$I$At$Ts$"
maskG18 = "$He$Ne$Ar$Kr$Xe$Rn$Og$"

blockNames  = ["P","D","F"]
periodNames = ["2p", "3p", "4p", "5p", "6p", "3d", "4d", "5d", "4f", "5f"]
periods = [mask2p, mask3p, mask4p, mask5p, mask6p, maskP,
           mask3d, mask4d, mask5d, maskD,
           mask4f, mask5f, maskF]
groups  = [maskG1,  maskG2,  maskG3,  maskG4,  maskG5,  maskG6,
           maskG7,  maskG8,  maskG9,  maskG10, maskG11, maskG12,
           maskG13, maskG14, maskG15, maskG16, maskG17, maskG18]
