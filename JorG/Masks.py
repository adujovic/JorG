from abc import ABC

class MaskTemplate:
    def __init__(self,mask):
        self.mask = mask
    def __contains__(self,element):
        return element in self.mask
    def __add__(self,other):
        return self.mask+other.mask
#
#
# MASKS
#
#


maskFull = MaskTemplate("$H$He$Li$Be$B$C$N$O$F$Ne$Na$Mg$Al$Si$P$S$Cl$Ar$K$Ca$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$Ga$Ge$As$Se$Br$Kr$Rb$Sr$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$In$Sn$Sb$Te$I$Xe$Cs$Ba$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$Tl$Pb$Bi$Po$At$Rn$Fr$Ra$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$Rf$Db$Sg$Bh$Hs$Mt$Ds$Rg$Cn$Nh$Fl$Mc$Lv$Ts$Og$")

mask2p = MaskTemplate("$B$C$N$O$F$")
mask3p = MaskTemplate("$Al$Si$P$S$Cl$")
mask4p = MaskTemplate("$Ga$Ge$As$Se$Br$")
mask5p = MaskTemplate("$In$Sn$Sb$Te$I$")
mask6p = MaskTemplate("$Tl$Pb$Bi$Po$At$")
maskP  = MaskTemplate("$B$C$N$O$F$Al$Si$P$S$Cl$Ga$Ge$As$Se$Br$In$Sn$Sb$Te$I$Tl$Pb$Bi$Po$At$")

mask3d = MaskTemplate("$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$")
mask4d = MaskTemplate("$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$")
mask5d = MaskTemplate("$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$")
maskD  = MaskTemplate("$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$")

mask4f = MaskTemplate("$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$")
mask5f = MaskTemplate("$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$")
maskF  = MaskTemplate("$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$")

maskG1  = MaskTemplate("$H$Li$Na$K$Rb$Cs$Fr$")
maskG2  = MaskTemplate("$Be$Mg$Ca$Sr$Ba$Ra$")
maskG3  = MaskTemplate("$Sc$Y$Lu$Lr$La$Ac$")
maskG4  = MaskTemplate("$Ti$Zr$Hf$Rf$Ce$Th$")
maskG5  = MaskTemplate("$V$Nb$Ta$Db$Pr$Pa$")
maskG6  = MaskTemplate("$Cr$Mo$W$Sg$Nd$U$")
maskG7  = MaskTemplate("$Mn$Tc$Re$Bh$Pm$Np$")
maskG8  = MaskTemplate("$Fe$Ru$Os$Hs$Sm$Pu$")
maskG9  = MaskTemplate("$Co$Rh$Ir$Mt$Eu$Am$")
maskG10 = MaskTemplate("$Ni$Pd$Pt$Ds$Gd$Cm$")
maskG11 = MaskTemplate("$Cu$Ag$Au$Rg$Tb$Bk$")
maskG12 = MaskTemplate("$Zn$Cd$Hg$Cn$Dy$Cf$")
maskG13 = MaskTemplate("$B$Al$Ga$In$Tl$Nh$Ho$Es$")
maskG14 = MaskTemplate("$C$Si$Ge$Sn$Pb$Fl$Er$Fm$")
maskG15 = MaskTemplate("$N$P$As$Sb$Bi$Mc$Tm$Md$")
maskG16 = MaskTemplate("$O$S$Se$Te$Po$Lv$Yb$No$")
maskG17 = MaskTemplate("$F$Cl$Br$I$At$Ts$")
maskG18 = MaskTemplate("$He$Ne$Ar$Kr$Xe$Rn$Og$")

blockNames  = ["P","D","F"]
periodNames = ["2p", "3p", "4p", "5p", "6p", "3d", "4d", "5d", "4f", "5f"]
periods = [mask2p, mask3p, mask4p, mask5p, mask6p, maskP,
           mask3d, mask4d, mask5d, maskD,
           mask4f, mask5f, maskF]
groups  = [maskG1,  maskG2,  maskG3,  maskG4,  maskG5,  maskG6,
           maskG7,  maskG8,  maskG9,  maskG10, maskG11, maskG12,
           maskG13, maskG14, maskG15, maskG16, maskG17, maskG18]
