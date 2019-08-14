class MaskTemplate:
    def __init__(self,mask):
        self.mask = mask
    def __contains__(self,element):
        return element in self.mask
    def __str__(self):
        return str(self.mask)
    def __add__(self,other):
        return self.mask+other.mask
    def count(self,what,begin=0,end=-1):
        if end < 0:
            end = len(self.mask)+1
        return self.mask.count(what,begin,end)
    def find(self,what):
        return self.mask.find(what)

maskFull = MaskTemplate("$H$He$Li$Be$B$C$N$O$F$Ne$Na$Mg$Al$Si$P$S$Cl$Ar$K$Ca$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$Ga$Ge$As$Se$Br\
$Kr$Rb$Sr$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$In$Sn$Sb$Te$I$Xe$Cs$Ba$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$Hf$Ta$W$Re\
$Os$Ir$Pt$Au$Hg$Tl$Pb$Bi$Po$At$Rn$Fr$Ra$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$Rf$Db$Sg$Bh$Hs$Mt$Ds$Rg$Cn$Nh$Fl\
$Mc$Lv$Ts$Og$")
blocks = {
  "P": MaskTemplate("$B$C$N$O$F$Al$Si$P$S$Cl$Ga$Ge$As$Se$Br$In$Sn$Sb$Te$I$Tl$Pb$Bi$Po$At$"),
  "D": MaskTemplate("$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$"),
  "F": MaskTemplate("$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$")}
periods = {
 "2p": MaskTemplate("$B$C$N$O$F$"),
 "3p": MaskTemplate("$Al$Si$P$S$Cl$"),
 "4p": MaskTemplate("$Ga$Ge$As$Se$Br$"),
 "5p": MaskTemplate("$In$Sn$Sb$Te$I$"),
 "6p": MaskTemplate("$Tl$Pb$Bi$Po$At$"),
 "3d": MaskTemplate("$Sc$Ti$V$Cr$Mn$Fe$Co$Ni$Cu$Zn$"),
 "4d": MaskTemplate("$Y$Zr$Nb$Mo$Tc$Ru$Rh$Pd$Ag$Cd$"),
 "5d": MaskTemplate("$Hf$Ta$W$Re$Os$Ir$Pt$Au$Hg$"),
 "4f": MaskTemplate("$La$Ce$Pr$Nd$Pm$Sm$Eu$Gd$Tb$Dy$Ho$Er$Tm$Yb$Lu$"),
 "5f": MaskTemplate("$Ac$Th$Pa$U$Np$Pu$Am$Cm$Bk$Cf$Es$Fm$Md$No$Lr$")}
groups  = [None,
 MaskTemplate("$H$Li$Na$K$Rb$Cs$Fr$"),
 MaskTemplate("$Be$Mg$Ca$Sr$Ba$Ra$"),
 MaskTemplate("$Sc$Y$Lu$Lr$La$Ac$"),
 MaskTemplate("$Ti$Zr$Hf$Rf$Ce$Th$"),
 MaskTemplate("$V$Nb$Ta$Db$Pr$Pa$"),
 MaskTemplate("$Cr$Mo$W$Sg$Nd$U$"),
 MaskTemplate("$Mn$Tc$Re$Bh$Pm$Np$"),
 MaskTemplate("$Fe$Ru$Os$Hs$Sm$Pu$"),
 MaskTemplate("$Co$Rh$Ir$Mt$Eu$Am$"),
 MaskTemplate("$Ni$Pd$Pt$Ds$Gd$Cm$"),
 MaskTemplate("$Cu$Ag$Au$Rg$Tb$Bk$"),
 MaskTemplate("$Zn$Cd$Hg$Cn$Dy$Cf$"),
 MaskTemplate("$B$Al$Ga$In$Tl$Nh$Ho$Es$"),
 MaskTemplate("$C$Si$Ge$Sn$Pb$Fl$Er$Fm$"),
 MaskTemplate("$N$P$As$Sb$Bi$Mc$Tm$Md$"),
 MaskTemplate("$O$S$Se$Te$Po$Lv$Yb$No$"),
 MaskTemplate("$F$Cl$Br$I$At$Ts$"),
 MaskTemplate("$He$Ne$Ar$Kr$Xe$Rn$Og$")]
