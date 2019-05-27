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
#
#
# DICTIONARIES
#
#
periodicTableElement = {  0:"H",    1:"He",   2:"Li",   3:"Be",   4:"B",    5:"C",    6:"N",    7:"O",    8:"F",    9:"Ne",
                         10:"Na",  11:"Mg",  12:"Al",  13:"Si",  14:"P",   15:"S",   16:"Cl",  17:"Ar",  18:"K",   19:"Ca",
                         20:"Sc",  21:"Ti",  22:"V",   23:"Cr",  24:"Mn",  25:"Fe",  26:"Co",  27:"Ni",  28:"Cu",  29:"Zn",
                         30:"Ga",  31:"Ge",  32:"As",  33:"Se",  34:"Br",  35:"Kr",  36:"Rb",  37:"Sr",  38:"Y",   39:"Zr",
                         40:"Nb",  41:"Mo",  42:"Tc",  43:"Ru",  44:"Rh",  45:"Pd",  46:"Ag",  47:"Cd",  48:"In",  49:"Sn",
                         50:"Sb",  51:"Te",  52:"I",   53:"Xe",  54:"Cs",  55:"Ba",  56:"La",  57:"Ce",  58:"Pr",  59:"Nd",
                         60:"Pm",  61:"Sm",  62:"Eu",  63:"Gd",  64:"Tb",  65:"Dy",  66:"Ho",  67:"Er",  68:"Tm",  69:"Yb",
                         70:"Lu",  71:"Hf",  72:"Ta",  73:"W",   74:"Re",  75:"Os",  76:"Ir",  77:"Pt",  78:"Au",  79:"Hg",
                         80:"Tl",  81:"Pb",  82:"Bi",  83:"Po",  84:"At",  85:"Rn",  86:"Fr",  87:"Ra",  88:"Ac",  89:"Th",
                         90:"Pa",  91:"U",   92:"Np",  93:"Pu",  94:"Am",  95:"Cm",  96:"Bk",  97:"Cf",  98:"Es",  99:"Fm",
                        100:"Md", 101:"No", 102:"Lr", 103:"Rf", 104:"Db", 105:"Sg", 106:"Bh", 107:"Hs", 108:"Mt", 109:"Ds",
                        110:"Rg", 111:"Cn", 112:"Nh", 113:"Fl", 114:"Mc", 115:"Lv", 116:"Ts", 117:"Og" }
periodicTableNumber = { "H"  :1  , "He" :2  , "Li" :3  , "Be" :4  , "B"  :5  , "C"  :6  , "N"  :7  , "O"  :8  , "F"  :9  ,
             "Ne" :10 , "Na" :11 , "Mg" :12 , "Al" :13 , "Si" :14 , "P"  :15 , "S"  :16 , "Cl" :17 , "Ar" :18 , "K"  :19 ,
             "Ca" :20 , "Sc" :21 , "Ti" :22 , "V"  :23 , "Cr" :24 , "Mn" :25 , "Fe" :26 , "Co" :27 , "Ni" :28 , "Cu" :29 ,
             "Zn" :30 , "Ga" :31 , "Ge" :32 , "As" :33 , "Se" :34 , "Br" :35 , "Kr" :36 , "Rb" :37 , "Sr" :38 , "Y"  :39 ,
             "Zr" :40 , "Nb" :41 , "Mo" :42 , "Tc" :43 , "Ru" :44 , "Rh" :45 , "Pd" :46 , "Ag" :47 , "Cd" :48 , "In" :49 ,
             "Sn" :50 , "Sb" :51 , "Te" :52 , "I"  :53 , "Xe" :54 , "Cs" :55 , "Ba" :56 , "La" :57 , "Ce" :58 , "Pr" :59 ,
             "Nd" :60 , "Pm" :61 , "Sm" :62 , "Eu" :63 , "Gd" :64 , "Tb" :65 , "Dy" :66 , "Ho" :67 , "Er" :68 , "Tm" :69 ,
             "Yb" :70 , "Lu" :71 , "Hf" :72 , "Ta" :73 , "W"  :74 , "Re" :75 , "Os" :76 , "Ir" :77 , "Pt" :78 , "Au" :79 , 
             "Hg" :80 , "Tl" :81 , "Pb" :82 , "Bi" :83 , "Po" :84 , "At" :85 , "Rn" :86 , "Fr" :87 , "Ra" :88 , "Ac" :89 ,
             "Th" :90 , "Pa" :91 , "U"  :92 , "Np" :93 , "Pu" :94 , "Am" :95 , "Cm" :96 , "Bk" :97 , "Cf" :98 , "Es" :99 ,
             "Fm" :100, "Md" :101, "No" :102, "Lr" :103, "Rf" :104, "Db" :105, "Sg" :106, "Bh" :107, "Hs" :108, "Mt" :109,
             "Ds" :110, "Rg" :111, "Cn" :112, "Nh" :113, "Fl" :114, "Mc" :115, "Lv" :116, "Ts" :117, "Og" :118 }
elementMagneticMoment = {"H"  : 1.00, "He" : 0.00, "Li" : 0.00, "Be" : 0.00, "B"  : 0.00, "C"  : 0.00, "N"  : 0.00, "O"  : 2.00, "F"  : 0.00,
            "Ne" : 0.00, "Na" : 0.00, "Mg" : 0.00, "Al" : 0.00, "Si" : 0.00, "P"  : 0.00, "S"  : 0.00, "Cl" : 0.00, "Ar" : 0.00, "K"  : 0.00,
            "Ca" : 0.00, "Sc" : 0.68, "Ti" : 0.00, "V"  : 0.00, "Cr" : 4.80, "Mn" : 4.50, "Fe" : 2.50, "Co" : 1.90, "Ni" : 0.73, "Cu" : 0.00,
            "Zn" : 0.00, "Ga" : 0.00, "Ge" : 0.00, "As" : 0.00, "Se" : 0.00, "Br" : 0.00, "Kr" : 0.00, "Rb" : 0.00, "Sr" : 0.00, "Y"  : 0.60,
            "Zr" : 0.00, "Nb" : 0.00, "Mo" : 0.00, "Tc" : 0.00, "Ru" : 0.00, "Rh" : 0.00, "Pd" : 1.00, "Ag" : 0.00, "Cd" : 0.00, "In" : 0.00,
            "Sn" : 0.00, "Sb" : 0.00, "Te" : 0.00, "I"  : 0.00, "Xe" : 0.00, "Cs" : 0.00, "Ba" : 0.00, "La" : 0.58, "Ce" : 0.88, "Pr" : 2.88,
            "Nd" : 4.10, "Pm" : 5.34, "Sm" : 0.40, "Eu" : 7.86, "Gd" : 7.68, "Tb" : 7.28, "Dy" : 7.00, "Ho" : 4.10, "Er" : 4.12, "Tm" : 1.90,
            "Yb" : 0.00, "Lu" : 0.00, "Hf" : 0.00, "Ta" : 0.00, "W"  : 0.00, "Re" : 0.00, "Os" : 0.00, "Ir" : 0.00, "Pt" : 1.00, "Au" : 0.00,
            "Hg" : 0.00, "Tl" : 0.00, "Pb" : 0.00, "Bi" : 0.00, "Po" : 0.00, "At" : 0.00, "Rn" : 0.00, "Fr" : 0.00, "Ra" : 0.00, "Ac" : 0.00,
            "Th" : 0.00, "Pa" : 0.00, "U"  : 1.00, "Np" : 1.00, "Pu" : 1.00, "Am" : 0.00, "Cm" : 0.00, "Bk" : 0.00, "Cf" : 0.00, "Es" : 0.00,
            "Fm" : 0.00, "Md" : 0.00, "No" : 0.00, "Lr" : 0.00, "Rf" : 0.00, "Db" : 0.00, "Sg" : 0.00, "Bh" : 0.00, "Hs" : 0.00, "Mt" : 0.00,
            "Ds" : 0.00, "Rg" : 0.00, "Cn" : 0.00, "Nh" : 0.00, "Fl" : 0.00, "Mc" : 0.00, "Lv" : 0.00, "Ts" : 0.00, "Og" : 0.0 } 
