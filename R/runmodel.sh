#!/bin/bash
# NOTE (after changing shell flag in modeloutcomes.R to TRUE)
# arg1: sensitivity analysis: none, base/lo/hi dscr, cdr (higher cdr for incidence), txd (completion of ATT/TPT included)

R --slave --vanilla --args <modeloutcomes.R hi & R --slave --vanilla --args <modeloutcomes.R lo &
R --slave --vanilla --args <modeloutcomes.R cfrltfu & R --slave --vanilla --args <modeloutcomes.R whocfr &
& R --slave --vanilla --args <modeloutcomes.R systRev &
R --slave --vanilla --args <modeloutcomes.R screen & R --slave --vanilla --args <modeloutcomes.R screenINT & 
R --slave --vanilla --args <modeloutcomes.R pooled_eff &
R --slave --vanilla --args <modeloutcomes.R none


