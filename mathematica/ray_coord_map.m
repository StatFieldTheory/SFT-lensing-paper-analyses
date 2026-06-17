{"kSrc_time_scalar" -> 
  (EE^2*(2*Hconf[tt]*(PhiF[tt, x1, x2, x3] + PsiF[tt, x1, x2, x3]) + 
     2*Cos[th]*Derivative[0, 0, 0, 1][PsiF][tt, x1, x2, x3] + 
     2*Sin[th]*(Sin[ph]*Derivative[0, 0, 1, 0][PsiF][tt, x1, x2, x3] + 
       Cos[ph]*Derivative[0, 1, 0, 0][PsiF][tt, x1, x2, x3]) + 
     Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] - 
     Derivative[1, 0, 0, 0][PsiF][tt, x1, x2, x3]))/a[tt]^2, 
 "kSrc_radial_scalar" -> 
  (EE^2*(Cos[th]*(Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, x3] - 
       Derivative[0, 0, 0, 1][PsiF][tt, x1, x2, x3]) + 
     Sin[ph]*Sin[th]*(Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, x3] - 
       Derivative[0, 0, 1, 0][PsiF][tt, x1, x2, x3]) + 
     Cos[ph]*Sin[th]*(Derivative[0, 1, 0, 0][PhiF][tt, x1, x2, x3] - 
       Derivative[0, 1, 0, 0][PsiF][tt, x1, x2, x3]) - 
     2*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3]))/a[tt]^2, 
 "kSrc_trans_scalar_pole" -> 
  {-((EE^2*(Derivative[0, 1, 0, 0][PhiF][tt, x1, x2, x3] + 
       Derivative[0, 1, 0, 0][PsiF][tt, x1, x2, x3]))/a[tt]^2), 
   -((EE^2*(Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, x3] + 
       Derivative[0, 0, 1, 0][PsiF][tt, x1, x2, x3]))/a[tt]^2), 0}, 
 "dchi_dlambda_bg" -> EE/a[tt], "det_JBg" -> EE/a[tt], 
 "JBg" -> {{EE/a[tt], 0, 0}, {0, 1, 0}, {0, 0, 1}}, 
 "dchi_dz_bg" -> a[tt]/Hconf[tt], "dzdlamBg_z_only" -> 
  E0*(1 + zz)^2*Hcosm[tt], "delta_k_RHS_Gam1" -> 
  {-1/2*(EE^2*(-4*Cos[th]*Derivative[0, 0, 0, 1][PsiF][tt, x1, x2, x3] + 
       2*Cos[th]^2*Derivative[0, 0, 0, 2][BF][tt, x1, x2, x3] - 
       4*Sin[ph]*Sin[th]*Derivative[0, 0, 1, 0][PsiF][tt, x1, x2, x3] + 
       2*Sin[ph]*Sin[2*th]*Derivative[0, 0, 1, 1][BF][tt, x1, x2, x3] + 
       2*Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 2, 0][BF][tt, x1, x2, x3] + 
       2*Hconf[tt]*(Cos[th]^2*hF33[tt, x1, x2, x3] - 2*PhiF[tt, x1, x2, x3] - 
         2*Cos[th]^2*PsiF[tt, x1, x2, x3] + Cos[ph]^2*hF11[tt, x1, x2, x3]*
          Sin[th]^2 - 2*Cos[ph]^2*PsiF[tt, x1, x2, x3]*Sin[th]^2 + 
         hF22[tt, x1, x2, x3]*Sin[ph]^2*Sin[th]^2 - 2*PsiF[tt, x1, x2, x3]*
          Sin[ph]^2*Sin[th]^2 + hF12[tt, x1, x2, x3]*Sin[2*ph]*Sin[th]^2 + 
         Cos[ph]*hF13[tt, x1, x2, x3]*Sin[2*th] + hF23[tt, x1, x2, x3]*
          Sin[ph]*Sin[2*th] + 2*Cos[th]*Derivative[0, 0, 0, 1][BF][tt, x1, 
           x2, x3] + 2*Sin[ph]*Sin[th]*Derivative[0, 0, 1, 0][BF][tt, x1, x2, 
           x3] + 2*Cos[ph]*Sin[th]*Derivative[0, 1, 0, 0][BF][tt, x1, x2, 
           x3]) - 4*Cos[ph]*Sin[th]*Derivative[0, 1, 0, 0][PsiF][tt, x1, x2, 
         x3] + 2*Cos[ph]*Sin[2*th]*Derivative[0, 1, 0, 1][BF][tt, x1, x2, 
         x3] + 2*Sin[2*ph]*Sin[th]^2*Derivative[0, 1, 1, 0][BF][tt, x1, x2, 
         x3] + 2*Cos[ph]^2*Sin[th]^2*Derivative[0, 2, 0, 0][BF][tt, x1, x2, 
         x3] + Cos[ph]^2*Sin[th]^2*Derivative[1, 0, 0, 0][hF11][tt, x1, x2, 
         x3] + Sin[2*ph]*Sin[th]^2*Derivative[1, 0, 0, 0][hF12][tt, x1, x2, 
         x3] + Cos[ph]*Sin[2*th]*Derivative[1, 0, 0, 0][hF13][tt, x1, x2, 
         x3] + Sin[ph]^2*Sin[th]^2*Derivative[1, 0, 0, 0][hF22][tt, x1, x2, 
         x3] + Sin[ph]*Sin[2*th]*Derivative[1, 0, 0, 0][hF23][tt, x1, x2, 
         x3] + Cos[th]^2*Derivative[1, 0, 0, 0][hF33][tt, x1, x2, x3] - 
       2*Cos[th]^2*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] - 
       2*Cos[ph]^2*Sin[th]^2*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] - 
       2*Sin[ph]^2*Sin[th]^2*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] + 
       2*Derivative[1, 0, 0, 0][PsiF][tt, x1, x2, x3]))/a[tt]^2, 
   -1/2*(EE^2*(Cos[ph]*Sin[2*th]*Derivative[0, 0, 0, 1][hF11][tt, x1, x2, 
         x3] + Sin[ph]*Sin[2*th]*Derivative[0, 0, 0, 1][hF12][tt, x1, x2, 
         x3] + 2*Cos[th]^2*Derivative[0, 0, 0, 1][hF13][tt, x1, x2, x3] - 
       4*Cos[ph]*Cos[th]*Sin[th]*Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, 
         x3] + Sin[2*ph]*Sin[th]^2*Derivative[0, 0, 1, 0][hF11][tt, x1, x2, 
         x3] + 2*Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][hF12][tt, x1, x2, 
         x3] + Sin[ph]*Sin[2*th]*Derivative[0, 0, 1, 0][hF13][tt, x1, x2, 
         x3] - 4*Cos[ph]*Sin[ph]*Sin[th]^2*Derivative[0, 0, 1, 0][PhiF][tt, 
         x1, x2, x3] + Cos[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][hF11][tt, 
         x1, x2, x3] - Sin[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][hF22][tt, 
         x1, x2, x3] - 2*Cos[th]*Sin[ph]*Sin[th]*Derivative[0, 1, 0, 0][hF23][
         tt, x1, x2, x3] - Cos[th]^2*Derivative[0, 1, 0, 0][hF33][tt, x1, x2, 
         x3] + 2*Cos[th]^2*Derivative[0, 1, 0, 0][PhiF][tt, x1, x2, x3] - 
       2*Cos[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][PhiF][tt, x1, x2, x3] + 
       2*Sin[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][PhiF][tt, x1, x2, x3] + 
       2*Derivative[0, 1, 0, 0][PsiF][tt, x1, x2, x3] - 
       2*Cos[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF11][tt, x1, x2, x3] - 
       2*Sin[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF12][tt, x1, x2, x3] - 
       2*Cos[th]*Derivative[1, 0, 0, 0][hF13][tt, x1, x2, x3] + 
       4*Cos[ph]*Sin[th]*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] - 
       2*Derivative[1, 1, 0, 0][BF][tt, x1, x2, x3]))/a[tt]^2, 
   (EE^2*(-2*Cos[ph]*Cos[th]*Sin[th]*Derivative[0, 0, 0, 1][hF12][tt, x1, x2, 
        x3] - 2*Cos[th]*Sin[ph]*Sin[th]*Derivative[0, 0, 0, 1][hF22][tt, x1, 
        x2, x3] - 2*Cos[th]^2*Derivative[0, 0, 0, 1][hF23][tt, x1, x2, x3] + 
      4*Cos[th]*Sin[ph]*Sin[th]*Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, 
        x3] + Cos[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][hF11][tt, x1, x2, 
        x3] + Cos[ph]*Sin[2*th]*Derivative[0, 0, 1, 0][hF13][tt, x1, x2, 
        x3] - Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][hF22][tt, x1, x2, 
        x3] + Cos[th]^2*Derivative[0, 0, 1, 0][hF33][tt, x1, x2, x3] - 
      2*Cos[th]^2*Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, x3] - 
      2*Cos[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, x3] + 
      2*Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, x3] - 
      2*Derivative[0, 0, 1, 0][PsiF][tt, x1, x2, x3] - 
      2*Cos[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][hF12][tt, x1, x2, x3] - 
      2*Cos[ph]*Sin[ph]*Sin[th]^2*Derivative[0, 1, 0, 0][hF22][tt, x1, x2, 
        x3] - 2*Cos[ph]*Cos[th]*Sin[th]*Derivative[0, 1, 0, 0][hF23][tt, x1, 
        x2, x3] + 4*Cos[ph]*Sin[ph]*Sin[th]^2*Derivative[0, 1, 0, 0][PhiF][
        tt, x1, x2, x3] + 2*Cos[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF12][tt, 
        x1, x2, x3] + 2*Sin[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF22][tt, x1, 
        x2, x3] + 2*Cos[th]*Derivative[1, 0, 0, 0][hF23][tt, x1, x2, x3] - 
      4*Sin[ph]*Sin[th]*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] + 
      2*Derivative[1, 0, 1, 0][BF][tt, x1, x2, x3]))/(2*a[tt]^2), 
   -1/2*(EE^2*(-(Cos[ph]^2*Sin[th]^2*Derivative[0, 0, 0, 1][hF11][tt, x1, x2, 
          x3]) - 2*Cos[ph]*Sin[ph]*Sin[th]^2*Derivative[0, 0, 0, 1][hF12][tt, 
         x1, x2, x3] - Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 0, 1][hF22][tt, 
         x1, x2, x3] + Cos[th]^2*Derivative[0, 0, 0, 1][hF33][tt, x1, x2, 
         x3] - 2*Cos[th]^2*Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, x3] + 
       2*Cos[ph]^2*Sin[th]^2*Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, x3] + 
       2*Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 0, 1][PhiF][tt, x1, x2, x3] + 
       2*Derivative[0, 0, 0, 1][PsiF][tt, x1, x2, x3] + 
       Sin[2*ph]*Sin[th]^2*Derivative[0, 0, 1, 0][hF13][tt, x1, x2, x3] + 
       2*Sin[ph]^2*Sin[th]^2*Derivative[0, 0, 1, 0][hF23][tt, x1, x2, x3] + 
       Sin[ph]*Sin[2*th]*Derivative[0, 0, 1, 0][hF33][tt, x1, x2, x3] - 
       4*Cos[th]*Sin[ph]*Sin[th]*Derivative[0, 0, 1, 0][PhiF][tt, x1, x2, 
         x3] + 2*Cos[ph]^2*Sin[th]^2*Derivative[0, 1, 0, 0][hF13][tt, x1, x2, 
         x3] + Sin[2*ph]*Sin[th]^2*Derivative[0, 1, 0, 0][hF23][tt, x1, x2, 
         x3] + Cos[ph]*Sin[2*th]*Derivative[0, 1, 0, 0][hF33][tt, x1, x2, 
         x3] - 4*Cos[ph]*Cos[th]*Sin[th]*Derivative[0, 1, 0, 0][PhiF][tt, x1, 
         x2, x3] - 2*Cos[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF13][tt, x1, x2, 
         x3] - 2*Sin[ph]*Sin[th]*Derivative[1, 0, 0, 0][hF23][tt, x1, x2, 
         x3] - 2*Cos[th]*Derivative[1, 0, 0, 0][hF33][tt, x1, x2, x3] + 
       4*Cos[th]*Derivative[1, 0, 0, 0][PhiF][tt, x1, x2, x3] - 
       2*Derivative[1, 0, 0, 1][BF][tt, x1, x2, x3]))/a[tt]^2}}
