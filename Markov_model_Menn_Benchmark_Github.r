
####################################################################################################
# Base R implementation of the original COPD Markov model by Petra Menn, translated from Excel.      #
# R implementation and optimization by Tobias Niedermaier. Original model design, parameterization,  #
# and implementation by Petra Menn. Supervisor: Rolf Holle.                                          #
#                                                                                                  #
# This benchmark model comprises GOLD stages 1-4, lung volume reduction surgery (LVR), lung          #
# transplantation (LTx), and death as absorbing state. Compared with the simplified GitHub model,    #
# it additionally retains mild acute exacerbations, the LVR state, and the LTx state. The goal of    #
# this script is to reproduce the original Excel model under identical assumptions before any        #
# structural simplification or parameter updating is introduced.                                    #
####################################################################################################

rm(list=ls(all=TRUE))
options(scipen=999) #disables scientific notation

# Update 2025-03-05: AE-related excess mortality differs by smoking status and stage, but it is constant and does not need to be redefined inside the probabilistic functions. It is therefore defined once outside the functions and accessed by cycle index, analogous to the simplified model. This makes the code shorter and faster; the probabilistic functions are also parallelised with foreach.

### Deterministic parameter definitions ###

#Var.name        value     meaning
prog1         = 0.01403  #probability progression stage 1
prog2         = 0.01031  #probability progression stage 2
prog3         = 0.02298  #probability progression stage 3
prog1_int     = 0.00247  #probability progression stage 1 (intervention group)
prog2_int     = 0.003    #probability progression stage 2 (intervention group)
prog3_int     = 0.00686  #probability progression stage 3 (intervention group)


# Additional original-model component not present in the simplified model:
# mild acute exacerbations are retained in GOLD stages 1-4 and affect costs and utilities.
p_mild1       = 0.07226  #probability mild exacerbation stage 1
d_mild1       = -0.0035  #change in utility mild exacerbation stage 1

p_mild2       = 0.08163  #probability mild exacerbation stage 2
d_mild2       = -0.00329 #change in utility mild exacerbation stage 2
p_mod2        = 0.18376  #probability moderate exacerbation stage 2
d_mod2        = -0.01975 #change in utility moderate exacerbation stage 2
p_sev2        = 0.03854  #probability severe exacerbation stage 2
d_sev2        = -0.03075 #change in utility severe exacerbation stage 2

p_mild3       = 0.14513  #probability mild exacerbation stage 3
d_mild3       = -0.00312 #change in utility mild exacerbation stage 3
p_mod3        = 0.26767  #probability moderate exacerbation stage 3
d_mod3        = -0.01875 #change in utility moderate exacerbation stage 3
p_sev3        = 0.03774  #probability severe exacerbation stage 3
d_sev3        = -0.02919 #change in utility severe exacerbation stage 3

p_mild4       = 0.14344  #probability mild exacerbation stage 4
d_mild4       = -0.00271 #change in utility mild exacerbation stage 4
p_mod4        = 0.28057  #probability moderate exacerbation stage 4
d_mod4        = -0.01625 #change in utility moderate exacerbation stage 4
p_sev4        = 0.04855  #probability severe exacerbation stage 4
d_sev4        = -0.03502 #change in utility severe exacerbation stage 4


# Additional original-model component not present in the simplified model:
# exacerbation probabilities and utility decrements are also specified for the post-LVR state.
p_mildL       = 0.10272  #probability mild exacerbation LVRS
d_mildL       = -0.00312 #change in utility mild exacerbation LVRS
p_modL        = 0.20587  #probability moderate exacerbation LVRS
d_modL        = -0.01875 #change in utility moderate exacerbation LVRS
p_sevL        = 0.03424  #probability severe exacerbation LVRS
d_sevL        = -0.04041 #change in utility severe exacerbation LVRS


# Additional original-model component not present in the simplified model:
# transitions from GOLD 4 to lung volume reduction surgery (LVR), lung transplantation (LTx),
# and transplantation after LVR are represented explicitly.
p_4ps         = 0.02729  #probability surgery (LVRS)
p_4pt         = 0.00294  #probability transplantation
p_pspt        = 0.00294  #probability transplantation after surgery
p_lvr_death   = 0.03     #probability to die with/after LVRS
p_death_atltx = 0.06     #probability to die with/after LTx
u_1           = 0.84     #utility stage 1
u_2           = 0.79     #utility stage 2
u_3           = 0.75     #utility stage 3
u_4           = 0.65     #utility stage 4
u_lvr         = 0.75     #utility LVRS
u_ltx         = 0.80     #utility LTx

c_1           = 102.64   #costs in stable stage 1
c_mild1       = 6.95     #costs mild ae, stage 1
c_2           = 184.69   #costs in stable stage 2
c_mild2       = 8.69     #costs mild ae, stage 2
c_2_mod       = 39.86    #costs mod ae, stage 2
c_mod2        = 39.86    #costs mod ae, stage 2 (under a different name)
c_2_sev       = 2923.94  #costs sev ae, stage 2
c_sev2        = 2923.94  #costs sev ae, stage 2 (under a different name)
c_3           = 366.57   #costs in stable stage 3
c_mild3       = 21.73    #costs mild ae, stage 3
c_mod3        = 52.90    #costs mod ae, stage 3
c_sev3        = 3109.05  #costs sev ae, stage 3
c_4           = 431.35   #costs in stable stage 4
c_mild4       = 20.35    #costs mild ae, stage 4
c_mod4        = 51.52    #costs mod ae, stage 4
c_sev4        = 3322.59  #costs sev ae, stage 4
c_int         = 595.92   #intervention costs
c_posttransp  = 8269.02  #costs in stable stage LTx
c_surgery     = 6413.88  #costs surgery
c_transp      = 65494.21 #costs transplantation 
c_atltx       = c_transp #costs transplantation (under a different name)
c_uc          = 0        #initial costs usual care
c_atlvr       = 6413.88  #costs lung volume reduction
c_ltx         = 8269.02  #Costs to live with lung transplantation (stable COPD after transplantation)

cycle_length  = 3        #[months]
disc_e        = 0.03     #discount rate
eff_int       = 0.22
eff_uc        = 0.06

## Additional notes on indirect-cost inputs:
# Excel worksheet "ind Kosten" (= indirect costs).
# Gross wages and salaries per employee per month: 3332 EUR in 2022 instead of 2221 EUR in 2004-2007; see Destatis URL.
# Labour costs per employee hour: 30.87 EUR in 2022 instead of 24.36 EUR; same source as above.

##Begin with smokers (just like in the Excel sheet), then repeat with former smokers. Excel sheet "det Modell" (!).
N <- 241 #number of cycles
dat_smo <- data.frame(matrix(nrow=N, ncol=7))
colnames(dat_smo) <- c("stage_1","stage_2","stage_3","stage_4","LVR","LTx","death")
dat_smo[] <- 0
dat_smo[1,"stage_1"] <- 100000

#Define mortality in the different cycles for the different stages for current smokers: (see Table 2.7 in Menn's dissertation)
dat_smo[,c("mort1","mort2","mort3","mort4")] <- NA
dat_smo[,c("mort1","mort2","mort3","mort4")][1:20,] <- rep(c(0.00122,0.0021,0.00497,0.04022), each=20) #Assign mortality according to disease stage (mort1 for stage_1 etc.) and cycle (higher mort. in later cycles)
dat_smo[,c("mort1","mort2","mort3","mort4")][21:40,] <- rep(c(0.0021,0.0036,0.00852,0.04022), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][41:60,] <- rep(c(0.00302,0.00518,0.01225,0.04022), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][61:80,] <- rep(c(0.0047,0.00804,0.01899,0.04022), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][81:100,] <- rep(c(0.00805,0.01377,0.03238,0.04022), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][101:120,] <- rep(c(0.01358,0.02318,0.05417,0.05417), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][121:140,] <- rep(c(0.0197,0.03354,0.07783,0.07783), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][141:160,] <- rep(c(0.034,0.05758,0.13138,0.13138), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][161:180,] <- rep(c(0.04729,0.07969,0.179,0.179), each=20)
dat_smo[,c("mort1","mort2","mort3","mort4")][181:241,] <- rep(c(0.0884,0.14672,0.31397,0.31397), each=61)

#Calculate number of individuals for dat_smo$stage_1 in each cycle (=all at the beginning, minus those who progress to stage 2, minus those who die):

dat_smo$stage_1[-1] <- (100000 * cumprod(1 - prog1 - dat_smo$mort1))[1:240]

#Stage 2, smokers:
dat_smo$stage_2 <- NA #First value stays NA by definition.
dat_smo$stage_2[2] <- dat_smo[1,"stage_1"]*prog1 #first cycle with N(stage 2) > 0 as starting point.
dat_smo$stage_2[3] <- dat_smo$stage_2[2]*(1 - prog2 - dat_smo$mort2[1] - p_sev2*0.13819) + dat_smo$stage_1[2]*prog1

ae_mort2 <- vector(mode="numeric", length=241)
ae_mort2[3:20] <- 0.13819
ae_mort2[21:40] <- 0.1369
ae_mort2[41:60] <- 0.13552
ae_mort2[61:80] <- 0.254
ae_mort2[81:100] <- 0.24967
ae_mort2[101:120] <- 0.24244
ae_mort2[121:140] <- 0.30675
ae_mort2[141:160] <- 0.28906
ae_mort2[161:180] <- 0.27198
ae_mort2[181:241] <- 0.21479

ae_mort2_fosmo <- vector(mode="numeric", length=241)
ae_mort2_fosmo[3:20] <- 0.13883
ae_mort2_fosmo[21:40] <- 0.13798
ae_mort2_fosmo[41:60] <- 0.13629
ae_mort2_fosmo[61:80] <- 0.25503
ae_mort2_fosmo[81:100] <- 0.25356
ae_mort2_fosmo[101:120] <- 0.24908
ae_mort2_fosmo[121:140] <- 0.31459
ae_mort2_fosmo[141:160] <- 0.30298
ae_mort2_fosmo[161:180] <- 0.28199
ae_mort2_fosmo[181:241] <- 0.23529

assign_stage_2_smo <- function(i){
 dat_smo$stage_2[i+1] <<- dat_smo$stage_2[i]*(1 - prog2 - dat_smo$mort2[i] - p_sev2*ae_mort2[i]) + dat_smo$stage_1[i]*prog1
}
for(i in 3:240){assign_stage_2_smo(i=i)}
dat_smo$stage_2 #OK.

#Stage 3, smokers:
dat_smo$stage_3 <- NA #First value stays NA by definition.
dat_smo$stage_3[2] <- 0 #by definition
dat_smo$stage_3[3] <- dat_smo[2,"stage_2"]*prog2 #first cycle with N(stage 3) > 0 as starting point.

ae_mort3 <- vector(mode="numeric", length=241)
ae_mort3[3:20] <- 0.1357
ae_mort3[21:40] <- 0.13261
ae_mort3[41:60] <- 0.12933
ae_mort3[61:80] <- 0.24568
ae_mort3[81:100] <- 0.23523
ae_mort3[101:120] <- 0.21762
ae_mort3[121:140] <- 0.27345
ae_mort3[141:160] <- 0.22866
ae_mort3[161:180] <- 0.18392
ae_mort3[181:241] <- 0.02336

ae_mort3_fosmo <- vector(mode="numeric", length=241)
ae_mort3_fosmo[3:20] <- 0.13721
ae_mort3_fosmo[21:40] <- 0.1352
ae_mort3_fosmo[41:60] <- 0.13116
ae_mort3_fosmo[61:80] <- 0.24813
ae_mort3_fosmo[81:100] <- 0.24462
ae_mort3_fosmo[101:120] <- 0.23379
ae_mort3_fosmo[121:140] <- 0.29281
ae_mort3_fosmo[141:160] <- 0.26404
ae_mort3_fosmo[161:180] <- 0.21031
ae_mort3_fosmo[181:241] <- 0.08281

assign_stage_3_smo <- function(i){
 dat_smo$stage_3[i+1] <<- dat_smo$stage_3[i]*(1 - prog3 - dat_smo$mort3[i] - p_sev3*ae_mort3[i]) + dat_smo$stage_2[i]*prog2
}
for(i in 3:240){assign_stage_3_smo(i=i)}
dat_smo$stage_3 #OK.

#Stage 4, smokers:
dat_smo$stage_4 <- NA #First value stays NA by definition.
dat_smo$stage_4[2:3] <- 0 #by definition
dat_smo$stage_4[4] <- dat_smo[4,"stage_3"]*prog3 #first cycle with N(stage 4) > 0 as starting point.

assign_stage_4_smo <- function(i){
 if (i >80){p_4pt <- 0} #Note the arrow/comment in the Excel sheet. The probability is set to 0 after 80 cycles, presumably because the model starts with 45-year-old COPD patients; after 80 cycles (=20 years), they are 65 years old and no longer eligible for lung transplantation.
 dat_smo$stage_4[i+1] <<- dat_smo$stage_4[i] * (1 - p_4ps - p_4pt - dat_smo$mort3[i] - p_sev4*ae_mort3[i]) + dat_smo$stage_3[i]*prog3
}
for(i in 3:240){assign_stage_4_smo(i=i)}
dat_smo$stage_4 #OK.


# Additional transition only relevant for the original model: state transition to lung volume reduction surgery (LVR).
#LVR, smokers:
dat_smo$LVR <- NA #First value stays NA by definition.
dat_smo$LVR[2:3] <- 0 #by definition
dat_smo$LVR[4] <- dat_smo[3,"stage_4"]*p_4ps*(1-p_lvr_death) #first cycle with N(LVR) > 0 as starting point.

assign_LVR_smo <- function(i){
 if (i >80){p_4pt <- p_pspt <- 0} #see comment at assign_stage_4_smo.
 dat_smo$LVR[i+1] <<- dat_smo$LVR[i] * (1 - p_pspt - dat_smo$mort3[i] - p_sevL*ae_mort3[i]) + dat_smo$stage_4[i]*p_4ps*(1-p_lvr_death)
}
for(i in 4:240){assign_LVR_smo(i=i)}
dat_smo$LVR #OK.



# Additional transition only relevant for the original model: state transition to lung transplantation (LTx).
#Calculate number of individuals for dat_smo$LTx in each cycle (=those who already are in stage LTx, plus those who progress from stage LVR, minus those who die):
 
dat_smo$LTx <- NA #First value stays NA by definition.
dat_smo$LTx[2:3] <- 0 #by definition
# dat_smo$LTx[4] <- dat_smo[3,"LTx"] * (1 - dat_smo$mort4[3]) + dat_smo[3,"stage_4"] * p_4pt * (1-p_death_atltx) + dat_smo[3,"LVR"] * p_pspt * (1 - p_death_atltx) # first cycle with N(LTx)>0 as starting point.
##                           N8      * (1 -         E8      ) +          L8          * p_4pt * (1-p_death_atltx) +         M8       * p_pspt * (1 - p_death_atltx)

dat_smo$LTx <- NA #First value stays NA by definition.
dat_smo$LTx[2:3] <- 0 #by definition

assign_LTx_smo <- function(i){
 if (i>80){p_4pt <- p_pspt <- 0} #see comment at assign_stage_4_smo.
 dat_smo$LTx[i+1] <<- dat_smo[i,"LTx"] * (1 - dat_smo$mort4[i]) + dat_smo[i,"stage_4"] * p_4pt * (1-p_death_atltx) + dat_smo[i,"LVR"] * p_pspt * (1 - p_death_atltx)

}
for(i in 3:240){assign_LTx_smo(i=i)}
dat_smo$LTx #OK.

## End of the smoker state-transition section. Next: deaths, LYG, QALYs, costs, and the analogous former-smoker calculations.

#Calculate number of individuals for dat_smo$death in each cycle (=those who died in stage 1, 2, 3, 4, LVR, LTx, or simpler (-> absorbing state): all (100000) minus those who are in one of the other stages).

#death, smokers:
dat_smo$death <- NA

# Deaths are calculated residually as initial cohort size minus all living states, rather than by the long original spreadsheet death formula.
dat_smo$death <- 100000 - rowSums(dat_smo[,c("stage_1","stage_2","stage_3","stage_4","LVR","LTx")]) #OK

#LYG, smokers, not discounted:
dat_smo$LYG <- rowSums(dat_smo[,c("stage_1","stage_2","stage_3","stage_4","LVR","LTx")]) * cycle_length / 12
dat_smo$LYG[c(1,241)] <- 0

# Calculate QALYs from utilities per stage (u_1, u_2 etc.), probabilities per stage (p_mild1 etc.), and "d_mild1", "d_mild2", "d_mod2" etc. (=?) #

# In contrast to the simplified model, QALYs include mild AEs and the additional LVR/LTx states.
#QALY, smokers, not discounted:
#(I6*(u_1+p_mild1*d_mild1)+J6*(u_2+p_mild2*d_mild2+p_mod2*d_mod2+p_sev2*d_sev2)+K6*(u_3+p_mild3*d_mild3+p_mod3*d_mod3+p_sev3*d_sev3)+L6*(u_4+p_mild4*d_mild4+p_mod4*d_mod4+p_sev4*d_sev4)+M6*(u_lvr+p_mildL*d_mildL+p_modL*d_modL+p_sevL*d_sevL)+N6*u_ltx)*cycle_length/12
# All constants have already been defined; values from Excel columns I, J, K, L, M, and N are entered manually here for checking.
(dat_smo[2,"stage_1"]*(u_1+p_mild1*d_mild1) + dat_smo[2,"stage_2"] *(u_2 + p_mild2*d_mild2 + p_mod2*d_mod2 + p_sev2*d_sev2) + dat_smo[2,"stage_3"]*(u_3 + p_mild3*d_mild3 + p_mod3*d_mod3 + p_sev3*d_sev3) + dat_smo[2,"stage_4"]*(u_4 + p_mild4*d_mild4 + p_mod4*d_mod4 + p_sev4*d_sev4) + dat_smo[2,"LVR"]*(u_lvr + p_mildL*d_mildL + p_modL*d_modL + p_sevL*d_sevL) + dat_smo[2,"LTx"]*u_ltx)*cycle_length/12

qaly <- function(data, i){
 (data[i,"stage_1"]*(u_1 + p_mild1*d_mild1) + data[i,"stage_2"] *(u_2 + p_mild2*d_mild2 + p_mod2*d_mod2 + p_sev2*d_sev2) + data[i,"stage_3"]*(u_3 + p_mild3*d_mild3 + p_mod3*d_mod3 + p_sev3*d_sev3) + data[i,"stage_4"]*(u_4 + p_mild4*d_mild4 + p_mod4*d_mod4 + p_sev4*d_sev4) + data[i,"LVR"]*(u_lvr+p_mildL*d_mildL+p_modL*d_modL + p_sevL*d_sevL) + data[i,"LTx"]*u_ltx)*cycle_length/12
 }
qaly(data=dat_smo, i=2:241)

dat_smo$QALY <- NA
dat_smo$QALY[2:241] <- qaly(data=dat_smo, i=2:241)
dat_smo$QALY #Almost OK. The initial value is still missing and is calculated differently from the remaining values:
dat_smo$QALY[1] <- 100000*p_mild1*d_mild1*cycle_length/12

# Direct costs, smokers, not discounted:

# In contrast to the simplified model, direct costs include mild AEs and procedure-related LVR/LTx costs.
#costs, smokers, not discounted:
dat_smo$Kosten <- NA
(dat_smo$Kosten[1] <- 100000 * p_mild1 * c_mild1) #50221, OK. Excel: =I5*p_mild1*c_mild1

costs_smo <- Vectorize(function(data, i){
 if (i>80){p_4pt <- 0}  
 if (i<=80){p_pspt <- 0.00294}
 else {p_pspt <- 0}
 #dat_smo$costs[i] <- 
 data$stage_1[i]*(c_1+p_mild1*c_mild1) + data$stage_2[i]*(c_2 + p_mild2*c_mild2 + p_mod2*c_mod2 + p_sev2*c_sev2) + data$stage_3[i]*(c_3 + p_mild3*c_mild3 + p_mod3*c_mod3 + p_sev3*c_sev3) + data$stage_4[i]*(c_4 + p_mild4*c_mild4 + p_mod4*c_mod4 + p_sev4*c_sev4 + p_4ps*c_atlvr) + data$LVR[i]*(c_4+p_mildL*c_mild4+p_modL*c_mod4+p_sevL*c_sev4) + (data$stage_4[i]*p_4pt + data$LVR[i]*p_pspt)*c_atltx + data$LTx[i]*c_ltx
}, vectorize.args="i")
(x <- costs_smo(data=dat_smo, i=1:241)) #OK! Only difference: the final value (cycle 240) was manually set to 0 in the Excel sheet.
dat_smo$Kosten <- x
dat_smo$Kosten[1] <- 100000 * p_mild1 * c_mild1 #Reset the initial value.
dat_smo$Kosten[241] <- 0 #Set the final value to 0 as in the Excel sheet.

## Next: discounting, the analogous former-smoker calculations, probabilistic sensitivity analysis, CEAC, and ICER over time.

# Discounted LYG:
Cycle <- 0:240
discount_factor <- 1/(((1+disc_e)^(cycle_length/12))^Cycle)
(dat_smo$LYG_discounted <- dat_smo$LYG * discount_factor) #OK!

# Discounted QALYs:
(dat_smo$QALY_discounted <- dat_smo$QALY * discount_factor) #OK!

# Discounted costs:
(dat_smo$Kosten_discounted <- dat_smo$Kosten * discount_factor) #OK!

#doi.org/10.1186/s12962-019-0196-1
# Note from the cited paper: effects may be discounted slightly less than costs. The present model keeps equal 3% discounting, as in the original model.

# Discounted LYG under the alternative discounting note above:
(dat_smo$LYG_discounted2 <- dat_smo$LYG * discount_factor)
# Discounted QALYs under the alternative discounting note above:
(dat_smo$QALY_discounted2 <- dat_smo$QALY * discount_factor)
# Discounted costs under the alternative discounting note above:
(dat_smo$Kosten_discounted2 <- dat_smo$Kosten * discount_factor) 



###############
# Former smokers

N <- 241 #number of cycles
dat_fosmo <- data.frame(matrix(nrow=N, ncol=7))
colnames(dat_fosmo) <- c("stage_1","stage_2","stage_3","stage_4","LVR","LTx","death")
dat_fosmo[] <- 0
dat_fosmo[1,"stage_1"] <- 100000 #start with 100000 people in stage 1 (and 0 in the other stages)

#Define mortality in the different cycles for the different stages for former smokers: (see Table 2.7 in Menn's dissertation)
dat_fosmo[,c("mort1","mort2","mort3","mort4")] <- NA
dat_fosmo[,c("mort1","mort2","mort3","mort4")][1:20,] <- rep(c(0.0008,0.00136,0.00323,0.04022), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][21:40,] <- rep(c(0.00137,0.00234,0.00555,0.04022), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][41:60,] <- rep(c(0.00251,0.0043,0.01018,0.04022), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][61:80,] <- rep(c(0.0039,0.00668,0.01578,0.04022), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][81:100,] <- rep(c(0.00504,0.00863,0.02037,0.04022), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][101:120,] <- rep(c(0.00851,0.01455,0.03421,0.03421), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][121:140,] <- rep(c(0.01318,0.02249,0.05258,0.05258), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][141:160,] <- rep(c(0.0228,0.03877,0.08963,0.08963), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][161:180,] <- rep(c(0.03957,0.06686,0.15156,0.15156), each=20)
dat_fosmo[,c("mort1","mort2","mort3","mort4")][181:241,] <- rep(c(0.07423,0.12386,0.2695,0.2695), each=61)

# Use the intervention progression rates for former smokers (prog1_int etc.).

dat_fosmo$stage_1[2:21] <- dat_fosmo$stage_1[1] * (1 - prog1_int - dat_fosmo$mort1[1])^(1:20)

dat_fosmo$stage_1[22:41] <- dat_fosmo$stage_1[21] * (1 - prog1_int - dat_fosmo$mort1[21])^(1:20)
dat_fosmo$stage_1[42:61] <- dat_fosmo$stage_1[41] * (1 - prog1_int - dat_fosmo$mort1[41])^(1:20)
dat_fosmo$stage_1[62:81] <- dat_fosmo$stage_1[61] * (1 - prog1_int - dat_fosmo$mort1[61])^(1:20)
dat_fosmo$stage_1[82:101] <- dat_fosmo$stage_1[81] * (1 - prog1_int - dat_fosmo$mort1[81])^(1:20)
dat_fosmo$stage_1[102:121] <- dat_fosmo$stage_1[101] * (1 - prog1_int - dat_fosmo$mort1[101])^(1:20)
dat_fosmo$stage_1[122:141] <- dat_fosmo$stage_1[121] * (1 - prog1_int - dat_fosmo$mort1[121])^(1:20)
dat_fosmo$stage_1[142:161] <- dat_fosmo$stage_1[141] * (1 - prog1_int - dat_fosmo$mort1[141])^(1:20)
dat_fosmo$stage_1[162:181] <- dat_fosmo$stage_1[161] * (1 - prog1_int - dat_fosmo$mort1[161])^(1:20)
dat_fosmo$stage_1[182:241] <- dat_fosmo$stage_1[181] * (1 - prog1_int - dat_fosmo$mort1[181])^(1:60)

dat_fosmo$stage_1 #OK.

#Stage 2, former smokers:
dat_fosmo$stage_2 <- NA #1st value is NA by definition
dat_fosmo$stage_2[2] <- dat_fosmo[1,"stage_1"]*prog1_int #first cycle with N(stage 2) > 0 as starting point.
dat_fosmo$stage_2[3] <- dat_fosmo$stage_2[2]*(1-prog2_int-dat_fosmo$mort2[1]-p_sev2*0.13883)+dat_fosmo$stage_1[2]*prog1_int

assign_stage_2_fosmo <- function(i){
 dat_fosmo$stage_2[i+1] <<- dat_fosmo$stage_2[i]*(1-prog2_int-dat_fosmo$mort2[i]-p_sev2*ae_mort2_fosmo[i])+dat_fosmo$stage_1[i]*prog1_int
}
for(i in 3:240){assign_stage_2_fosmo(i=i)}
dat_fosmo$stage_2 #OK.

#Stage 3, former smokers:

dat_fosmo$stage_3 <- NA #1st value is NA by definition
#=AG5*(1-prog3_int-Z5-p_sev3*AD5)+AF5*prog2_int
dat_fosmo$stage_3[2] <- 0 #dat_fosmo$stage_3[1] * prog2_int #OK.
dat_fosmo$stage_3[3] <- dat_fosmo$stage_2[2] * prog2_int #OK.

assign_stage_3_fosmo <- function(i){ #assign different mortality rates according to the cycle
 dat_fosmo$stage_3[i+1] <<- dat_fosmo$stage_3[i]*(1-prog3_int-dat_fosmo$mort3[i]-p_sev3*ae_mort3_fosmo[i])+dat_fosmo$stage_2[i]*prog2_int
}
for(i in 3:240){assign_stage_3_fosmo(i=i)}
dat_fosmo$stage_3 #OK. 

#Stage 4, former smokers:
#=AH5*(1-p_4ps-p_4pt-Z5-p_sev4*AD5)+AG5*prog3_int

dat_fosmo$stage_4 <- NA

#=AH5*(1-p_4ps-p_4pt-Z5-p_sev4*AD5)+AG5*prog3_int
dat_fosmo$stage_4[2] <- 0
dat_fosmo$stage_4[3] <- dat_fosmo$stage_3[2] * prog3_int

assign_stage_4_fosmo <- function(i){ #assign different mortality rates according to the cycle
 if (i >80){p_4pt <- 0} #Note the arrow/comment in the Excel sheet.
 dat_fosmo$stage_4[i+1] <<- dat_fosmo$stage_4[i]*(1-p_4ps-p_4pt-dat_fosmo$mort3[i]-p_sev4*ae_mort3_fosmo[i])+dat_fosmo$stage_3[i]*prog3_int
}
for(i in 3:240){assign_stage_4_fosmo(i=i)}
dat_fosmo$stage_4 #OK!


# Additional transition only relevant for the original model: LVR in former smokers.
#Stage LVR, former smokers:
#=AI5*(1-p_pspt-Z5-p_sevL*AD5)+AH5*p_4ps*(1-p_lvr_death)

dat_fosmo$LVR <- NA #LVR[1] stays NA by definition
dat_fosmo$LVR[2] <- 0 #by definition
dat_fosmo$LVR[3] <- dat_fosmo$stage_4[3] * p_4ps*(1-p_lvr_death) #still 0.

assign_LVR_fosmo <- function(i){ #assign different mortality rates according to the cycle
 if (i >80){p_pspt <- 0} #Note the arrow/comment in the Excel sheet.
 dat_fosmo$LVR[i+1] <<- dat_fosmo$LVR[i]*(1-p_pspt-dat_fosmo$mort3[i]-p_sevL*ae_mort3_fosmo[i])+dat_fosmo$stage_4[i]*p_4ps*(1-p_lvr_death)
}
for(i in 3:240){assign_LVR_fosmo(i=i)}
dat_fosmo$LVR #OK!


# Additional transition only relevant for the original model: LTx in former smokers.
dat_fosmo$LTx <- NA #1st value is NA by definition
dat_fosmo$LTx[2:3] <- 0 #by definition

assign_LTx_fosmo <- function(i){
 if (i>80){p_4pt <- p_pspt <- 0}  
 dat_fosmo$LTx[i+1] <<- dat_fosmo[i,"LTx"] * (1 - dat_fosmo$mort4[i]) + dat_fosmo[i,"stage_4"] * p_4pt * (1-p_death_atltx) + dat_fosmo[i,"LVR"] * p_pspt * (1 - p_death_atltx)

}
for(i in 3:240){assign_LTx_fosmo(i=i)}
dat_fosmo$LTx #OK!

#Calculate number of individuals for dat_fosmo$death in each cycle (=those who died in stage 1, 2, 3, 4, LVR, LTx, or simpler (-> absorbing state): all (100000) minus those who are in one of the other stages).

dat_fosmo$death <- 100000 - rowSums(dat_fosmo[,c("stage_1","stage_2","stage_3","stage_4","LVR","LTx")]) #OK

#LYG, former smokers, not discounted:
(dat_fosmo$LYG <- rowSums(dat_fosmo[,c("stage_1","stage_2","stage_3","stage_4","LVR","LTx")]) * cycle_length / 12) #OK!
dat_fosmo$LYG[c(1,241)] <- 0

#QALYs:
dat_fosmo$QALY <- NA
dat_fosmo$QALY[2:241] <- qaly(data=dat_fosmo, i=2:241)
dat_fosmo$QALY[241] <- 0
dat_fosmo$QALY #OK! The initial value is still missing and is calculated differently from the remaining values:
dat_fosmo$QALY[1] <- 100000*p_mild1*d_mild1*cycle_length/12

#dir. costs:
dat_fosmo$Kosten <- costs_smo(data=dat_fosmo, i=1:241)
dat_fosmo$Kosten[1] <- 100000 * p_mild1 * c_mild1 #Reset the initial value.
dat_fosmo$Kosten[241] <- 0 #Set the final value to 0 as in the Excel sheet.


# Discounted LYG:
Cycle <- 0:240
(dat_fosmo$LYG_discounted <- dat_fosmo$LYG * discount_factor) #OK!

# Discounted QALYs:
(dat_fosmo$QALY_discounted <- dat_fosmo$QALY * discount_factor) #OK!

# Discounted costs:
(dat_fosmo$Kosten_discounted <- dat_fosmo$Kosten * discount_factor) #OK!


# Next: probabilistic sensitivity analysis, CEAC, and ICER over time.
# The probabilistic model uses formulas from the Excel worksheet "distributions" rather than manually copied random values from "prob Modell".

# Example from Excel: pb_prog1 was generated with BETAINV(RAND(); alpha; beta), yielding values around 0.014, close to the mean.

## Create parameter distributions using beta distributions and corresponding parameters.

# Note: alpha and beta can be calculated from the parameter mean and standard error.
# alpha = mean^2 * (1 - mean) / SE^2 - mean
# beta = mean * (1 - mean)^2 / SE^2 - (1 - mean)


# Compared with the simplified model, indirect-cost vectors also include mild AEs and LTx.
# Create vectors with indirect costs by stage and exacerbation type. They are used for both dat_smo and dat_fosmo and are therefore not attached to either data frame.
costs_stage2_stabil <- c(0, rep(112.50, 19), rep(376.37, 20), rep(667.24, 20), rep(237.22, 20), rep(0, 161))
costs_stage3_stabil <- c(0, rep(268.06, 19), rep(896.81, 20), rep(1589.9, 20), rep(565.24, 20), rep(0, 161))
costs_stage4_stabil <- c(0, rep(281.24, 19), rep(940.92, 20), rep(1668.1, 20), rep(593.04, 20), rep(0, 161))
costs_LTx_stabil <- c(0, rep(5059.01, 19), rep(4588.15, 20), rep(3746.88, 20), rep(1185.72, 20), rep(0, 161))
costs_stage1_leichte_AE <- c(rep(122.99, 20), rep(111.54, 20), rep(91.09, 20), rep(28.83, 20), rep(0, 161))
costs_stage2_leichte_AE <- costs_stage3_leichte_AE <- c(rep(245.97, 20), rep(223.08, 20), rep(182.18, 20), rep(57.65, 20), rep(0, 161)) #Interestingly, the Excel sheet uses the same costs for stage 2 and stage 3.
costs_stage4_leichte_AE <- c(rep(199.12, 20), rep(180.59, 20), rep(147.48, 20), rep(46.67, 20), rep(0, 161))
costs_stage2_mod_AE <- c(rep(491.94, 20), rep(446.16, 20), rep(364.35, 20), rep(115.30, 20), rep(0, 161))
costs_stage3_mod_AE <- c(rep(550.51, 20), rep(499.27, 20), rep(407.73, 20), rep(129.03, 20), rep(0, 161))
costs_stage4_mod_AE <- c(rep(491.94, 20), rep(446.16, 20), rep(364.35, 20), rep(115.30, 20), rep(0, 161))
costs_stage2_schwere_AE <- c(rep(367, 20), rep(332.85, 20), rep(271.82, 20), rep(86.02, 20), rep(0, 161))
costs_stage3_schwere_AE <- c(rep(663.73, 20), rep(601.96, 20), rep(491.58, 20), rep(155.56, 20), rep(0,161))
costs_stage4_schwere_AE <- c(rep(648.12, 20), rep(587.79, 20), rep(480.02, 20), rep(151.9, 20), rep(0,161))

##########################################################
# Next: indirect costs for smokers (undiscounted and discounted).
##########################################################

# Smokers, indirect costs, stage 1 (not discounted):
#=      I5             * p_mild1 * sum(AX5:AX84)
dat_smo$stage_1_ind_costs <- NA
dat_smo$stage_1_ind_costs[1:20] <- dat_smo$stage_1[1:20]   * p_mild1 * 122.99 #Indirect costs decrease in later cycles because a larger share of the cohort has retired.
dat_smo$stage_1_ind_costs[21:40] <- dat_smo$stage_1[21:40] * p_mild1 * 111.54
dat_smo$stage_1_ind_costs[41:60] <- dat_smo$stage_1[41:60] * p_mild1 * 91.09
dat_smo$stage_1_ind_costs[61:80] <- dat_smo$stage_1[61:80] * p_mild1 * 28.83
dat_smo$stage_1_ind_costs[81:241] <- 0
dat_smo$stage_1_ind_costs #OK!

# Continue with stages 2-4, LVR, and LTx.


# Smokers, indirect costs, stage 2:
#=J6*(AT6+p_mild2*AY6+p_mod2*BB6+p_sev2*BE6)
(dat_smo$stage_2_ind_costs <- dat_smo$stage_2 * (costs_stage2_stabil + p_mild2*costs_stage2_leichte_AE + p_mod2*costs_stage2_mod_AE + p_sev2*costs_stage2_schwere_AE)) #OK.

# Smokers, indirect costs, stage 3:
#=K6*(AU6+p_mild3*AZ6+p_mod3*BC6+p_sev3*BF6)
(dat_smo$stage_3_ind_costs <- dat_smo$stage_3 * (costs_stage3_stabil + p_mild3*costs_stage3_leichte_AE + p_mod3*costs_stage3_mod_AE + p_sev3*costs_stage3_schwere_AE)) #OK.

# Smokers, indirect costs, stage 4:
#=L6*(AV6+p_mild4*BA6+p_mod4*BD6+p_sev4*BG6)
(dat_smo$stage_4_ind_costs <- dat_smo$stage_4 * (costs_stage4_stabil + p_mild4*costs_stage4_leichte_AE + p_mod4*costs_stage4_mod_AE + p_sev4*costs_stage4_schwere_AE)) #OK.

# Smokers, indirect costs, LVR:
#=M6*(AV6+p_mildL*BA6+p_modL*BD6+p_sevL*BG6)
(dat_smo$LVR_ind_costs <- dat_smo$LVR * (costs_stage4_stabil + p_mildL*costs_stage4_leichte_AE + p_modL*costs_stage4_mod_AE + p_sevL*costs_stage4_schwere_AE)) #OK.

# Smokers, indirect costs, LTx:
#=N6*AW6
(dat_smo$LTx_ind_costs <- dat_smo$LTx * costs_LTx_stabil) #OK.

#indirect costs smokers, not discounted:
(dat_smo$ind_costs <- rowSums(dat_smo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs","LVR_ind_costs","LTx_ind_costs")], na.rm=T)) #OK.

#indirect costs smokers, discounted:
(dat_smo$ind_costs_discounted <- dat_smo$ind_costs * discount_factor) #OK.

#indirect costs, former smokers, stage 1:
(dat_fosmo$stage_1_ind_costs <- dat_fosmo$stage_1 * p_mild1 * costs_stage1_leichte_AE) #OK.

#indirect costs, former smokers, stage 2:
#=AF6*(AT6+p_mild2*AY6+p_mod2*BB6+p_sev2*BE6)
(dat_fosmo$stage_2_ind_costs <- dat_fosmo$stage_2 * (costs_stage2_stabil + p_mild2*costs_stage2_leichte_AE + p_mod2*costs_stage2_mod_AE+p_sev2*costs_stage2_schwere_AE)) #OK.

#indirect costs, former smokers, stage 3:
#=AG6*(AU6+p_mild3*AZ6+p_mod3*BC6+p_sev3*BF6)
(dat_fosmo$stage_3_ind_costs <- dat_fosmo$stage_3 * (costs_stage3_stabil + p_mild3*costs_stage3_leichte_AE + p_mod3*costs_stage3_mod_AE+p_sev3*costs_stage3_schwere_AE)) #OK.

#indirect costs, former smokers, stage 4:
#=AH6*(AV6+p_mild4*BA6+p_mod4*BD6+p_sev4*BG6)
(dat_fosmo$stage_4_ind_costs <- dat_fosmo$stage_4 * (costs_stage4_stabil + p_mild4*costs_stage4_leichte_AE + p_mod4*costs_stage4_mod_AE+p_sev4*costs_stage4_schwere_AE)) #OK.

#indirect costs, former smokers, LVR:
#=AI6*(AV6+p_mildL*BA6+p_modL*BD6+p_sevL*BG6)
(dat_fosmo$LVR_ind_costs <- dat_fosmo$LVR*(costs_stage4_stabil + p_mildL*costs_stage4_leichte_AE + p_modL*costs_stage4_mod_AE + p_sevL*costs_stage4_schwere_AE)) #OK.

#indirect costs, former smokers, LTx:
#=AJ6*AW6
(dat_fosmo$LTx_ind_costs <- dat_fosmo$LTx * costs_LTx_stabil) #OK.

#indirect costs, former smokers, total:
(dat_fosmo$ind_costs <- rowSums(dat_fosmo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs","LVR_ind_costs","LTx_ind_costs")], na.rm=TRUE)) #OK.

#indirect costs, former smokers, total, discounted:
(dat_fosmo$ind_costs_discounted <- dat_fosmo$ind_costs * discount_factor) #OK.

### ICER over time
## Direct + indirect costs

# Smokers:

#=(V5+BO5)/$J$1
(dat_smo$IKER <- (dat_smo$Kosten_discounted + dat_smo$ind_costs_discounted)/100000) #OK.

# Former smokers:
(dat_fosmo$IKER <- (dat_fosmo$Kosten_discounted + dat_fosmo$ind_costs_discounted)/100000) #OK.

# Costs for the intervention arm ["Int"]:
c_int + dat_smo$IKER * (1-eff_int) + dat_fosmo$IKER*eff_int
dat_fosmo$Costs_Int[1] <- (c_int + dat_smo$IKER * (1-eff_int) + dat_fosmo$IKER*eff_int)[1] #OK.

dat_fosmo$Costs_Int[2] <- dat_fosmo$Costs_Int[1] + dat_smo$IKER[2]*(1-eff_int) + dat_fosmo$IKER[2]*eff_int #OK.
#Remaining cycles via helper function:
assign_ICER_int <- function(i){dat_fosmo$Costs_Int[i] <<- dat_fosmo$Costs_Int[i-1] + dat_smo$IKER[i]*(1-eff_int) + dat_fosmo$IKER[i]*eff_int}
for(i in 3:241){assign_ICER_int(i)}
dat_fosmo$Costs_Int #OK.

# Costs for the standard-care/usual-care arm ["ST"]:
dat_fosmo$ST[1] <- dat_smo$IKER[1]*(1-eff_uc) + dat_fosmo$IKER[1]*eff_uc #CG5*(1-eff_uc)+CH5*eff_uc
dat_fosmo$ST[2] <- dat_fosmo$ST[1] + dat_smo$IKER[2]*(1-eff_uc) + dat_fosmo$IKER[2]*eff_uc #OK.

assign_costs_ST <- function(i){
 dat_fosmo$ST[i+1] <<- dat_fosmo$ST[i] + dat_smo$IKER[i+1]*(1-eff_uc)+dat_fosmo$IKER[i+1]*eff_uc 
}
for(i in 2:240){assign_costs_ST(i)}
dat_fosmo$ST #OK.

# Incremental costs: intervention costs minus standard-care costs.
dat_fosmo$inkr_Kosten <- dat_fosmo$Costs_Int - dat_fosmo$ST #OK

## Effects:

#smokers:

#=U5/$J$1, =U6/$J$1 u.s.w. = QALY_discounted smokers / 100000
(dat_smo$effects <- dat_smo$QALY_discounted / 100000) #OK

#former smokers:
#=AQ5/$J$1, =AQ6/$J$1 u.s.a. = QALY_discounted former smokers / 100000
(dat_fosmo$effects <- dat_fosmo$QALY_discounted / 100000) #OK

# Intervention arm:

#=CL5*(1-eff_int)+CM5*eff_int #Startwert
#=CN5+CL6*(1-eff_int)+CM6*eff_int #first regular value
#=CN6+CL7*(1-eff_int)+CM7*eff_int #same pattern thereafter

cycle_eff_int <- NA
(cycle_eff_int[1] <- dat_smo$effects[1] * (1-eff_int) + dat_fosmo$effects[1] * eff_int)
(cycle_eff_int[2] <- cycle_eff_int[1] + dat_smo$effects[2]*(1-eff_int) + dat_fosmo$effects[2]*eff_int) #Second value also OK. Same pattern from here onward.

assign_cycle_eff_int <- function(i){
 cycle_eff_int[i+1] <<- cycle_eff_int[i] + dat_smo$effects[i+1]*(1-eff_int) + dat_fosmo$effects[i+1]*eff_int
}
for(i in 2:240){assign_cycle_eff_int(i)}
cycle_eff_int #OK

cycle_eff_st <- NA
cycle_eff_st[1] <- dat_smo$effects[1]*(1-eff_uc) + dat_fosmo$effects[1]*eff_uc #=CL5*(1-eff_uc)+CM5*eff_uc #OK
cycle_eff_st[2] <- cycle_eff_st[1] + dat_smo$effects[2]*(1-eff_uc) + dat_fosmo$effects[2]*eff_uc  #=CN5+CL6*(1-eff_int)+CM6*eff_int #OK. Same pattern from here onward.

assign_cycle_eff_st <- function(i){
 cycle_eff_st[i+1] <<- cycle_eff_st[i] + dat_smo$effects[i+1]*(1-eff_uc) + dat_fosmo$effects[i+1]*eff_uc
}
for(i in 2:240){assign_cycle_eff_st(i)}
cycle_eff_st #OK.

(incr_eff <- cycle_eff_int - cycle_eff_st) #OK

# ICER: incremental costs divided by incremental effects.
# Excel-specific ICER rounding: up to cycle 28, effects are rounded to four decimals; from cycle 29 onward, effects are rounded to two decimals.
# (ICER_exact <- dat_fosmo$inkr_costs[-1] / incr_eff)
# (ICER_round_effects_3_decimals <- dat_fosmo$inkr_costs[-1] / round(incr_eff,3))
# (ICER_round_effects_4_decimals <- dat_fosmo$inkr_costs[-1] / round(incr_eff,4))
# (ICER_round_effects_5_decimals <- dat_fosmo$inkr_costs[-1] / round(incr_eff,5))
# (ICER_round_effects_6_decimals <- dat_fosmo$inkr_costs[-1] / round(incr_eff,6))
(ICER_round_Menn_1 <- dat_fosmo$inkr_Kosten[1:29] / round(incr_eff[1:29],4)) #OK!
(ICER_round_Menn_2 <- dat_fosmo$inkr_Kosten[30:89] / round(incr_eff[30:89],2)) #OK!
(ICER_round_Menn_3 <- dat_fosmo$inkr_Kosten[90:240] / round(incr_eff[90:240],2)) #OK!
(ICER <- ICER_round_Menn <- c(ICER_round_Menn_1, ICER_round_Menn_2, ICER_round_Menn_3)) #OK!

# ICER over time:
plot_x <- ICER
plot(plot_x, type="l", ylim=c(-20000,50000), xlab="Zyklus", ylab="IKER (€/QALY)", main="IKER im Zeitverlauf", lwd=2, xaxt="n", xaxs="i", yaxs="i")
axis(1, seq(0,220,20))
abline(h=seq(-10000, 50000, 10000)) #OK

# Reproduction of Excel worksheet "Ergebnis_xls" with LYG, direct costs, indirect costs, and summary results:

### Smokers ###

## Total LYG, "Standard":
sum(dat_smo$LYG) #2186478, matches exactly.
## LYG per patient, "Standard":
sum(dat_smo$LYG) / 100000 #21.86478, OK (trivial)
## Total QALYs, "Standard":
sum(dat_smo$QALY) #1777840, matches exactly.
# Per patient:
sum(dat_smo$QALY) / 100000 #17.7784, OK (trivial)
## Total costs, "Standard":
sum(dat_smo$Kosten) #1874785113, matches exactly.
# Per patient:
sum(dat_smo$Kosten) / 100000 #18748, OK (trivial)
## Total LYG, discounted:
sum(dat_smo$LYG_discounted) #1531246, matches exactly.
# Per patient:
sum(dat_smo$LYG_discounted) / 100000 #15.31246, OK (trivial)
## Total QALYs, discounted:
sum(dat_smo$QALY_discounted) #1249561, matches exactly.
# Per patient:
sum(dat_smo$QALY_discounted) / 100000 #12.49561, OK (trivial)
## Total costs, discounted:
sum(dat_smo$Kosten_discounted) #1234340696, matches exactly.
# Per patient:
sum(dat_smo$Kosten_discounted) / 100000 #12343.41, OK.
## Indirect costs, "Standard":
sum(dat_smo$ind_costs) #1406894457, matches exactly.
# Per patient:
sum(dat_smo$ind_costs) / 100000 #14068.94, OK (trivial)
## Indirect costs, discounted:
sum(dat_smo$ind_costs_discounted) #1003696787, OK.
# Per patient:
sum(dat_smo$ind_costs_discounted) / 100000 #10036.97, OK (trivial)


### Former smokers ###

## Total LYG, "Standard":
sum(dat_fosmo$LYG, na.rm=T) #3020713, matches exactly.
## LYG per patient, "Standard":
sum(dat_fosmo$LYG, na.rm=T) / 100000 #30.20713, OK (trivial)
## Total QALYs, "Standard":
sum(dat_fosmo$QALY) #2519515, matches exactly.
# Per patient:
sum(dat_fosmo$QALY) / 100000 #25.19515, OK (trivial)
## Total costs, "Standard":
sum(dat_fosmo$Kosten) #1504653567, OK. Corrected on 2022-11-17; the Excel sheet fixes the final value at 0, which had not been done in R before.
# Per patient:
sum(dat_fosmo$Kosten) / 100000 #15046.62, OK (trivial)
## Total LYG, discounted:
sum(dat_fosmo$LYG_discounted, na.rm=T) #1903967, matches exactly.
# Per patient:
sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 #19.03967, OK (trivial)
## Total QALYs, discounted:
sum(dat_fosmo$QALY_discounted) #1589293, matches exactly.
# Per patient:
sum(dat_fosmo$QALY_discounted) / 100000 #15.89293, OK (trivial)
## Total costs, discounted:
sum(dat_fosmo$Kosten_discounted) #929724839, OK. Corrected on 2022-11-17; the Excel sheet fixes the final value at 0, which had not been done in R before.
# Per patient:
sum(dat_fosmo$Kosten_discounted) / 100000 #9297, OK.
## Indirect costs, "Standard":
sum(dat_fosmo$ind_costs) #331883115, matches exactly.
# Per patient:
sum(dat_fosmo$ind_costs) / 100000 #3318.831, OK (trivial)
## Indirect costs, discounted:
sum(dat_fosmo$ind_costs_discounted) #238700095, OK
# Per patient:
sum(dat_fosmo$ind_costs_discounted) / 100000 #2387.001, OK (trivial)

## "Light-blue boxes" with summary results:

# Direct costs, standard care (ST):
#=V$248*(1-eff_uc)+AR$248*eff_uc+c_uc
(total_dir_costs_standard <- sum(dat_smo$Kosten_discounted)/100000 * (1-eff_uc) + sum(dat_fosmo$Kosten_discounted)/100000 * eff_uc + c_uc) #12161, OK.

# Direct costs, intervention (Int):
#=V$248*(1-eff_int)+AR$248*eff_int+c_int
(total_dir_costs_intervention <- sum(dat_smo$Kosten_discounted)/100000 * (1-eff_int) + sum(dat_fosmo$Kosten_discounted)/100000 * eff_int + c_int) #12269, OK.

# LYG, standard care (ST):
#=T$248*(1-Werte!$C$3)+AP$248*Werte!$C$3
sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.06 #15.54, OK.

# LYG, intervention:
#=T$248*(1-Werte!$C$2)+AP$248*Werte!$C$2
sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.22 #16.13, OK.

(total_dir_costs_standard - total_dir_costs_intervention) / ((sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.06)-(sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.22))
#182 EUR/LYG, almost identical to Excel (183 EUR).

# QALYs, standard care (ST):
#=U$248*(1-Werte!$C$3)+AQ$248*Werte!$C$3
(total_QALY_standard <- sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.06) #12.70, OK.

# QALYs, intervention:
#=U$248*(1-Werte!$C$2)+AQ$248*Werte!$C$2
(total_QALY_intervention <- sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.22) #13.24, OK.

# Indirect costs, standard care (ST):
#=(V$248+BO$248)*(1-eff_uc)+(AR$248+BX$248)*eff_uc+c_uc
(total_indir_costs_standard <- (sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc) #21739, OK.

# Indirect costs, intervention (Int):
(total_indir_costs_intervention <- (sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_int) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int)
#20623, OK.

#Overall ICER (delta costs / delta effects) for direct costs:
(total_dir_costs_intervention - total_dir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #199.6697, OK.

#Overall ICER (delta costs / delta effects) for indirect costs:
(total_indir_costs_intervention - total_indir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #-2052.097, OK.


#########################################
# Create probabilistic parameters:
# Needed later for probabilistic sensitivity analyses; see QALYs+costs_probab_full.r.

## Create parameter distributions using beta distributions and corresponding parameters.
# Then check if the 0.025 and 0.975 quantiles match the lower and upper 95% CI in the Excel sheet.

# Note: alpha and beta can be calculated from the parameter mean and standard error.
# alpha = mean^2 * (1 - mean) / SE^2 - mean
# beta = mean * (1 - mean)^2 / SE^2 - (1 - mean)

# For reproducibility and future updating, derive alpha and beta from the mean and SE rather than copying them from Excel.

#pb_prog1:
#Attempt to reproduce alpha (xlsx) for pb_prog1 exactly using the formula:
(alpha_pb_prog1 <- prog1^2*(1-prog1)/0.00053^2-prog1) #690.9052
calc_alpha <- function(mu, se){mu^2*(1-mu)/se^2-mu}
calc_alpha(mu=prog1, se=0.00053) #OK

# Calculate beta analogously with a helper function:
calc_beta <- function(x, se){x*(1-x)^2/se^2-(1-x)}
calc_beta(x=prog1, se=0.00053) #OK
(beta_pb_prog1 <- calc_beta(x=prog1, se=0.00053))

set.seed(1)
pb_prog1 <- rbeta(100000, calc_alpha(mu=prog1, se=0.00053), calc_beta(x=prog1, se=0.00053))
quantile(pb_prog1, probs=c(0.025,0.975)) #=the CIs from Excel (lower CI = 0.0130 = 1.30%, upper CI = 0.0151 = 1.51%)


## Calculate gamma-distribution parameters for costs in a helper function:
#c_1 = 102.64
#c_1_sd = 20.00

#mean = k * phi
#variance = k * phi^2

#102.64 = k * phi
#20.00 = k * phi^2

#(I) in (II): 102.64/phi * phi^2 = 20

#k = 102.64^2/400 = 26.33742 = alpha = first parameter of the gamma distribution for pb_c_1.
#phi = 400/102.64 = 3.897116 = beta(xls) = second parameter of the gamma distribution for pb_c_1.

gamma_params <- function(mean, sd){
 k <- mean^2/sd^2
 phi_xls <- mean/k
 phi_tre <- mean/sd^2
 res <- list(k, phi_xls, phi_tre)
 names(res) <- c("alpha","beta_xls", "beta_tre")
 return(res)
}
gamma_params(mean=102.64, sd=20.00) #alpha=26.33742, beta_xls=3.897116, beta_tre=0.2566
# This yields alpha and beta_xls.
# alpha is identical in the Excel and R parameterisations; beta_xls is the Excel scale parameter.
# beta_tre is the R rate parameter used for the actual distribution.

#pb_eff_int:
set.seed(1)
pb_eff_int <- rbeta(100000, calc_alpha(mu=eff_int, se=0.0700), calc_beta(x=eff_int, se=0.0700))
quantile(pb_eff_int, probs=c(0.025,0.975)) #OK

#pb_eff_uc:
set.seed(1)
pb_eff_uc <- rbeta(100000, calc_alpha(mu=eff_uc, se=0.0200), calc_beta(x=eff_uc, se=0.0200))
quantile(pb_eff_uc, probs=c(0.025,0.975)) #OK

#pb_c_int:
x <- gamma_params(mean=595.92, sd=50.00)
set.seed(1)
pb_c_int <- rgamma(100000, x$alpha, x$beta_tre)
quantile(pb_c_int, probs=c(0.025,0.975)) #OK


# How to create a CEAC in R: https://search.r-project.org/CRAN/refmans/BCEA/html/bcea.html
library(BCEA)
data(Vaccine)
head(eff, 20) #Effects (QALYs) for status quo (left) and intervention (right), in matrix form. Each row corresponds to one simulation.
head(cost, 20) #Costs for status quo (left) and intervention (right), in matrix form. Each row corresponds to the associated costs from one simulation.
m <- bcea(e=eff, c=cost, ref=2) #ref=2 if the first column is the status quo and the second column is the intervention (usually more costly and more effective than the status quo).
plot(m)
dev.off()
summary(m)
# This suggests a reduced probabilistic calculation with all steps up to cycle_eff_int, cycle_eff_st, Costs_Int, and ST inside one function.

# Compared with the simplified model, the PSA includes uncertainty for mild AEs, LVR, LTx,
# and procedure-related mortality and costs.
############ Probabilistic sensitivity analysis ############

Pb_prog1 <- function(){rbeta(1, calc_alpha(mu=prog1, se=0.00053), calc_beta(x=prog1, se=0.00053))} #replaces pb_prog1
Pb_prog2 <- function(){rbeta(1, calc_alpha(mu=prog2, se=0.00034), calc_beta(x=prog2, se=0.00034))} #replaces pb_prog2
Pb_prog3 <- function(){rbeta(1, calc_alpha(mu=prog3, se=0.00162), calc_beta(x=prog3, se=0.00162))} #replaces pb_prog3
Pb_p_LVR <- function(){rbeta(1, calc_alpha(mu=0.02729, se=0.01), calc_beta(x=0.02729, se=0.01))} #replaces pb_p_LVR
Pb_p_LTx <- function(){rbeta(1, calc_alpha(mu=0.00294, se=0.001), calc_beta(x=0.00294, se=0.001))} #replaces pb_p_LTx
Pb_p_death_atltx <- function(){rbeta(1, calc_alpha(mu=0.06, se=0.02), calc_beta(x=0.06, se=0.02))}
Pb_prog1_int <- function(){rbeta(1, calc_alpha(mu=prog1_int, se=0.00036), calc_beta(x=prog1_int, se=0.00036))} #replaces pb_prog1_int
Pb_prog2_int <- function(){rbeta(1, calc_alpha(mu=prog2_int, se=0.00024), calc_beta(x=prog2_int, se=0.00024))} #replaces pb_prog2_int
Pb_prog3_int <- function(){rbeta(1, calc_alpha(mu=prog3_int, se=0.00036), calc_beta(x=prog3_int, se=0.00036))} #replaces pb_prog3_int
Pb_p_mild1 <- function(){rbeta(1, calc_alpha(mu=p_mild1, se=0.0361), calc_beta(x=p_mild1, se=0.0361))} #replaces pb_p_mild1
Pb_p_mild2 <- function(){rbeta(1, calc_alpha(mu=p_mild2, se=0.00511), calc_beta(x=p_mild2, se=0.00511))} #replaces pb_p_mild2
Pb_p_mild3 <- function(){rbeta(1, calc_alpha(mu=p_mild3, se=0.00762), calc_beta(x=p_mild3, se=0.00762))} #replaces pb_p_mild3
Pb_p_mild4 <- function(){rbeta(1, calc_alpha(mu=p_mild4, se=0.00531), calc_beta(x=p_mild4, se=0.00531))} #replaces pb_p_mild4
Pb_p_mildL <- function(){rbeta(1, calc_alpha(mu=p_mildL, se=0.00557), calc_beta(x=p_mildL, se=0.00557))} #replaces pb_p_mildL
Pb_p_mod2 <- function(){rbeta(1, calc_alpha(mu=p_mod2, se=0.01090), calc_beta(x=p_mod2, se=0.01090))} #replaces pb_p_mod2
Pb_p_mod3 <- function(){rbeta(1, calc_alpha(mu=p_mod3, se=0.01412), calc_beta(x=p_mod3, se=0.01412))} #replaces pb_p_mod3
Pb_p_mod4 <- function(){rbeta(1, calc_alpha(mu=p_mod4, se=0.00956), calc_beta(x=p_mod4, se=0.00956))} #replaces pb_p_mod4
Pb_p_modL <- function(){rbeta(1, calc_alpha(mu=p_modL, se=0.01055), calc_beta(x=p_modL, se=0.01055))} #replaces pb_p_modL
Pb_p_sev2 <- function(){rbeta(1, calc_alpha(mu=p_sev2, se=0.00246), calc_beta(x=p_sev2, se=0.00246))} #replaces pb_p_sev2
Pb_p_sev3 <- function(){rbeta(1, calc_alpha(mu=p_sev3, se=0.00225), calc_beta(x=p_sev3, se=0.00225))} #replaces pb_p_sev3
Pb_p_sev4 <- function(){rbeta(1, calc_alpha(mu=p_sev4, se=0.00189), calc_beta(x=p_sev4, se=0.00189))} #replaces pb_p_sev4
Pb_p_sevL <- function(){rbeta(1, calc_alpha(mu=p_sevL, se=0.00192), calc_beta(x=p_sevL, se=0.00192))} #replaces pb_p_sevL
Pb_p_4ps <- function(){rbeta(1, calc_alpha(mu=p_4ps, se=0.0100), calc_beta(x=p_4ps, se=0.0100))} #replaces pb_p_4ps
Pb_p_4pt <- function(i){if(i<=80){rbeta(1, calc_alpha(mu=p_4pt, se=0.0010), calc_beta(x=p_4pt, se=0.0010))} else 0} #replaces pb_p_4pt
Pb_p_pspt <- function(i){if(i<=80){rbeta(1, calc_alpha(mu=p_4pt, se=0.0010), calc_beta(x=p_4pt, se=0.0010))} else 0} #replaces pb_p_pspt
Pb_u_1 <- function(){rbeta(1, calc_alpha(mu=u_1, se=0.0300), calc_beta(x=u_1, se=0.0300))} #replaces pb_u_1
Pb_u_2 <- function(){rbeta(1, calc_alpha(mu=u_2, se=0.0100), calc_beta(x=u_2, se=0.0100))} #replaces pb_u_2
Pb_u_3 <- function(){rbeta(1, calc_alpha(mu=u_3, se=0.0100), calc_beta(x=u_3, se=0.0100))} #replaces pb_u_3
Pb_u_4 <- function(){rbeta(1, calc_alpha(mu=u_4, se=0.0300), calc_beta(x=u_4, se=0.0300))} #replaces pb_u_4
Pb_u_lvr <- function(){rbeta(1, calc_alpha(mu=u_lvr, se=0.0500), calc_beta(x=u_lvr, se=0.0500))} #replaces pb_u_lvr
Pb_u_1tx <- function(){rbeta(1, calc_alpha(mu=u_ltx, se=0.0500), calc_beta(x=u_ltx, se=0.0500))} #replaces pb_u_ltx
Pb_p_lvr_death <- function(){rbeta(1, calc_alpha(mu=p_lvr_death, se=0.01), calc_beta(x=0.03, se=0.01))} #replaces pb_p_lvr_death
Pb_p_death_atltx <- function(){rbeta(1, calc_alpha(mu=p_death_atltx, se=0.02), calc_beta(x=0.06, se=0.02))} #replaces pb_p_death_atltx
Pb_eff_int <- function(){rbeta(1, calc_alpha(mu=eff_int, se=0.0700), calc_beta(x=eff_int, se=0.0700))} #replaces pb_eff_int
Pb_eff_uc <- function(){rbeta(1, calc_alpha(mu=eff_uc, se=0.0200), calc_beta(x=eff_uc, se=0.0200))} #replaces pb_eff_uc

# Cost parameters:
Pb_c_1 <- function(){
 x <- gamma_params(mean=102.64, sd=20.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_1

Pb_c_2 <- function(){
 x <- gamma_params(mean=184.69, sd=25.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_2

Pb_c_3 <- function(){
 x <- gamma_params(mean=366.57, sd=30.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_3

Pb_c_4 <- function(){
 x <- gamma_params(mean=431.35, sd=35.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_4

Pb_c_mild1 <- function(){
 x <- gamma_params(mean=6.950, sd=5.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mild1

Pb_c_mild2 <- function(){
 x <- gamma_params(mean=8.690, sd=5.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mild2

Pb_c_mild3 <- function(){
 x <- gamma_params(mean=21.73, sd=5.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mild3

Pb_c_mild4 <- function(){
 x <- gamma_params(mean=20.35, sd=5.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mild4

Pb_c_mod2 <- function(){
 x <- gamma_params(mean=39.86, sd=10.009)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mod2

Pb_c_mod3 <- function(){
 x <- gamma_params(mean=52.90, sd=10.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mod3

Pb_c_mod4 <- function(){
 x <- gamma_params(mean=51.52, sd=10.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_mod4

Pb_c_sev2 <- function(){
 x <- gamma_params(mean=2923.94, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_sev2

Pb_c_sev3 <- function(){
 x <- gamma_params(mean=3109.05, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_sev3

Pb_c_sev4 <- function(){
 x <- gamma_params(mean=3322.59, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_sev4

Pb_c_atlvr <- function(){
 x <- gamma_params(mean=6413.88, sd=1000.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_atlvr

Pb_c_atltx <- function(){
 x <- gamma_params(mean=65494.21, sd=5000.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_atltx

Pb_c_ltx <- function(){
 x <- gamma_params(mean=8269.02, sd=1000.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_ltx

Pb_c_int <- function(){
 x <- gamma_params(mean=595.92, sd=50.00)
 rgamma(1, x$alpha, x$beta_tre)} #replaces pb_c_int

pb_c_uc <- rep(0, 100000)

discount <- function(x, disc_e=0.03){x * (1/(((1+disc_e)^(cycle_length/12))^Cycle))}

# Recreate qaly_smo_probab_new etc. using the functions just defined:
library(doParallel)

n_sim <- 1000
n_cycles <- 240
cl <- makeCluster(round(0.4*parallel::detectCores())) # Uses 10 cores on the workstation and 3 on the work laptop. Could be optimized further: when many cores are available, proportionally fewer should be used because of overhead.
#rbind(1:24, round((1:24)^0.7))

registerDoParallel(cl)

# Parallel loop for simulations.
# Compared with the simplified model, each probabilistic simulation carries additional LVR and LTx state vectors.
qaly_smo_probab <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- LVR <- LTx <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_smo$stage_1[1]
  for(j in 1:n_cycles){
    stage_1[j+1] <- stage_1[j]*(1 - Pb_prog1() - dat_smo$mort1[j])
  }
  
  # Stage 2 initialization and loop:
  stage_2[2] <- stage_1[1]*Pb_prog1()
  stage_2[3] <- stage_2[2]*(1 - Pb_prog2() - dat_smo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2]*Pb_prog1()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j]*(1 - Pb_prog2() - dat_smo$mort2[j] - Pb_p_sev2()*ae_mort2[j]) + stage_1[j]*Pb_prog1()
  }
  
  # Stage 3 initialization and loop:
  stage_3[3] <- stage_2[2]*Pb_prog2()
  for(j in 3:n_cycles){
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3() - dat_smo$mort3[j] - Pb_p_sev3()*ae_mort3[j]) + stage_2[j]*Pb_prog2()
  }
  
  # Stage 4 initialization and loop:
  stage_4[4] <- stage_3[4]*Pb_prog3()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j]*(1 - dat_smo$mort3[j] - Pb_p_sev4()*ae_mort3[j]) + stage_3[j]*Pb_prog3()
  }
  
  # LVR and LTx calculations:
  LVR[4] <- stage_4[3] * Pb_p_4ps() * (1 - Pb_p_lvr_death())
  for(j in 4:n_cycles){
    LVR[j+1] <- LVR[j] * (1 - Pb_p_pspt(j) - dat_smo$mort3[j] - Pb_p_sevL()*ae_mort3[j]) + stage_4[j] * Pb_p_4ps() * (1 - Pb_p_lvr_death())
  }
  
  for(j in 3:n_cycles){
    LTx[j+1] <- LTx[j] * (1 - dat_smo$mort4[j]) + stage_4[j] * Pb_p_4pt(j) * (1 - Pb_p_death_atltx()) + LVR[j] * Pb_p_pspt(j) * (1 - Pb_p_death_atltx())
  }
  
# Final QALY calculation.
  QALY <- ((stage_1*(Pb_u_1() + Pb_p_mild1()*d_mild1) +
            stage_2*(Pb_u_2() + Pb_p_mild2()*d_mild2 + Pb_p_mod2()*d_mod2 + Pb_p_sev2()*d_sev2)) +
           stage_3*(Pb_u_3() + Pb_p_mild3()*d_mild3 + Pb_p_mod3()*d_mod3 + Pb_p_sev3()*d_sev3) +
           stage_4*(Pb_u_4() + Pb_p_mild4()*d_mild4 + Pb_p_mod4()*d_mod4 + Pb_p_sev4()*d_sev4) +
           LVR*(Pb_u_lvr() + Pb_p_mildL()*d_mildL + Pb_p_modL()*d_modL + Pb_p_sevL()*d_sevL) +
           LTx*Pb_u_1tx()) * cycle_length/12
  
  QALY
}
qaly_smo_probab

for(i in 1:1000){qaly_smo_probab[1,i] <- 100000 * Pb_p_mild1() * c_mild1}

quantile(colSums(qaly_smo_probab), probs=c(0.025,0.975))
hist(colSums(qaly_smo_probab)) #loOKs OK (largely normally distributed).
dev.off()

qaly_smo_probab_discounted <- discount(qaly_smo_probab)

# Next: do the analogous calculations for former-smoker QALYs and for smoker/former-smoker costs.
# Goal: inspect the distributions of differences in QALYs and costs.

n_sim <- 1000
n_cycles <- 240

# Parallel loop for simulations.
qaly_fosmo_probab <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- LVR <- LTx <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_fosmo$stage_1[1]
  for(j in 1:n_cycles){
    stage_1[j+1] <- stage_1[j]*(1 - Pb_prog1_int() - dat_fosmo$mort1[j])
  }
  
  # Stage 2 initialization and loop:
  stage_2[2] <- stage_1[1]*Pb_prog1_int()
  stage_2[3] <- stage_2[2]*(1 - Pb_prog2_int() - dat_fosmo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2]*Pb_prog1_int()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j]*(1 - Pb_prog2_int() - dat_fosmo$mort2[j] - Pb_p_sev2()*ae_mort2_fosmo[j]) + stage_1[j]*Pb_prog1_int()
  }
  
  # Stage 3 initialization and loop:
  stage_3[3] <- stage_2[2]*Pb_prog2_int()
  for(j in 3:n_cycles){
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3_int() - dat_fosmo$mort3[j] - Pb_p_sev3()*ae_mort3_fosmo[j]) + stage_2[j]*Pb_prog2_int()
  }
  
  # Stage 4 initialization and loop:
  stage_4[4] <- stage_3[4]*Pb_prog3_int()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j]*(1 - dat_fosmo$mort3[j] - Pb_p_sev4()*ae_mort3_fosmo[j]) + stage_3[j]*Pb_prog3_int()
  }
  
  # LVR and LTx calculations:
  LVR[4] <- stage_4[3] * Pb_p_4ps() * (1 - Pb_p_lvr_death())
  for(j in 4:n_cycles){
    LVR[j+1] <- LVR[j] * (1 - Pb_p_pspt(j) - dat_fosmo$mort3[j] - Pb_p_sevL()*ae_mort3_fosmo[j]) + stage_4[j] * Pb_p_4ps() * (1 - Pb_p_lvr_death())
  }
  
  for(j in 3:n_cycles){
    LTx[j+1] <- LTx[j] * (1 - dat_fosmo$mort4[j]) + stage_4[j] * Pb_p_4pt(j) * (1 - Pb_p_death_atltx()) + LVR[j] * Pb_p_pspt(j) * (1 - Pb_p_death_atltx())
  }
  
# Final QALY calculation.
  QALY <- ((stage_1*(Pb_u_1() + Pb_p_mild1()*d_mild1) +
            stage_2*(Pb_u_2() + Pb_p_mild2()*d_mild2 + Pb_p_mod2()*d_mod2 + Pb_p_sev2()*d_sev2)) +
           stage_3*(Pb_u_3() + Pb_p_mild3()*d_mild3 + Pb_p_mod3()*d_mod3 + Pb_p_sev3()*d_sev3) +
           stage_4*(Pb_u_4() + Pb_p_mild4()*d_mild4 + Pb_p_mod4()*d_mod4 + Pb_p_sev4()*d_sev4) +
           LVR*(Pb_u_lvr() + Pb_p_mildL()*d_mildL + Pb_p_modL()*d_modL + Pb_p_sevL()*d_sevL) +
           LTx*Pb_u_1tx()) * cycle_length/12
  
  QALY
}

dim(qaly_fosmo_probab) #241 1000, OK

for(i in 1:1000){qaly_fosmo_probab[1,i] <- 100000 * Pb_p_mild1() * c_mild1}

quantile(colSums(qaly_fosmo_probab), probs=c(0.025,0.975))
hist(colSums(qaly_fosmo_probab)) #loOKs OK.

hist(colSums(qaly_fosmo_probab) - colSums(qaly_smo_probab)) #loOKs OK.
dev.off()

qaly_fosmo_probab_discounted <- discount(qaly_fosmo_probab)

# Next: analogous calculations for costs_smo_probab and costs_fosmo_probab.
# Goal: inspect the distributions of differences in QALYs and costs.


n_sim <- 1000
n_cycles <- 240

costs_smo_probab <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- LVR <- LTx <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_smo$stage_1[1]
 for(i in 1:n_cycles){
   stage_1[i+1] <- stage_1[i]*(1 - Pb_prog1() - dat_smo$mort1[i])
 }
 # Stage 2 initialization and loop:
 
 stage_2[2] <- stage_1[1]*Pb_prog1()
 stage_2[3] <- stage_2[2]*(1 - Pb_prog2() - dat_smo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2]*Pb_prog1()
 for(i in 3:n_cycles){
   stage_2[i+1] <- stage_2[i]*(1 - Pb_prog2() - dat_smo$mort2[i] - Pb_p_sev2()*ae_mort2[i]) + stage_1[i]*Pb_prog1()
 }

  # Stage 3 initialization and loop:
 stage_3[3] <- stage_2[2]*Pb_prog2()
  for(j in 3:n_cycles){
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3() - dat_smo$mort3[j] - Pb_p_sev3()*ae_mort3[j]) + stage_2[j]*Pb_prog2()
  }
  
  # Stage 4 initialization and loop:
 stage_4[4] <- stage_3[4]*Pb_prog3()
 for(j in 3:n_cycles){
   stage_4[j+1] <- stage_4[j]*(1 - dat_smo$mort3[j] - Pb_p_sev3()*ae_mort3[j]) + stage_3[j]*Pb_prog3()
 }
 
  # LVR and LTx calculations:
 LVR[4] <- stage_4[3] * Pb_p_4ps()*(1 - Pb_p_lvr_death())
 
 for(i in 4:n_cycles){
   LVR[i+1] <- LVR[i] * (1 - Pb_p_pspt(i) - dat_smo$mort3[i] - Pb_p_sevL()*ae_mort3[i]) + stage_4[i]*Pb_p_4ps()*(1-Pb_p_lvr_death())
 }
 
 for(i in 3:n_cycles){
   #if (i>80){pb_p_4pt <- pb_p_pspt <- 0}  
   LTx[i+1] <- LTx[i] * (1 - dat_smo$mort4[i]) + stage_4[i] * Pb_p_4pt(i) * (1 - Pb_p_death_atltx()) + LVR[i] * Pb_p_pspt(i) * (1 - Pb_p_death_atltx())
 }
 
 costs <- (stage_1*(Pb_c_1() + Pb_p_mild1()*Pb_c_mild1()) + stage_2*(Pb_c_2() + Pb_p_mild2()*Pb_c_mild2() + Pb_p_mod2()*Pb_c_mod2() + Pb_p_sev2()*Pb_c_sev2()) + stage_3*(Pb_c_3() + Pb_p_mild3()*Pb_c_mild3() + Pb_p_mod3()*Pb_c_mod3() + Pb_p_sev3()*Pb_c_sev3()) + stage_4*(Pb_c_4() + Pb_p_mild4()*Pb_c_mild4() + Pb_p_mod4()*Pb_c_mod4() + Pb_p_sev4()*Pb_c_sev4()) + Pb_p_4ps()*Pb_c_atlvr() + LVR*(Pb_c_sev4() + Pb_p_mildL()*Pb_c_mild4() + Pb_p_modL()*Pb_c_mod4() + Pb_p_sevL() * Pb_c_sev4()) + (stage_4 * Pb_p_4pt(i) + LVR*Pb_p_pspt(i)) * Pb_c_atltx() + LTx * Pb_c_ltx())
 
 costs
}
costs_smo_probab

dim(costs_smo_probab) #OK soweit.

for(i in 1:1000){costs_smo_probab[1,i] <- 100000 * Pb_p_mild1() * Pb_c_mild1()}
costs_smo_probab[241,] <- 0  #Set the final value to 0 as in the Excel sheet.

quantile(colSums(costs_smo_probab), probs=c(0.025,0.975))
hist(colSums(costs_smo_probab)) #loOKs ....
dev.off()


#costs_smo_probab_discounted:
costs_smo_probab_discounted <- discount(costs_smo_probab)

# Discounted QALYs:


# Repeat the same calculation for costs_fosmo_probab, then calculate differences and inspect the distributions.
# Goal: inspect the distributions of differences in QALYs and costs.

n_sim <- 1000
n_cycles <- 240

costs_fosmo_probab <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- LVR <- LTx <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_fosmo$stage_1[1]
 for(i in 1:n_cycles){
 stage_1[i+1] <- stage_1[i]*(1 - Pb_prog1_int() - dat_fosmo$mort1[i])
 }
 # Stage 2 initialization and loop:
 
 stage_2[2] <- stage_1[1]*Pb_prog1_int()
 stage_2[3] <- stage_2[2]*(1 - Pb_prog2_int() - dat_fosmo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2]*Pb_prog1_int()
 for(i in 3:n_cycles){
  stage_2[i+1] <- stage_2[i]*(1 - Pb_prog2_int() - dat_fosmo$mort2[i] - Pb_p_sev2()*ae_mort2_fosmo[i]) + stage_1[i]*Pb_prog1_int()
 }

  # Stage 3 initialization and loop:
 stage_3[3] <- stage_2[2]*Pb_prog2_int()
  for(j in 3:n_cycles){
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3_int() - dat_fosmo$mort3[j] - Pb_p_sev3()*ae_mort3_fosmo[j]) + stage_2[j]*Pb_prog2_int()
  }
  
  # Stage 4 initialization and loop:
 stage_4[4] <- stage_3[4]*Pb_prog3_int()
 for(j in 3:n_cycles){
   stage_4[j+1] <- stage_4[j]*(1 - dat_fosmo$mort3[j] - Pb_p_sev3()*ae_mort3_fosmo[j]) + stage_3[j]*Pb_prog3_int()
 }
 
  # LVR and LTx calculations:
 LVR[4] <- stage_4[3] * Pb_p_4ps()*(1 - Pb_p_lvr_death())
 
 for(i in 4:n_cycles){
   LVR[i+1] <- LVR[i] * (1 - Pb_p_pspt(i) - dat_fosmo$mort3[i] - Pb_p_sevL()*ae_mort3_fosmo[i]) + stage_4[i]*Pb_p_4ps()*(1-Pb_p_lvr_death())
 }
 
 for(i in 3:n_cycles){
   #if (i>80){pb_p_4pt <- pb_p_pspt <- 0}  
   LTx[i+1] <- LTx[i] * (1 - dat_fosmo$mort4[i]) + stage_4[i] * Pb_p_4pt(i) * (1 - Pb_p_death_atltx()) + LVR[i] * Pb_p_pspt(i) * (1 - Pb_p_death_atltx())
 }
 
 costs <- (stage_1*(Pb_c_1() + Pb_p_mild1()*Pb_c_mild1()) + stage_2*(Pb_c_2() + Pb_p_mild2()*Pb_c_mild2() + Pb_p_mod2()*Pb_c_mod2() + Pb_p_sev2()*Pb_c_sev2()) + stage_3*(Pb_c_3() + Pb_p_mild3()*Pb_c_mild3() + Pb_p_mod3()*Pb_c_mod3() + Pb_p_sev3()*Pb_c_sev3()) + stage_4*(Pb_c_4() + Pb_p_mild4()*Pb_c_mild4() + Pb_p_mod4()*Pb_c_mod4() + Pb_p_sev4()*Pb_c_sev4()) + Pb_p_4ps()*Pb_c_atlvr() + LVR*(Pb_c_sev4() + Pb_p_mildL()*Pb_c_mild4() + Pb_p_modL()*Pb_c_mod4() + Pb_p_sevL() * Pb_c_sev4()) + (stage_4 * Pb_p_4pt(i) + LVR*Pb_p_pspt(i)) * Pb_c_atltx() + LTx * Pb_c_ltx())
 
 costs
}
costs_fosmo_probab
stopCluster(cl)

dim(costs_fosmo_probab) #OK soweit.

for(i in 1:1000){costs_fosmo_probab[1,i] <- 100000 * Pb_p_mild1() * Pb_c_mild1()}
costs_fosmo_probab[241,] <- 0  #Set the final value to 0 as in the Excel sheet.

summary(colSums(costs_fosmo_probab))
quantile(colSums(costs_fosmo_probab), probs=c(0.025,0.975))
hist(colSums(costs_fosmo_probab)) #looks OKs ....
dev.off()

#costs_fosmo_probab_discounted:
costs_fosmo_probab_discounted <- discount(costs_fosmo_probab)


quantile(colSums(costs_fosmo_probab), probs=c(0.025,0.975))
hist(colSums(costs_fosmo_probab)) #looks OK. Warning probably occurs because the numbers are very large for R. Dividing everything by 100 avoids the warning.

quantile(colSums(costs_fosmo_probab), probs=c(0.025,0.975))
hist(colSums(costs_fosmo_probab)) #looks OK.
hist(colSums(costs_smo_probab) - colSums(costs_fosmo_probab)) #looks OK.

# Next: plots of differences in QALYs versus differences in costs.
delta_qaly_probab <- qaly_smo_probab - qaly_fosmo_probab
colsums_delta_qaly_probab <- colSums(delta_qaly_probab)

delta_costs_probab <- costs_smo_probab - costs_fosmo_probab
colsums_delta_costs_probab <- colSums(delta_costs_probab)

plot(colsums_delta_qaly_probab, colsums_delta_costs_probab)
#Tada! Optional next step: repeat this after dividing by the number of simulated individuals.


## Probabilistic version of "direct costs, intervention (Int)":
total_dir_costs_intervention_probab_discounted <- colSums(costs_smo_probab_discounted)/100000 * (1-pb_eff_int[1:1000]) + colSums(costs_fosmo_probab_discounted)/100000 * pb_eff_int[1:1000] + pb_c_int[1:1000] #OK

## Probabilistic version of "direct costs, standard care":
total_dir_costs_standard_probab_discounted <- colSums(costs_smo_probab_discounted)/100000 * (1-pb_eff_uc[1:1000]) + colSums(costs_fosmo_probab_discounted)/100000 * pb_eff_uc[1:1000] + pb_c_uc[1:1000] #OK

#Overall ICER (delta costs / delta effects) for direct costs:
x <- (total_dir_costs_intervention_probab_discounted - total_dir_costs_standard_probab_discounted) / (qaly_smo_probab_discounted - qaly_fosmo_probab_discounted)
x[1,] <- NA
#colSums(x, na.rm=T)
plot(colSums(x, na.rm=T))
summary(colSums(x, na.rm=T))


######################
# For comparison: plot using the BCEA package:

library(BCEA)
costs <- as.matrix(cbind(total_dir_costs_standard_probab_discounted, total_dir_costs_intervention_probab_discounted))
effects <- as.matrix(cbind(colSums(qaly_fosmo_probab_discounted), colSums(qaly_smo_probab_discounted)))
m <- bcea(e=effects, c=costs, ref=1)
plot(m)
dev.off()
summary(m)
# Check whether both plots are plausible and where differences may arise.

#Overall ICER (delta costs / delta effects) for indirect costs:
## First calculate indirect costs probabilistically.
### First express the deterministic calculation as formulas.
(total_indir_costs_intervention - total_indir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #see around line 650

((sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_int) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int) #total_indir_costs_intervention

((sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc) #total_indir_costs_standard

(sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.22) #total_QALY_intervention

(sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.06) #total_QALY_standard

(((sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_int) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int) - ((sum(dat_smo$Kosten_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$Kosten_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc)) / ((sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_int) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_int) - (sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_uc) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_uc)) #OK.

# Now implement this probabilistically.
sample(colSums(costs_smo_probab_discounted), 1, replace=TRUE) / 100000 #replaces sum(dat_smo$costs_discounted) / 100000 #OK so far.

# For ind_costs_discounted, dat_smo$ind_costs and dat_fosmo$ind_costs are needed.
# rowSums(dat_smo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs")], na.rm=T)
# rowSums(dat_fosmo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs")], na.rm=T)
# This in turn requires the stage-specific indirect-cost calculations.
# dat_fosmo$stage_2 * (costs_stage2_stabil + p_mod2*costs_stage2_mod_AE + p_sev2*costs_stage2_schwere_AE)
#etc.
# It is unclear how to model costs_stage2_stabil etc. probabilistically.
# These parameters do not appear in the Excel probabilistic calculations.

##################
# Full model, including probabilistic sensitivity analyses
###################
