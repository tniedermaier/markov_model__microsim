
########################
# Base R implementation of a simplified Markov model for COPD. R implementation and optimization by Tobias Niedermaier. Original model design, parametrization, and implementation by Petra Menn. Supervisor: Rolf Holle.

#The model comprises 5 stages (GOLD 1, 2, 3, 4, + death as absorbing state), 2 groups (intervention + usual care) and has a cycle length of 3 months. The intervention (smoking cessation) slows down disease progression. Because more advanced stages are associated with higher costs and mortality and lower quality of life, the intervention (indirectly) also lowers costs and quality of life while lowering mortalitiy. During each cycle in stages 2-4, acute exacerbations (moderate or severe) may occur which reduce the utility and increase the costs and mortality in that cycle. Exacerbation risk differs between smokers and former smokers, by GOLD stage, and changes over time (according to cycle). Costs and utilities are discounted equally with a rate of 3% per year.
#We will hereafter start with the deterministic Markov model, followed by a probabilistic implementation.

rm(list=ls(all=TRUE))
options(scipen=999) #disables scientific notation
#Var.name        value     meaning
prog1         = 0.01403  #probability progression stage 1
prog2         = 0.01031  #probability progression stage 2
prog3         = 0.02298  #probability progression stage 3
prog1_int     = 0.00247  #probability progression stage 1, intervention group
prog2_int     = 0.003    #probability progression stage 2, intervention group
prog3_int     = 0.00686  #probability progression stage 3, intervention group

## Global switch for acute exacerbations:
# 1 = active (source-compatible simplified model; default);
# 0 = deactivated (simplification, corresponds to updated COSYCONET base case) ##

ae_switch <- 1

p_mod2        = 0.18376  #probability moderate exacerbation stage 2
d_mod2        = -0.01975 #change in utility moderate exacerbation stage 2
p_sev2        = 0.03854  #probability severe exacerbation stage 2
d_sev2        = -0.03075 #change in utility severe exacerbation stage 2

p_mod3        = 0.26767  #probability moderate exacerbation stage 3
d_mod3        = -0.01875 #change in utility moderate exacerbation stage 3
p_sev3        = 0.03774  #probability severe exacerbation stage 3
d_sev3        = -0.02919 #change in utility severe exacerbation stage 3

p_mod4        = 0.28057  #probability moderate exacerbation stage 4
d_mod4        = -0.01625 #change in utility moderate exacerbation stage 4
p_sev4        = 0.04855  #probability severe exacerbation stage 4
d_sev4        = -0.03502 #change in utility severe exacerbation stage 4

if (ae_switch == 0){
 p_mod2 = p_sev2 = p_mod3 = p_sev3 = p_mod4 = p_sev4 = 0
}

u_1           = 0.84     #utility stage 1
u_2           = 0.79     #utility stage 2
u_3           = 0.75     #utility stage 3
u_4           = 0.65     #utility stage 4

c_1           = 102.64   #costs in stable stage 1

c_2           = 184.69   #costs in stable stage 2
c_mod2        = 39.86    #costs mod ae, stage 2
c_sev2        = 2923.94  #costs sev ae, stage 2
c_3           = 366.57   #costs in stable stage 3
c_mod3        = 52.90    #costs mod ae, stage 3
c_sev3        = 3109.05  #costs sev ae, stage 3
c_4           = 431.35   #costs in stable stage 4
c_mod4        = 51.52    #costs mod ae, stage 4
c_sev4        = 3322.59  #costs sev ae, stage 4
c_int         = 595.92   #Interventionskosten
c_uc          = 0        #initial costs usual care

cycle_length  = 3        #[months]
disc_e        = 0.03     #discount rate
eff_int       = 0.22
eff_uc        = 0.06

##Begin with smokers (dat_smo), then repeat with former smokers (dat_fosmo).
N <- 241 #Number of cycles.
dat_smo <- data.frame(matrix(nrow=N, ncol=5))
colnames(dat_smo) <- c("stage_1","stage_2","stage_3","stage_4","death")
dat_smo[] <- 0
dat_smo[1,"stage_1"] <- 100000

##Define mortality rates per cycle (age), according to stage. The same age-specific mortality rates are assumed for stages 3 and 4.
dat_smo[,c("mort1","mort2","mort3")] <- NA
dat_smo[,c("mort1","mort2","mort3")][1:20,] <- rep(c(0.00122,0.0021,0.00497), each=20)
dat_smo[,c("mort1","mort2","mort3")][21:40,] <- rep(c(0.0021,0.0036,0.00852), each=20)
dat_smo[,c("mort1","mort2","mort3")][41:60,] <- rep(c(0.00302,0.00518,0.01225), each=20)
dat_smo[,c("mort1","mort2","mort3")][61:80,] <- rep(c(0.0047,0.00804,0.01899), each=20)
dat_smo[,c("mort1","mort2","mort3")][81:100,] <- rep(c(0.00805,0.01377,0.03238), each=20)
dat_smo[,c("mort1","mort2","mort3")][101:120,] <- rep(c(0.01358,0.02318,0.05417), each=20)
dat_smo[,c("mort1","mort2","mort3")][121:140,] <- rep(c(0.0197,0.03354,0.07783), each=20)
dat_smo[,c("mort1","mort2","mort3")][141:160,] <- rep(c(0.034,0.05758,0.13138), each=20)
dat_smo[,c("mort1","mort2","mort3")][161:180,] <- rep(c(0.04729,0.07969,0.179), each=20)
dat_smo[,c("mort1","mort2","mort3")][181:241,] <- rep(c(0.0884,0.14672,0.31397), each=61)

dat_smo$stage_1[-1] <- (100000 * cumprod(1 - prog1 - dat_smo$mort1))[1:240]

#Stage 2, smokers:
dat_smo$stage_2 <- NA
dat_smo$stage_2[2] <- dat_smo[1,"stage_1"]*prog1 #1st cycle with N(stage 2)>0 as starting point. 1st value remains NA.
dat_smo$stage_2[3] <- dat_smo$stage_2[2]*(1 - prog2 - dat_smo$mort2[1] - p_sev2*0.13819) + dat_smo$stage_1[2]*prog1

#Excess mortality in cycles with acute exacerbations in stage 2, separately for smokers and former smokers:
ae_mort2 <- ae_mort2_fosmo <- vector(mode="numeric", length=241)
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

for(i in 3:240){dat_smo$stage_2[i+1] <- dat_smo$stage_2[i] * (1 - prog2 - dat_smo$mort2[i] - p_sev2 * ae_mort2[i]) + dat_smo$stage_1[i] * prog1}
dat_smo$stage_2

#Stage 3, smokers:
dat_smo$stage_3 <- NA
dat_smo$stage_3[2] <- 0 #by definition
dat_smo$stage_3[3] <- dat_smo[2,"stage_2"]*prog2 #1st cycle with N(stage 3)>0 as starting point.

#Excess mortality in cycles with acute exacerbations in stage 3, separately for smokers and former smokers:
ae_mort3 <- ae_mort3_fosmo <- vector(mode="numeric", length=241)
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

for(i in 3:240){dat_smo$stage_3[i+1] <- dat_smo$stage_3[i] * (1 - prog3 - dat_smo$mort3[i] - p_sev3 * ae_mort3[i]) + dat_smo$stage_2[i] * prog2}

#Stage 4, smokers:
dat_smo$stage_4 <- NA
#dat_smo$stage_4[1] <- NA #by definition
dat_smo$stage_4[2:3] <- 0
dat_smo$stage_4[4] <- dat_smo[4,"stage_3"]*prog3 #1st cycle with N(stage 4)>0 as starting point.

for(i in 3:240){dat_smo$stage_4[i+1] <- dat_smo$stage_4[i] * (1 - dat_smo$mort3[i] - p_sev4 * ae_mort3[i]) + dat_smo$stage_3[i] * prog3}
dat_smo$stage_4


#death, smokers:
dat_smo$death <- NA
dat_smo$death <- 100000 - rowSums(dat_smo[,c("stage_1","stage_2","stage_3","stage_4")])

#Life-years gained (LYG), smokers, not discounted:
dat_smo$LYG <- rowSums(dat_smo[,c("stage_1","stage_2","stage_3","stage_4")]) * cycle_length / 12
dat_smo$LYG[c(1,241)] <- 0

Cycle <- 0:240
discount_factor <- 1/(((1+disc_e)^(cycle_length/12))^Cycle)
(dat_smo$LYG_discounted <- dat_smo$LYG * discount_factor)

#Quality-adjusted life-years (QALYs): We take the number of people per stage and cycle, multiply it with the utility in each stage and sum it up over all cycles. Note that d_mod2-4 and d_sev2-4 are negative. Thus, the utility is diminished if moderate or severe exacerbations occur.
qaly <- function(data, i){
 (data[i,"stage_1"]*(u_1) + data[i,"stage_2"] *(u_2 + p_mod2*d_mod2 + p_sev2*d_sev2) + data[i,"stage_3"]*(u_3 + p_mod3*d_mod3 + p_sev3*d_sev3) + data[i,"stage_4"]*(u_4 + p_mod4*d_mod4 + p_sev4*d_sev4))*cycle_length/12
}

dat_smo$QALY <- qaly(dat_smo, i=1:241)
dat_smo$QALY[1] <- 0
sum(dat_smo$QALY) #1775641

(dat_smo$QALY_discounted <- dat_smo$QALY * discount_factor)
sum(dat_smo$QALY_discounted) #1248523


#Costs:

costs_smo_new2 <- Vectorize(function(data, i){
 data$stage_1[i]*(c_1) + data$stage_2[i]*(c_2 + p_mod2*c_mod2 + p_sev2*c_sev2) + data$stage_3[i]*(c_3 + p_mod3*c_mod3 + p_sev3*c_sev3) + data$stage_4[i]*(c_4 + p_mod4*c_mod4 + p_sev4*c_sev4)
}, vectorize.args="i")
costs_new2 <- costs_smo_new2(data=dat_smo, i=1:241)
#Starting value is 0:
costs_new2[1] <- 0

dat_smo$costs <- costs_new2

sum(dat_smo$costs) #1757475693.

#discounted costs:
dat_smo$costs_discounted <- dat_smo$costs * discount_factor
sum(dat_smo$costs_discounted) #1164532502.


#dat_fosmo: This is practically identical to the codes above, but using different parameters because the progression and mortality rates are lower in former smokers compared to smokers.

N <- 241 #Number of cycles.
dat_fosmo <- data.frame(matrix(nrow=N, ncol=5))
colnames(dat_fosmo) <- c("stage_1","stage_2","stage_3","stage_4","death")
dat_fosmo[] <- 0
dat_fosmo[1,"stage_1"] <- 100000 #start with 100000 people in stage 1 (and 0 in the other stages)

#Define mortality in the different cycles for the different stages for former smokers:
dat_fosmo[,c("mort1","mort2","mort3")] <- NA
dat_fosmo[,c("mort1","mort2","mort3")][1:20,] <- rep(c(0.0008,0.00136,0.00323), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][21:40,] <- rep(c(0.00137,0.00234,0.00555), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][41:60,] <- rep(c(0.00251,0.0043,0.01018), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][61:80,] <- rep(c(0.0039,0.00668,0.01578), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][81:100,] <- rep(c(0.00504,0.00863,0.02037), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][101:120,] <- rep(c(0.00851,0.01455,0.03421), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][121:140,] <- rep(c(0.01318,0.02249,0.05258), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][141:160,] <- rep(c(0.0228,0.03877,0.08963), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][161:180,] <- rep(c(0.03957,0.06686,0.15156), each=20)
dat_fosmo[,c("mort1","mort2","mort3")][181:241,] <- rep(c(0.07423,0.12386,0.2695), each=61)

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

dat_fosmo$stage_1

#Stage 2, former smokers:
dat_fosmo$stage_2 <- dat_fosmo$stage_3 <- dat_fosmo$stage_4 <- dat_fosmo$death <- NA
dat_fosmo$stage_2[2] <- dat_fosmo[1,"stage_1"]*prog1_int #1st cycle with N(stage 2)>0 as starting point.
dat_fosmo$stage_2[3] <- dat_fosmo$stage_2[2]*(1 - prog2_int - dat_fosmo$mort2[1] - p_sev2*0.13883) + dat_fosmo$stage_1[2]*prog1_int

for(i in 3:240){dat_fosmo$stage_2[i+1] <- dat_fosmo$stage_2[i] * (1 - prog2_int - dat_fosmo$mort2[i] - p_sev2 * ae_mort2_fosmo[i]) + dat_fosmo$stage_1[i] * prog1_int}
dat_fosmo$stage_2


#Stage 3, former smokers:

#=AG5*(1-prog3_int-Z5-p_sev3*AD5)+AF5*prog2_int
dat_fosmo$stage_3[2] <- 0 #dat_fosmo$stage_3[1] * prog2_int
dat_fosmo$stage_3[3] <- dat_fosmo$stage_2[2] * prog2_int

for(i in 3:240){dat_fosmo$stage_3[i+1] <- dat_fosmo$stage_3[i] * (1 - prog3_int - dat_fosmo$mort3[i] - p_sev3 * ae_mort3_fosmo[i]) + dat_fosmo$stage_2[i] * prog2_int}
dat_fosmo$stage_3


#Stage 4, former smokers:

dat_fosmo$stage_4[2] <- 0
dat_fosmo$stage_4[3] <- dat_fosmo$stage_3[2] * prog3_int

for(i in 3:240){dat_fosmo$stage_4[i+1] <- dat_fosmo$stage_4[i] * (1 - dat_fosmo$mort3[i] - p_sev4 * ae_mort3_fosmo[i]) + dat_fosmo$stage_3[i] * prog3_int}
dat_fosmo$stage_4

#death, former smokers:
dat_fosmo$death <- 100000 - rowSums(dat_fosmo[,c("stage_1","stage_2","stage_3","stage_4")])


dat_fosmo$QALY <- NA
(dat_fosmo$QALY[2:241] <- qaly(data=dat_fosmo, i=2:241))

dat_fosmo$QALY #+initial value which is calculated differently from the others:
dat_fosmo$QALY[1] <- 0
sum(dat_fosmo$QALY) #2520182
(dat_fosmo$QALY_discounted <- dat_fosmo$QALY * discount_factor)
sum(dat_fosmo$QALY_discounted) #1589728

#LYG, fosmo, not discounted:
dat_fosmo$LYG <- rowSums(dat_fosmo[,c("stage_1","stage_2","stage_3","stage_4")]) * cycle_length / 12
dat_fosmo$LYG[c(1,241)] <- 0
sum(dat_fosmo$LYG) #3020696

(dat_fosmo$LYG_discounted <- dat_fosmo$LYG * discount_factor)

#costs, fosmo:
dat_fosmo$costs <- costs_smo_new2(data=dat_fosmo, i=1:241)
#Provide value for 1st cycle:
dat_fosmo$costs[1] <- 0 #=0 because c_mild1 is one of the factors for the initial value, which is now 0.

sum(dat_fosmo$costs) #1494553273. Original model: 1504653567.

#discounted costs:
dat_fosmo$costs_discounted <- dat_fosmo$costs * discount_factor
sum(dat_fosmo$costs_discounted) #923630944. Original model: 929724839.


#Indirect costs: Note that GOLD stage 1 is not assumed to cause inability to work. Furthermore, people retire (by the latest) at age 65 years, with increasing proportions of inability to work from GOLD 2 (15%) over GOLD 3 (35%) to GOLD 4 (37%).

#Vectors with the costs per stage and exacerbations. (To be used for both, dat_smo and dat_fosmo. Thus, it is better to keep them as vectors and not attach them to dat_smo or dat_fosmo as data.frame.):
costs_stage2_stable <- c(0, rep(112.50, 19), rep(376.37, 20), rep(667.24, 20), rep(237.22, 20), rep(0, 161)) #ff: columns AT - BF in copd08.xlsx.
costs_stage3_stable <- c(0, rep(268.06, 19), rep(896.81, 20), rep(1589.9, 20), rep(565.24, 20), rep(0, 161))
costs_stage4_stable <- c(0, rep(281.24, 19), rep(940.92, 20), rep(1668.1, 20), rep(593.04, 20), rep(0, 161))

costs_stage2_mod_AE <- c(rep(491.94, 20), rep(446.16, 20), rep(364.35, 20), rep(115.30, 20), rep(0, 161))
costs_stage3_mod_AE <- c(rep(550.51, 20), rep(499.27, 20), rep(407.73, 20), rep(129.03, 20), rep(0, 161))
costs_stage4_mod_AE <- c(rep(491.94, 20), rep(446.16, 20), rep(364.35, 20), rep(115.30, 20), rep(0, 161))
costs_stage2_severe_AE <- c(rep(367, 20), rep(332.85, 20), rep(271.82, 20), rep(86.02, 20), rep(0, 161))
costs_stage3_severe_AE <- c(rep(663.73, 20), rep(601.96, 20), rep(491.58, 20), rep(155.56, 20), rep(0,161))
costs_stage4_severe_AE <- c(rep(648.12, 20), rep(587.79, 20), rep(480.02, 20), rep(151.9, 20), rep(0,161))

#indirect costs, smokers, not discounted:
dat_smo$stage_1_ind_costs <- 0 #= dat_smo$stage_1  * p_mild1 * x, but p_mild1=0

#indirect costs, smokers stage 2:
(dat_smo$stage_2_ind_costs <- dat_smo$stage_2 * (costs_stage2_stable + p_mod2*costs_stage2_mod_AE + p_sev2*costs_stage2_severe_AE)) #OK.

#indirect costs, smokers stage 3:
(dat_smo$stage_3_ind_costs <- dat_smo$stage_3 * (costs_stage3_stable + p_mod3*costs_stage3_mod_AE + p_sev3*costs_stage3_severe_AE)) #OK.

#indirect costs, smokers stage 4:
(dat_smo$stage_4_ind_costs <- dat_smo$stage_4 * (costs_stage4_stable + p_mod4*costs_stage4_mod_AE + p_sev4*costs_stage4_severe_AE)) #OK.

(dat_smo$ind_costs <- rowSums(dat_smo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs")], na.rm=T))

#indirect costs, smokers, discounted:
(dat_smo$ind_costs_discounted <- dat_smo$ind_costs * discount_factor)

#indirect costs, former smokers, stage 1:
(dat_fosmo$stage_1_ind_costs <- 0)

#indirect costs, former smokers, stage 2:
(dat_fosmo$stage_2_ind_costs <- dat_fosmo$stage_2 * (costs_stage2_stable + p_mod2*costs_stage2_mod_AE + p_sev2*costs_stage2_severe_AE))

#indirect costs, former smokers, stage 3:
(dat_fosmo$stage_3_ind_costs <- dat_fosmo$stage_3 * (costs_stage3_stable + p_mod3*costs_stage3_mod_AE + p_sev3*costs_stage3_severe_AE))

#indirect costs, former smokers, stage 4:
(dat_fosmo$stage_4_ind_costs <- dat_fosmo$stage_4 * (costs_stage4_stable + p_mod4*costs_stage4_mod_AE + p_sev4*costs_stage4_severe_AE))


##ICER calculations:
(dat_smo$effects <- dat_smo$QALY_discounted / 100000)

#former smokers:
(dat_fosmo$effects <- dat_fosmo$QALY_discounted / 100000)

#Int:
cycle_eff_int <- NA
(cycle_eff_int[1] <- dat_smo$effects[1] * (1-eff_int) + dat_fosmo$effects[1] * eff_int)
(cycle_eff_int[2] <- cycle_eff_int[1] + dat_smo$effects[2]*(1-eff_int) + dat_fosmo$effects[2]*eff_int)

for(i in 2:240){cycle_eff_int[i+1] <- cycle_eff_int[i] + dat_smo$effects[i+1]*(1-eff_int) + dat_fosmo$effects[i+1]*eff_int}
cycle_eff_int #This is the cumulative number of QALYs that are generated over time in the intervention group.

cycle_eff_st <- NA
cycle_eff_st[1] <- dat_smo$effects[1]*(1-eff_uc) + dat_fosmo$effects[1]*eff_uc
cycle_eff_st[2] <- cycle_eff_st[1] + dat_smo$effects[2]*(1-eff_uc) + dat_fosmo$effects[2]*eff_uc

for(i in 2:240){cycle_eff_st[i+1] <- cycle_eff_st[i] + dat_smo$effects[i+1]*(1-eff_uc) + dat_fosmo$effects[i+1]*eff_uc}
cycle_eff_st #This is the cumulative number of QALYs that are generated over time in the standard treatment group.

(incr_eff <- cycle_eff_int - cycle_eff_st) #Incremental effectiveness = N(QALYs with intervention) - N(QALYs with standard treatment)
(dat_fosmo$ind_costs <- rowSums(dat_fosmo[,c("stage_1_ind_costs","stage_2_ind_costs","stage_3_ind_costs","stage_4_ind_costs")], na.rm=TRUE))
(dat_fosmo$ind_costs_discounted <- dat_fosmo$ind_costs * discount_factor) 
(dat_smo$IKER <- (dat_smo$costs_discounted + dat_smo$ind_costs_discounted)/100000)
(dat_fosmo$IKER <- (dat_fosmo$costs_discounted + dat_fosmo$ind_costs_discounted)/100000)

# Costs for the Intervention (?) ["Int"]:
c_int + dat_smo$IKER * (1-eff_int) + dat_fosmo$IKER*eff_int
dat_fosmo$Costs_Int[1] <- (c_int + dat_smo$IKER * (1-eff_int) + dat_fosmo$IKER*eff_int)[1]

dat_fosmo$Costs_Int[2] <- dat_fosmo$Costs_Int[1] + dat_smo$IKER[2]*(1-eff_int) + dat_fosmo$IKER[2]*eff_int #ok.
#For the remaining cycles:
for(i in 2:240){dat_fosmo$Costs_Int[i+1] <- dat_fosmo$Costs_Int[i] + dat_smo$IKER[i+1]*(1-eff_int) + dat_fosmo$IKER[i+1]*eff_int}
dat_fosmo$Costs_Int #ok.

# Costs for "standard care"/"usual care" ["ST"]:
dat_fosmo$ST[1] <- dat_smo$IKER[1]*(1-eff_uc) + dat_fosmo$IKER[1]*eff_uc #CG5*(1-eff_uc)+CH5*eff_uc
dat_fosmo$ST[2] <- dat_fosmo$ST[1] + dat_smo$IKER[2]*(1-eff_uc) + dat_fosmo$IKER[2]*eff_uc

for(i in 2:240){dat_fosmo$ST[i+1] <- dat_fosmo$ST[i] + dat_smo$IKER[i+1]*(1-eff_uc) + dat_fosmo$IKER[i+1]*eff_uc}
dat_fosmo$ST #ok.

dat_fosmo$inkr_Kosten <- dat_fosmo$Costs_Int - dat_fosmo$ST
(ICER_Menn <- dat_fosmo$inkr_Kosten / incr_eff)
plot(ICER_Menn, type="l", ylim=c(-20000,50000), xlab="Zyklus", ylab="IKER (€/QALY)", main="ICER over time", lwd=2, xaxt="n", xaxs="i", yaxs="i")
axis(1, seq(0,220,20))
abline(h=seq(-10000, 50000, 10000)) 
dev.off()


#Stage distribution over time in smokers and former smokers:

layout(matrix(c(1,2), ncol=2, byrow=TRUE))
#layout.show(n=2)
plot(dat_smo$stage_1, type="l", col="green", xlab="3-month cycle", main="Smokers")
lines(dat_smo$stage_2, type="l", col="blue")
lines(dat_smo$stage_3, type="l", col="orange")
lines(dat_smo$stage_4, type="l", col="brown")
lines(dat_smo$death, type="l", col="black")
legend(150,60000, legend=c("Stage I","Stage II","Stage III","Stage IV","Dead"), lty=1, col=c("green","blue","orange","brown","black"), bty="n")

plot(dat_fosmo$stage_1, type="l", col="green", xlab="3-month cycle", main="Former smokers")
lines(dat_fosmo$stage_2, type="l", col="blue")
lines(dat_fosmo$stage_3, type="l", col="orange")
lines(dat_fosmo$stage_4, type="l", col="brown")
lines(dat_fosmo$death, type="l", col="black")
legend(150,60000, legend=c("Stage I","Stage II","Stage III","Stage IV","Dead"), lty=1, col=c("green","blue","orange","brown","black"), bty="n")

mtext(expression(paste(bold("COPD stage distribution over time"))), side=3, line=-1, outer=TRUE)


dat_smo2 <- dat_smo[1:120,]
plot(dat_smo2$stage_1, type="l", col="green", xlab="3-month cycle", main="Stages over time (smokers)", ylim=c(0,100000), ylab="Participants")
lines(dat_smo2$stage_2, type="l", col="blue")
lines(dat_smo2$stage_3, type="l", col="orange")
lines(dat_smo2$stage_4, type="l", col="brown")
lines(dat_smo2$death, type="l", col="black")
legend(90,50000, legend=c("Stage I","Stage II","Stage III","Stage IV","Dead"), lty=1, col=c("green","blue","orange","brown","black"), bty="n")


##########################################
# Main results:


### smokers ###

## total LYG, "Standard":
sum(dat_smo$LYG) #2186346
## LYG per patient, "Standard":
sum(dat_smo$LYG) / 100000 #21.86346
## total QALYs, "Standard":
sum(dat_smo$QALY) #1775641
#per patient:
sum(dat_smo$QALY) / 100000 #17.75641
## total costs, "Standard":
sum(dat_smo$costs) #1757475693.
#per patient:
sum(dat_smo$costs) / 100000 #17574.76
## total LYG, discounted:
sum(dat_smo$LYG_discounted) #1531264
#per patient:
sum(dat_smo$LYG_discounted) / 100000 #15.31264
## total QALYs, discounted:
sum(dat_smo$QALY_discounted) #1248523
#per patient:
sum(dat_smo$QALY_discounted) / 100000 #12.48523
## total costs, discounted:
sum(dat_smo$costs_discounted) #1164532502
#per patient:
sum(dat_smo$costs_discounted) / 100000 #11645.33
## indirect costs, "Standard":
sum(dat_smo$ind_costs) #1338300535
#per patient:
sum(dat_smo$ind_costs) / 100000 #13383.01
## indirect costs, discounted:
sum(dat_smo$ind_costs_discounted) #949410662
#per patient:
sum(dat_smo$ind_costs_discounted) / 100000 #9494.107


### Former smokers ###

## total LYG, "Standard":
sum(dat_fosmo$LYG, na.rm=T) #3020696
## LYG per patient, "Standard":
sum(dat_fosmo$LYG, na.rm=T) / 100000 #30.20696
## total QALYs, "Standard":
sum(dat_fosmo$QALY) #2520182
#per patient:
sum(dat_fosmo$QALY) / 100000 #25.20182
## total costs, "Standard":
sum(dat_fosmo$costs) #1500596409.
#per patient:
sum(dat_fosmo$costs) / 100000 #15005.96
## total LYG, discounted:
sum(dat_fosmo$LYG_discounted, na.rm=T) #1903963
#per patient:
sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 #19.03963
## total QALYs, discounted:
sum(dat_fosmo$QALY_discounted) #1589728
#per patient:
sum(dat_fosmo$QALY_discounted) / 100000 #15.89728
## total costs, discounted:
sum(dat_fosmo$costs_discounted) #923630944.
#per patient:
sum(dat_fosmo$costs_discounted) / 100000 #9236.309
## indirect costs, "Standard":
sum(dat_fosmo$ind_costs) #278945553.
#per patient:
sum(dat_fosmo$ind_costs) / 100000 #2789.456
## indirect costs, discounted:
sum(dat_fosmo$ind_costs_discounted) #196146845.
#per patient:
sum(dat_fosmo$ind_costs_discounted) / 100000 #1961.468

## "Light-blue boxes" with summary results:

#direct costs, standard care (ST):
#=V$248*(1-eff_uc)+AR$248*eff_uc+c_uc
(total_dir_costs_standard <- sum(dat_smo$costs_discounted)/100000 * (1-eff_uc) + sum(dat_fosmo$costs_discounted)/100000 * eff_uc + c_uc) #11500.78.

#direct costs, Intervention (Int):
#=V$248*(1-eff_int)+AR$248*eff_int+c_int
(total_dir_costs_intervention <- sum(dat_smo$costs_discounted)/100000 * (1-eff_int) + sum(dat_fosmo$costs_discounted)/100000 * eff_int + c_int) #11711.26.

#LYG, standard care (ST):
#=T$248*(1-Werte!$C$3)+AP$248*Werte!$C$3
sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.06 #15.53626

#LYG, Intervention:
#=T$248*(1-Werte!$C$2)+AP$248*Werte!$C$2
sum(dat_smo$LYG_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$LYG_discounted, na.rm=T) / 100000 * 0.22 #16.13258

#QALY, standard care (ST):
#=U$248*(1-Werte!$C$3)+AQ$248*Werte!$C$3
(total_QALY_standard <- sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.06) #12.68995

#QALY, Intervention:
#=U$248*(1-Werte!$C$2)+AQ$248*Werte!$C$2
(total_QALY_intervention <- sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.22) #13.23588

#indirect costs, standard care (ST):
#=(V$248+BO$248)*(1-eff_uc)+(AR$248+BX$248)*eff_uc+c_uc
(total_indir_costs_standard <- (sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc) #20542.93.

#indirect costs, Intervention (Int):
(total_indir_costs_intervention <- (sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_int) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int) #19548.19.

#Overall ICER (delta costs / delta effects) for direct costs:
(total_dir_costs_intervention - total_dir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #385.5405.

#Overall ICER (delta costs / delta effects) for indirect costs:
(total_indir_costs_intervention - total_indir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #-1822.115.

###########################
# Consistency checks:

#1.) Sum of all states = 100000 in all cycles?
rowSums(dat_smo[,c("stage_1", "stage_2", "stage_3", "stage_4", "death")], na.rm=TRUE) #OK 

#Note: This is trivial/tautological because the number of deaths per cycle was calculated as difference between the initial number of subjects and the number of subjects that transitioned to other states. This is a simplified implementation; the original version calculated the number of deaths from the case numbers per state and the mortalitiy rates with and without AEs. 

alive_cols <- c("stage_1", "stage_2", "stage_3", "stage_4")

dat_smo$death_resid <- 100000 - rowSums(dat_smo[, alive_cols], na.rm = TRUE)

dat_smo$death_cycle_formula <- 0
dat_smo$death_formula <- 0

# Transition from row 1 to row 2:
# In stage_2[2], progressors from stage_1 are added,
# but stage_2 itself is not yet subject to mortality in this step.
dat_smo$death_cycle_formula[2] <-
  dat_smo$stage_1[1] * dat_smo$mort1[1]

dat_smo$death_formula[2] <-
  dat_smo$death_formula[1] + dat_smo$death_cycle_formula[2]


# Transition from row 2 to row 3:
# Special case due to your formula:
# stage_2[3] <- stage_2[2] * (1 - prog2 - mort2[1] - p_sev2 * 0.13819) + ...
dat_smo$death_cycle_formula[3] <-
  dat_smo$stage_1[2] * dat_smo$mort1[2] +
  dat_smo$stage_2[2] * (dat_smo$mort2[1] + p_sev2 * 0.13819)

dat_smo$death_formula[3] <-
  dat_smo$death_formula[2] + dat_smo$death_cycle_formula[3]


# From row 3 to row 4 onward: regular recursion
for (i in 3:(N - 1)) {
  
  dat_smo$death_cycle_formula[i + 1] <-
    dat_smo$stage_1[i] * dat_smo$mort1[i] +
    dat_smo$stage_2[i] * (dat_smo$mort2[i] + p_sev2 * ae_mort2[i]) +
    dat_smo$stage_3[i] * (dat_smo$mort3[i] + p_sev3 * ae_mort3[i]) +
    dat_smo$stage_4[i] * (dat_smo$mort3[i] + p_sev4 * ae_mort3[i])
  
  dat_smo$death_formula[i + 1] <-
    dat_smo$death_formula[i] + dat_smo$death_cycle_formula[i + 1]
}

# Check: should only be numerical noise
summary(dat_smo$death_resid - dat_smo$death_formula)
max(abs(dat_smo$death_resid - dat_smo$death_formula), na.rm = TRUE) ########## OK ##########


#Former smokers:
rowSums(dat_fosmo[,c("stage_1", "stage_2", "stage_3", "stage_4", "death")], na.rm=TRUE) #OK 

start_pop_fosmo <- dat_fosmo$stage_1[1]

dat_fosmo$death_resid <- start_pop_fosmo - rowSums(
  dat_fosmo[, alive_cols],
  na.rm = TRUE
)

dat_fosmo$death_cycle_formula <- 0
dat_fosmo$death_formula <- 0


# Transition from row 1 to row 2:
# In stage_2[2], progressors from stage_1 are added,
# but stage_2 itself is not yet subject to mortality in this step.
dat_fosmo$death_cycle_formula[2] <-
  dat_fosmo$stage_1[1] * dat_fosmo$mort1[1]

dat_fosmo$death_formula[2] <-
  dat_fosmo$death_formula[1] + dat_fosmo$death_cycle_formula[2]


# Transition from row 2 to row 3:
# Special case due to the formula used for stage_2[3]:
# stage_2[3] <- stage_2[2] * 
#   (1 - prog2_int - mort2[1] - p_sev2 * 0.13883) + ...
#
# IMPORTANT:
# Use 0.13883 here, not 0.13819.
# 0.13819 is the smoker value; 0.13883 is the former-smoker value.
dat_fosmo$death_cycle_formula[3] <-
  dat_fosmo$stage_1[2] * dat_fosmo$mort1[2] +
  dat_fosmo$stage_2[2] * (dat_fosmo$mort2[1] + p_sev2 * 0.13883)

dat_fosmo$death_formula[3] <-
  dat_fosmo$death_formula[2] + dat_fosmo$death_cycle_formula[3]


# From row 3 to row 4 onward: regular recursion
for (i in 3:(N - 1)) {
  
  dat_fosmo$death_cycle_formula[i + 1] <-
    dat_fosmo$stage_1[i] * dat_fosmo$mort1[i] +
    dat_fosmo$stage_2[i] * (dat_fosmo$mort2[i] + p_sev2 * ae_mort2_fosmo[i]) +
    dat_fosmo$stage_3[i] * (dat_fosmo$mort3[i] + p_sev3 * ae_mort3_fosmo[i]) +
    dat_fosmo$stage_4[i] * (dat_fosmo$mort3[i] + p_sev4 * ae_mort3_fosmo[i])
  
  dat_fosmo$death_formula[i + 1] <-
    dat_fosmo$death_formula[i] + dat_fosmo$death_cycle_formula[i + 1]
}


# Check: should only be numerical noise
summary(dat_fosmo$death_resid - dat_fosmo$death_formula)
max(abs(dat_fosmo$death_resid - dat_fosmo$death_formula), na.rm = TRUE) ########## OK ##########

#We just demonstrated that the number of individuals is 100000 for all cycles when summing up the states, as it should be.

# Continue with consistency checks:
#2.) Number of subjects non-negative in all states?

summary(sign(dat_smo[,c("stage_1","stage_2","stage_3","stage_4","death")])) ########## OK ##########
summary(sign(dat_fosmo[,c("stage_1","stage_2","stage_3","stage_4","death")])) ########## OK ##########
#If there were any negative values, the minimum in the columns would be -1, which is not the case.


#3.) Number of deaths is (weakly) monotonically increasing?
for(i in 1:length(dat_smo$death)){print(dat_smo$death[i+1] > dat_smo$death[i])} ########## OK ##########
for(i in 1:length(dat_fosmo$death)){print(dat_fosmo$death[i+1] > dat_fosmo$death[i])} ########## OK ##########





###########################
# probabilistic sensitivity analysis:

library(doParallel) #This is the only dependency in the codes. Parallel computing drastically speeds up the probabilistic sensitivity analysis on modern computers with several CPU cores. If you really have only 1 logical CPU core, you can skip this line and also the registerDoParallel(cl) function call later.

gamma_params <- function(mean, sd){
 k <- mean^2/sd^2
 phi_xls <- mean/k
 phi_tre <- mean/sd^2
 res <- list(k, phi_xls, phi_tre)
 names(res) <- c("alpha","beta_xls", "beta_tre")
 return(res)
}
gamma_params(mean=102.64, sd=20.00)

#pb_prog1:
(alpha_pb_prog1 <- prog1^2*(1-prog1)/0.00053^2-prog1) #690.9052
calc_alpha <- function(mu, se){mu^2*(1-mu)/se^2-mu}
calc_alpha(mu=prog1, se=0.00053) #ok

#Calculate beta analogously in a function:
calc_beta <- function(x, se){x*(1-x)^2/se^2-(1-x)}
calc_beta(x=prog1, se=0.00053) #ok
(beta_pb_prog1 <- calc_beta(x=prog1, se=0.00053))

set.seed(1) #The following functions [Pb_...()] have no arguments. They return a value from a beta distribution (for progression rates and utilities) or gamma distribution (for costs) every time they are called, using the defined distribution parameters. The seed ensures that the results are replicable.
Pb_prog1 <- function(){rbeta(1, calc_alpha(mu=prog1, se=0.00053), calc_beta(x=prog1, se=0.00053))} #replaces pb_prog1
Pb_prog2 <- function(){rbeta(1, calc_alpha(mu=prog2, se=0.00034), calc_beta(x=prog2, se=0.00034))} #replaces pb_prog2
Pb_prog3 <- function(){rbeta(1, calc_alpha(mu=prog3, se=0.00162), calc_beta(x=prog3, se=0.00162))} #replaces pb_prog3
Pb_prog1_int <- function(){rbeta(1, calc_alpha(mu=prog1_int, se=0.00036), calc_beta(x=prog1_int, se=0.00036))} #replaces pb_prog1_int
Pb_prog2_int <- function(){rbeta(1, calc_alpha(mu=prog2_int, se=0.00024), calc_beta(x=prog2_int, se=0.00024))} #replaces pb_prog2_int
Pb_prog3_int <- function(){rbeta(1, calc_alpha(mu=prog3_int, se=0.00036), calc_beta(x=prog3_int, se=0.00036))} #replaces pb_prog3_int
Pb_p_mod2 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_mod2, se=0.01090), calc_beta(x=p_mod2, se=0.01090))
} #replaces pb_p_mod2
Pb_p_sev2 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_sev2, se=0.00246), calc_beta(x=p_sev2, se=0.00246))
} #replaces pb_p_sev2
Pb_p_mod3 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_mod3, se=0.01412), calc_beta(x=p_mod3, se=0.01412))
} #replaces pb_p_mod3
Pb_p_sev3 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_sev3, se=0.002250), calc_beta(x=p_sev3, se=0.002250))
} #replaces pb_p_sev3
Pb_p_mod4 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_mod4, se=0.009560), calc_beta(x=p_mod4, se=0.009560))
} #replaces pb_p_mod4
Pb_p_sev4 <- function(){
if (ae_switch == 0) return(0) else 
rbeta(1, calc_alpha(mu=p_sev4, se=0.001890), calc_beta(x=p_sev4, se=0.001890))
} #replaces pb_p_sev4
Pb_u_1 <- function(){rbeta(1, calc_alpha(mu=u_1, se=0.0300), calc_beta(x=u_1, se=0.0300))} #replaces pb_u_1
Pb_u_2 <- function(){rbeta(1, calc_alpha(mu=u_2, se=0.0100), calc_beta(x=u_2, se=0.0100))} #replaces pb_u_2
Pb_u_3 <- function(){rbeta(1, calc_alpha(mu=u_3, se=0.0100), calc_beta(x=u_3, se=0.0100))} #replaces pb_u_3
Pb_u_4 <- function(){rbeta(1, calc_alpha(mu=u_4, se=0.0300), calc_beta(x=u_4, se=0.0300))} #replaces pb_u_4

Pb_c_1 <- function(){
 x <- gamma_params(mean=c_1, sd=20.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_1

Pb_c_2 <- function(){
 x <- gamma_params(mean=c_2, sd=25.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_2

Pb_c_3 <- function(){
 x <- gamma_params(mean=c_3, sd=30.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_3

Pb_c_4 <- function(){
 x <- gamma_params(mean=c_4, sd=35.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_4

Pb_c_mod2 <- function(){
 x <- gamma_params(mean=c_mod2, sd=10.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_mod2

Pb_c_mod3 <- function(){
 x <- gamma_params(mean=c_mod3, sd=10.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_mod3

Pb_c_sev2 <- function(){
 x <- gamma_params(mean=c_sev2, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_sev2

Pb_c_sev3 <- function(){
 x <- gamma_params(mean=c_sev3, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_sev3

Pb_c_mod4 <- function(){
 x <- gamma_params(mean=c_mod4, sd=10.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_mod4

Pb_c_sev4 <- function(){
 x <- gamma_params(mean=c_sev4, sd=500.00)
 rgamma(1, x$alpha, x$beta_tre)
} #replaces pb_c_sev4

x <- gamma_params(mean=c_int, sd=50.00)
pb_c_int <- rgamma(1000, x$alpha, x$beta_tre)
pb_eff_int <- rbeta(1000, calc_alpha(mu=0.22, se=0.07), calc_beta(x=0.22, se=0.07))
pb_eff_uc <- rbeta(1000, calc_alpha(mu=0.06, se=0.02), calc_beta(x=0.06, se=0.02))

n_sim <- 1000
n_cycles <- 240
cl <- makeCluster(round(parallel::detectCores()^0.8)) #This will use e.g. 13 cores if 24 logical cores are available. Uses proportionally more cores if fewer cores are available. Because of the overhead, using all 24 cores probably takes longer than using 13 cores.

registerDoParallel(cl) #Starts the multicore cluster.

costs_smo_probab_new <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% { #replace %dopar% with %do% if you really have only 1 CPU core.
set.seed(i) #Uses a seed of 1 in the first simulation, 2 in the second etc. -> Results are different in every simulation, but still replicable (i.e., identical every time the function is called).
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- numeric(n_cycles + 1)
  
  # Stage 1: initialize and update over cycles
  stage_1[1] <- dat_smo$stage_1[1]
  stage_1[2:n_cycles] <- 100000 * cumprod(1 - Pb_prog1() - dat_smo$mort1)
  
  # Stages 2+3: use initial conditions as in the efficient QALY code
  stage_2[2] <- stage_1[1] * Pb_prog1()  # first cycle where stage 2 > 0
  stage_2[3] <- stage_2[2] * (1 - Pb_prog2() - dat_smo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2] * Pb_prog1()
  stage_3[3] <- stage_2[2] * Pb_prog2()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j] * (1 - Pb_prog2() - dat_smo$mort2[j] - Pb_p_sev2()*ae_mort2[j]) + stage_1[j] * Pb_prog1()
    stage_3[j+1] <- stage_3[j] * (1 - Pb_prog3() - dat_smo$mort3[j] - Pb_p_sev3()*ae_mort3[j]) + stage_2[j] * Pb_prog2()
  }
    
  # Stage 4
  stage_4[4] <- stage_3[4] * Pb_prog3()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j] * (1 - dat_smo$mort3[j] - Pb_p_sev4()*ae_mort3[j]) + stage_3[j] * Pb_prog3()
  }
  
  # Calculate costs at each cycle
  stage_1 * Pb_c_1() +
    stage_2 * (Pb_c_2() + Pb_p_mod2() * Pb_c_mod2() + Pb_p_sev2() * Pb_c_sev2()) +
    stage_3 * (Pb_c_3() + Pb_p_mod3() * Pb_c_mod3() + Pb_p_sev3() * Pb_c_sev3()) +
    stage_4 * (Pb_c_4() + Pb_p_mod4() * Pb_c_mod4() + Pb_p_sev4() * Pb_c_sev4())
}
dim(costs_smo_probab_new) #241 1000
quantile(colSums(costs_smo_probab_new), probs=c(0.025,0.975)) #8 seconds. 95%-CI ca. 150 billion - 205 billion.
#95%-CI ca. 1483487064-2055904427
hist(colSums(costs_smo_probab_new)) #roughly normally distributed
dev.off()

#Analogously for former smokers (costs_fosmo_probab_new):

costs_fosmo_probab_new <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
  # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- numeric(n_cycles + 1)
  
  # Stage 1: initialize and update over cycles
  stage_1[1] <- dat_fosmo$stage_1[1]
  stage_1[2:n_cycles] <- 100000 * cumprod(1 - Pb_prog1_int() - dat_fosmo$mort1)
  
  # Stage 2: use initial conditions as in the efficient QALY code
  stage_2[2] <- stage_1[1] * Pb_prog1_int()  # first cycle where stage 2 > 0
  stage_2[3] <- stage_2[2] * (1 - Pb_prog2_int() - dat_fosmo$mort2[1] - Pb_p_sev2()*0.13883) + stage_1[2] * Pb_prog1_int()
  
  # Stage 3
  stage_3[3] <- stage_2[2] * Pb_prog2_int()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j]*(1 - Pb_prog2_int() - dat_fosmo$mort2[j] - Pb_p_sev2()*ae_mort2_fosmo[j]) + stage_1[j]*Pb_prog1_int()
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3_int() - dat_fosmo$mort3[j] - Pb_p_sev3()*ae_mort3_fosmo[j]) + stage_2[j]*Pb_prog2_int()
  }
  
  # Stage 4
  stage_4[4] <- stage_3[4] * Pb_prog3_int()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j]*(1 - dat_fosmo$mort3[j] - Pb_p_sev4()*ae_mort3_fosmo[j]) + stage_3[j] * Pb_prog3_int()
  }
  
  # Calculate costs at each cycle
  stage_1 * Pb_c_1() +
    stage_2 * (Pb_c_2() + Pb_p_mod2() * Pb_c_mod2() + Pb_p_sev2() * Pb_c_sev2()) +
    stage_3 * (Pb_c_3() + Pb_p_mod3() * Pb_c_mod3() + Pb_p_sev3() * Pb_c_sev3()) +
    stage_4 * (Pb_c_4() + Pb_p_mod4() * Pb_c_mod4() + Pb_p_sev4() * Pb_c_sev4())
}

costs_fosmo_probab_new_backup <- costs_fosmo_probab_new
costs_fosmo_probab_new <- as.numeric(colSums(costs_fosmo_probab_new))
hist(costs_fosmo_probab_new) #roughly normally distributed
quantile(costs_fosmo_probab_new, probs=c(0.025,0.975)) #95% CI ca. 108 billion - 201 billion.
dev.off()

costs_smo_probab_new_backup <- costs_smo_probab_new
costs_smo_probab_new <- as.numeric(colSums(costs_smo_probab_new))
#costs_smo_probab_new - costs_fosmo_probab_new
plot(costs_smo_probab_new - costs_fosmo_probab_new)
abline(h=0) #The intervention is cost-saving in most simulations, because costs_smo_probab_new > costs_fosmo_probab_new.
dev.off()
table(sign((costs_smo_probab_new - costs_fosmo_probab_new))) #More specifically: cost-saving in 84% of the simulations. Which looks plausible.

#QALYs_probab:

qaly_smo_probab_new <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
 # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_smo$stage_1[1]
  stage_1[2:n_cycles] <- 100000 * cumprod(1 - Pb_prog1() - dat_smo$mort1)

  stage_2[2] <- stage_1[1] * Pb_prog1()  # first cycle where stage 2 > 0
  stage_2[3] <- stage_2[2] * (1 - Pb_prog2() - dat_smo$mort2[1] - Pb_p_sev2()*0.13819) + stage_1[2] * Pb_prog1()
  stage_3[3] <- stage_2[2] * Pb_prog2()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j] * (1 - Pb_prog2() - dat_smo$mort2[j] - Pb_p_sev2()*ae_mort2[j]) + stage_1[j] * Pb_prog1()
    stage_3[j+1] <- stage_3[j] * (1 - Pb_prog3() - dat_smo$mort3[j] - Pb_p_sev3()*ae_mort3[j]) + stage_2[j] * Pb_prog2()
  }

  stage_4[4] <- stage_3[4]*Pb_prog3()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j]*(1 - dat_smo$mort3[j] - Pb_p_sev4()*ae_mort3[j]) + stage_3[j]*Pb_prog3()
  }
  (QALY <- ((stage_1*(Pb_u_1()) +
            stage_2*(Pb_u_2() + Pb_p_mod2()*d_mod2 + Pb_p_sev2()*d_sev2)) +
           stage_3*(Pb_u_3() + Pb_p_mod3()*d_mod3 + Pb_p_sev3()*d_sev3) +
           stage_4*(Pb_u_4() + Pb_p_mod4()*d_mod4 + Pb_p_sev4()*d_sev4)) * cycle_length/12)
}
dim(qaly_smo_probab_new) #241 1000, ok. I.e., there are 241 rows / values (=cycles) per column and 1000 columns (=replications/simulations).

for(i in 1:1000){qaly_smo_probab_new[1,i] <- 0} #Calculated as 100000 * p_mild1 * c_mild1 in the original model. But since p_mild1 = 0 (no more mild AEs in the simplified model) -> First row = 0.

quantile(colSums(qaly_smo_probab_new), probs=c(0.025,0.975)) #Almost the same range as in the original model. Looks OK.
hist(colSums(qaly_smo_probab_new)) #Looks OK (largely normally distributed).
dev.off()

#Even faster alternative according to LLM:

# # helper function solving the recursion y[t] = a[t-1] * y[t-1] + u[t]
# compute_stage <- function(u, a) {
  # # u: input vector (length n)
  # # a: decay factors (length n)
  # # build c[k] = prod(a[1:(k-1)]), with c[1] = 1
  # c <- cumprod(c(1, head(a, -1)))
  # # then y[k] = c[k] * sum_{i=1..k}(u[i] / c[i])
  # c * cumsum(u / c)
# }

# qaly_smo_probab_new <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
  # set.seed(i)
  # n <- n_cycles + 1

  # # 1) One-time pre-computations of probabilities & rates
  # p1   <- Pb_prog1()
  # p2   <- Pb_prog2()
  # p3   <- Pb_prog3()
  # sev2 <- Pb_p_sev2()
  # sev3 <- Pb_p_sev3()
  # sev4 <- Pb_p_sev4()
  # mort1 <- dat_smo$mort1
  # mort2 <- dat_smo$mort2
  # mort3 <- dat_smo$mort3
  # ae2   <- ae_mort2
  # ae3   <- ae_mort3

  # # 2) Stage 1 using cumprod (unchanged)
  # stage_1 <- numeric(n)
  # stage_1[1] <- dat_smo$stage_1[1]
  # stage_1[2:n] <- 100000 * cumprod(1 - p1 - mort1)

  # # 3) Stage 2 vectorized
  # a2 <- 1 - p2 - mort2 - sev2 * ae2
  # u2 <- c(0, head(stage_1, -1) * p1)
  # stage_2 <- compute_stage(u2, a2)

  # # 4) Stage 3 vectorized
  # a3 <- 1 - p3 - mort3 - sev3 * ae3
  # u3 <- c(0, head(stage_2, -1) * p2)
  # stage_3 <- compute_stage(u3, a3)

  # # 5) Stage 4 vectorized
  # a4 <- 1 - mort3 - sev4 * ae3
  # u4 <- c(0, head(stage_3, -1) * p3)
  # stage_4 <- compute_stage(u4, a4)

  # # 6) calculate QALY costs per cycle
  # qalys <- stage_1 * Pb_u_1() +
           # stage_2 * (Pb_u_2() + Pb_p_mod2() * d_mod2 + sev2 * d_sev2) +
           # stage_3 * (Pb_u_3() + Pb_p_mod3() * d_mod3 + sev3 * d_sev3) +
           # stage_4 * (Pb_u_4() + Pb_p_mod4() * d_mod4 + sev4 * d_sev4)

  # qalys * cycle_length / 12
# }
## It is indeed much faster. The results are not exactly identical. But the difference is small (~0.1%). Thus, it could be adapted to analogous function calls for costs and former smokers.


#The same for former smokers (fosmo):

qaly_fosmo_probab <- foreach(i = 1:n_sim, .combine = 'cbind') %dopar% {
set.seed(i)
 # Preallocate stage vectors (length n_cycles+1, where n_cycles is 240)
  stage_1 <- stage_2 <- stage_3 <- stage_4 <- numeric(n_cycles + 1)
  
  stage_1[1] <- dat_fosmo$stage_1[1]
  stage_1[2:n_cycles] <- 100000 * cumprod(1 - Pb_prog1_int() - dat_fosmo$mort1)

  # Stage 2: use initial conditions as in the efficient QALY code
  stage_2[2] <- stage_1[1] * Pb_prog1_int()  # first cycle where stage 2 > 0
  stage_2[3] <- stage_2[2] * (1 - Pb_prog2_int() - dat_fosmo$mort2[1] - Pb_p_sev2()*0.13883) + stage_1[2] * Pb_prog1_int()
  
  # Stage 3
  stage_3[3] <- stage_2[2] * Pb_prog2_int()
  for(j in 3:n_cycles){
    stage_2[j+1] <- stage_2[j]*(1 - Pb_prog2_int() - dat_fosmo$mort2[j] - Pb_p_sev2()*ae_mort2_fosmo[j]) + stage_1[j]*Pb_prog1_int()
    stage_3[j+1] <- stage_3[j]*(1 - Pb_prog3_int() - dat_fosmo$mort3[j] - Pb_p_sev3()*ae_mort3_fosmo[j]) + stage_2[j]*Pb_prog2_int()
  }

  stage_4[4] <- stage_3[4]*Pb_prog3_int()
  for(j in 3:n_cycles){
    stage_4[j+1] <- stage_4[j]*(1 - dat_fosmo$mort3[j] - Pb_p_sev4()*ae_mort3_fosmo[j]) + stage_3[j]*Pb_prog3_int()
  }
  (QALY <- ((stage_1*(Pb_u_1()) +
            stage_2*(Pb_u_2() + Pb_p_mod2()*d_mod2 + Pb_p_sev2()*d_sev2)) +
           stage_3*(Pb_u_3() + Pb_p_mod3()*d_mod3 + Pb_p_sev3()*d_sev3) +
           stage_4*(Pb_u_4() + Pb_p_mod4()*d_mod4 + Pb_p_sev4()*d_sev4)) * cycle_length/12)
}
dim(qaly_fosmo_probab) #241 1000, ok
stopCluster(cl) #Stops the multicore cluster.
qaly_fosmo_probab[1,] <- 0 #Calculated as 100000 * p_mild1 * c_mild1 in the original model. But since p_mild1 = 0 (no more mild AEs in the simplified model) -> First row = 0.
 #2352985-2659840. Ok.

## Comparison of the simulation restults, analogously to Fig. 3.3 in Petra Menn's dissertation: ##

#total_QALY_intervention <- sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_int) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_int #13.23588 (OK. Original model: 13.24)

#total_QALY_standard <- sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_uc) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_uc #12.68995 (OK. Original model: 12.70)

#probabilistic:
qaly_smo_probab_new #To be multiplied with (1/((1+disc_e)^(cycle_length/12))^Cycle) for discounting.
qaly_smo_probab_new[,1] * discount_factor
discount_qaly_smo_probab <- function(i){qaly_smo_probab_new[,i] * discount_factor}
discount_qaly_smo_probab(i=1)
cbind(qaly_smo_probab_new[,1], discount_qaly_smo_probab(1)) #Looks OK.
qaly_smo_probab_discounted <- discount_qaly_smo_probab(i=1:1000)

discount_qaly_fosmo_probab <- function(i){qaly_fosmo_probab[,i] * discount_factor}
qaly_fosmo_probab_discounted <- discount_qaly_fosmo_probab(i=1:1000)

#plot(colSums(qaly_smo_probab_discounted)/100000)
costs_smo_probab_new/10000
costs_fosmo_probab_new/10000

plot(colSums(qaly_smo_probab_discounted)/100000, costs_smo_probab_new/10000, xlab="Effects (QALYs)", ylab="Costs (€)")
dev.off()
plot(colSums(qaly_fosmo_probab_discounted)/100000, costs_fosmo_probab_new/10000, xlab="Effects (QALYs)", ylab="Costs (€)")
dev.off()

delta_qalys <- colSums(qaly_fosmo_probab_discounted)/100000 - colSums(qaly_smo_probab_discounted)/100000
delta_costs <- (costs_fosmo_probab_new/100000) - (costs_smo_probab_new/100000)
prop.table(table(sign(delta_costs)))

plot(delta_qalys, delta_costs)
dev.off()

#ICER:
summary(delta_costs / delta_qalys)
quantile(delta_costs / delta_qalys, probs=seq(0,1,0.05))
CEAC <- function(wtp){sum((delta_costs / delta_qalys) < wtp) / length((delta_costs / delta_qalys))}
sapply(seq(100, 1000, 100), CEAC) #OK. As expected, a higher willingness-to-pay goes along with a higher probability of the intervention being cost-effective, with the latter converging towards 1.
plot(sapply(seq(0, 1000, 1), CEAC), type="l", ylim=c(0.5,1), xlab="Willingness-to-pay (€/QALY)", ylab="Probability of cost-effectiveness", main="Cost-effectiveness acceptability curve (CEAC)")
dev.off()

#direct costs, intervention (Int):
(total_dir_costs_intervention <- sum(dat_smo$costs_discounted)/100000 * (1 - eff_int) + sum(dat_fosmo$costs_discounted)/100000 * eff_int + c_int)

#probabilistic:
##First, calculate costs_smo_probab_discounted from costs_smo_probab.

(a <- costs_smo_probab_new_backup)
discount <- function(x, disc_e=0.03){x * discount_factor}
costs_smo_probab_discounted <- discount(a) #Looks OK. Costs increase until ~cycle 32 and decrease thereafter. See plot(apply((costs_smo_probab_discounted), 1, mean))

(b <- costs_fosmo_probab_new_backup)
costs_fosmo_probab_discounted <- discount(b)

#The following is calculated as 100000 * p_mild1 * c_mild1 in the original model. But since p_mild1 = 0 (no more mild AEs in the simplified model) -> First row = 0.
costs_fosmo_probab_discounted[1,] <- costs_smo_probab_discounted[1,] <- 0

## Now make the calculation "direct costs, Intervention (Int):" probabilistic:
(total_dir_costs_intervention_probab <- colSums(costs_smo_probab_discounted)/100000 * (1-pb_eff_int) + colSums(costs_fosmo_probab_discounted)/100000 * pb_eff_int + pb_c_int)

plot(total_dir_costs_intervention_probab)
abline(h=total_dir_costs_intervention)
dev.off()

###########
c_uc <- 0
pb_c_uc <- rep(0, 1000)

#direct costs, standard care (ST):
(total_dir_costs_standard <- sum(dat_smo$costs_discounted)/100000 * (1-eff_uc) + sum(dat_fosmo$costs_discounted)/100000 * eff_uc + c_uc)

(total_dir_costs_standard_probab <- colSums(costs_smo_probab_discounted)/100000 * (1-pb_eff_uc) + colSums(costs_fosmo_probab_discounted)/100000 * pb_eff_uc + pb_c_uc)

plot(total_dir_costs_standard_probab)
abline(h=total_dir_costs_standard) 

#Overall ICER (delta costs / delta effects) for direct costs:
x <- (total_dir_costs_intervention_probab - total_dir_costs_standard_probab) / (qaly_smo_probab_discounted - qaly_fosmo_probab_discounted)
x[1,] <- NA
colSums(x, na.rm=T)
plot(colSums(x, na.rm=T))
summary(colSums(x, na.rm=T))

#Overall ICER (delta costs / delta effects) für indirect costs:
##First calculate indir. costs probabilistically:
###There first express everything in a calculation:
(total_indir_costs_intervention - total_indir_costs_standard) / (total_QALY_intervention - total_QALY_standard) #see ~ line 550.

((sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1 - eff_int) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int) #total_indir_costs_intervention

((sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc) #total_indir_costs_standard

(sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.22) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.22) #total_QALY_intervention

(sum(dat_smo$QALY_discounted) / 100000 * (1 - 0.06) + sum(dat_fosmo$QALY_discounted) / 100000 * 0.06) #total_QALY_standard

(((sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1 - eff_int) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_int + c_int) - ((sum(dat_smo$costs_discounted) / 100000 + sum(dat_smo$ind_costs_discounted) / 100000) * (1-eff_uc) + (sum(dat_fosmo$costs_discounted)/100000 + sum(dat_fosmo$ind_costs_discounted)/100000) * eff_uc + c_uc)) / ((sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_int) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_int) - (sum(dat_smo$QALY_discounted) / 100000 * (1 - eff_uc) + sum(dat_fosmo$QALY_discounted) / 100000 * eff_uc)) #ok.

#probabilistic:
sample(colSums(costs_smo_probab_discounted), 1, replace=TRUE) / 100000 #instead of sum(dat_smo$costs_discounted) / 100000 #ok.
