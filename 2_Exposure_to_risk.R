#######################################################
# Empirical exposure to risk of each subnational region
#######################################################

## Population ratios: Area/Australia_Total (bottom level)

# Row 1
pop_ratio_A_F_to_T = c("pop_ratio_A1_F_to_T", "pop_ratio_A2_F_to_T", "pop_ratio_A3_F_to_T", 
                       "pop_ratio_A4_F_to_T", "pop_ratio_A5_F_to_T", "pop_ratio_A6_F_to_T",
                       "pop_ratio_A7_F_to_T", "pop_ratio_A8_F_to_T", "pop_ratio_A9_F_to_T",
                       "pop_ratio_A10_F_to_T", "pop_ratio_A11_F_to_T", "pop_ratio_A12_F_to_T",
                       "pop_ratio_A13_F_to_T", "pop_ratio_A14_F_to_T", "pop_ratio_A15_F_to_T",
                       "pop_ratio_A16_F_to_T", "pop_ratio_A17_F_to_T", "pop_ratio_A18_F_to_T",
                       "pop_ratio_A19_F_to_T", "pop_ratio_A20_F_to_T", "pop_ratio_A21_F_to_T",
                       "pop_ratio_A22_F_to_T", "pop_ratio_A23_F_to_T", "pop_ratio_A24_F_to_T",
                       "pop_ratio_A25_F_to_T", "pop_ratio_A26_F_to_T", "pop_ratio_A27_F_to_T",
                       "pop_ratio_A28_F_to_T", "pop_ratio_A29_F_to_T", "pop_ratio_A30_F_to_T",
                       "pop_ratio_A31_F_to_T", "pop_ratio_A32_F_to_T", "pop_ratio_A33_F_to_T",
                       "pop_ratio_A34_F_to_T", "pop_ratio_A35_F_to_T", "pop_ratio_A36_F_to_T",
                       "pop_ratio_A37_F_to_T", "pop_ratio_A38_F_to_T", "pop_ratio_A39_F_to_T",
                       "pop_ratio_A40_F_to_T", "pop_ratio_A41_F_to_T", "pop_ratio_A42_F_to_T",
                       "pop_ratio_A43_F_to_T", "pop_ratio_A44_F_to_T", "pop_ratio_A45_F_to_T",
                       "pop_ratio_A46_F_to_T", "pop_ratio_A47_F_to_T")

pop_ratio_A_M_to_T = c("pop_ratio_A1_M_to_T",  "pop_ratio_A2_M_to_T",  "pop_ratio_A3_M_to_T", 
                       "pop_ratio_A4_M_to_T",  "pop_ratio_A5_M_to_T",  "pop_ratio_A6_M_to_T",
                       "pop_ratio_A7_M_to_T",  "pop_ratio_A8_M_to_T",  "pop_ratio_A9_M_to_T",
                       "pop_ratio_A10_M_to_T", "pop_ratio_A11_M_to_T", "pop_ratio_A12_M_to_T",
                       "pop_ratio_A13_M_to_T", "pop_ratio_A14_M_to_T", "pop_ratio_A15_M_to_T",
                       "pop_ratio_A16_M_to_T", "pop_ratio_A17_M_to_T", "pop_ratio_A18_M_to_T",
                       "pop_ratio_A19_M_to_T", "pop_ratio_A20_M_to_T", "pop_ratio_A21_M_to_T",
                       "pop_ratio_A22_M_to_T", "pop_ratio_A23_M_to_T", "pop_ratio_A24_M_to_T",
                       "pop_ratio_A25_M_to_T", "pop_ratio_A26_M_to_T", "pop_ratio_A27_M_to_T",
                       "pop_ratio_A28_M_to_T", "pop_ratio_A29_M_to_T", "pop_ratio_A30_M_to_T",
                       "pop_ratio_A31_M_to_T", "pop_ratio_A32_M_to_T", "pop_ratio_A33_M_to_T",
                       "pop_ratio_A34_M_to_T", "pop_ratio_A35_M_to_T", "pop_ratio_A36_M_to_T",
                       "pop_ratio_A37_M_to_T", "pop_ratio_A38_M_to_T", "pop_ratio_A39_M_to_T",
                       "pop_ratio_A40_M_to_T", "pop_ratio_A41_M_to_T", "pop_ratio_A42_M_to_T",
                       "pop_ratio_A43_M_to_T", "pop_ratio_A44_M_to_T", "pop_ratio_A45_M_to_T",
                       "pop_ratio_A46_M_to_T", "pop_ratio_A47_M_to_T")


for(iw in 2:48)
{
  assign(pop_ratio_A_F_to_T[iw-1], (get(Area[iw])$pop$female/Australia$pop$total)[,20:24])
  assign(pop_ratio_A_M_to_T[iw-1], (get(Area[iw])$pop$male/Australia$pop$total)[,20:24])
}

# Row 2

pop_ratio_A_F_to_F = c("pop_ratio_A1_F_to_F", "pop_ratio_A2_F_to_F", "pop_ratio_A3_F_to_F", 
                       "pop_ratio_A4_F_to_F", "pop_ratio_A5_F_to_F", "pop_ratio_A6_F_to_F",
                       "pop_ratio_A7_F_to_F", "pop_ratio_A8_F_to_F", "pop_ratio_A9_F_to_F",
                       "pop_ratio_A10_F_to_F", "pop_ratio_A11_F_to_F", "pop_ratio_A12_F_to_F",
                       "pop_ratio_A13_F_to_F", "pop_ratio_A14_F_to_F", "pop_ratio_A15_F_to_F",
                       "pop_ratio_A16_F_to_F", "pop_ratio_A17_F_to_F", "pop_ratio_A18_F_to_F",
                       "pop_ratio_A19_F_to_F", "pop_ratio_A20_F_to_F", "pop_ratio_A21_F_to_F",
                       "pop_ratio_A22_F_to_F", "pop_ratio_A23_F_to_F", "pop_ratio_A24_F_to_F",
                       "pop_ratio_A25_F_to_F", "pop_ratio_A26_F_to_F", "pop_ratio_A27_F_to_F",
                       "pop_ratio_A28_F_to_F", "pop_ratio_A29_F_to_F", "pop_ratio_A30_F_to_F",
                       "pop_ratio_A31_F_to_F", "pop_ratio_A32_F_to_F", "pop_ratio_A33_F_to_F",
                       "pop_ratio_A34_F_to_F", "pop_ratio_A35_F_to_F", "pop_ratio_A36_F_to_F",
                       "pop_ratio_A37_F_to_F", "pop_ratio_A38_F_to_F", "pop_ratio_A39_F_to_F",
                       "pop_ratio_A40_F_to_F", "pop_ratio_A41_F_to_F", "pop_ratio_A42_F_to_F",
                       "pop_ratio_A43_F_to_F", "pop_ratio_A44_F_to_F", "pop_ratio_A45_F_to_F",
                       "pop_ratio_A46_F_to_F", "pop_ratio_A47_F_to_F")

# Row 3

pop_ratio_A_M_to_M = c("pop_ratio_A1_M_to_M",  "pop_ratio_A2_M_to_M",  "pop_ratio_A3_M_to_M", 
                       "pop_ratio_A4_M_to_M",  "pop_ratio_A5_M_to_M",  "pop_ratio_A6_M_to_M",
                       "pop_ratio_A7_M_to_M",  "pop_ratio_A8_M_to_M",  "pop_ratio_A9_M_to_M",
                       "pop_ratio_A10_M_to_M", "pop_ratio_A11_M_to_M", "pop_ratio_A12_M_to_M",
                       "pop_ratio_A13_M_to_M", "pop_ratio_A14_M_to_M", "pop_ratio_A15_M_to_M",
                       "pop_ratio_A16_M_to_M", "pop_ratio_A17_M_to_M", "pop_ratio_A18_M_to_M",
                       "pop_ratio_A19_M_to_M", "pop_ratio_A20_M_to_M", "pop_ratio_A21_M_to_M",
                       "pop_ratio_A22_M_to_M", "pop_ratio_A23_M_to_M", "pop_ratio_A24_M_to_M",
                       "pop_ratio_A25_M_to_M", "pop_ratio_A26_M_to_M", "pop_ratio_A27_M_to_M",
                       "pop_ratio_A28_M_to_M", "pop_ratio_A29_M_to_M", "pop_ratio_A30_M_to_M",
                       "pop_ratio_A31_M_to_M", "pop_ratio_A32_M_to_M", "pop_ratio_A33_M_to_M",
                       "pop_ratio_A34_M_to_M", "pop_ratio_A35_M_to_M", "pop_ratio_A36_M_to_M",
                       "pop_ratio_A37_M_to_M", "pop_ratio_A38_M_to_M", "pop_ratio_A39_M_to_M",
                       "pop_ratio_A40_M_to_M", "pop_ratio_A41_M_to_M", "pop_ratio_A42_M_to_M",
                       "pop_ratio_A43_M_to_M", "pop_ratio_A44_M_to_M", "pop_ratio_A45_M_to_M",
                       "pop_ratio_A46_M_to_M", "pop_ratio_A47_M_to_M")

for(iw in 2:48)
{
  assign(pop_ratio_A_F_to_F[iw-1], (get(Area[iw])$pop$female/Australia$pop$female)[,20:24])
  assign(pop_ratio_A_M_to_M[iw-1], (get(Area[iw])$pop$male/Australia$pop$male)[,20:24])
}

## Population ratios: Region/Australia_Total (middle row)

# Row 4

pop_ratio_A_F_to_R1_T = "pop_ratio_A1_F_to_R1_T"
pop_ratio_A_F_to_R2_T = c("pop_ratio_A2_F_to_R2_T", "pop_ratio_A3_F_to_R2_T", "pop_ratio_A4_F_to_R2_T")
pop_ratio_A_F_to_R3_T = "pop_ratio_A5_F_to_R3_T"
pop_ratio_A_F_to_R4_T = c("pop_ratio_A6_F_to_R4_T", "pop_ratio_A7_F_to_R4_T", "pop_ratio_A8_F_to_R4_T",
                          "pop_ratio_A9_F_to_R4_T", "pop_ratio_A10_F_to_R4_T")
pop_ratio_A_F_to_R5_T = "pop_ratio_A11_F_to_R5_T"
pop_ratio_A_F_to_R6_T = "pop_ratio_A12_F_to_R6_T"
pop_ratio_A_F_to_R7_T = "pop_ratio_A13_F_to_R7_T"
pop_ratio_A_F_to_R8_T = "pop_ratio_A14_F_to_R8_T"
pop_ratio_A_F_to_R9_T = "pop_ratio_A15_F_to_R9_T"
pop_ratio_A_F_to_R10_T = c("pop_ratio_A16_F_to_R10_T", "pop_ratio_A17_F_to_R10_T", "pop_ratio_A18_F_to_R10_T", 
                           "pop_ratio_A19_F_to_R10_T", "pop_ratio_A20_F_to_R10_T", "pop_ratio_A21_F_to_R10_T", 
                           "pop_ratio_A22_F_to_R10_T", "pop_ratio_A23_F_to_R10_T", "pop_ratio_A24_F_to_R10_T", 
                           "pop_ratio_A25_F_to_R10_T", "pop_ratio_A26_F_to_R10_T", "pop_ratio_A27_F_to_R10_T", 
                           "pop_ratio_A28_F_to_R10_T", "pop_ratio_A29_F_to_R10_T", "pop_ratio_A30_F_to_R10_T", 
                           "pop_ratio_A31_F_to_R10_T", "pop_ratio_A32_F_to_R10_T", "pop_ratio_A33_F_to_R10_T", 
                           "pop_ratio_A34_F_to_R10_T")
pop_ratio_A_F_to_R11_T = c("pop_ratio_A35_F_to_R11_T", "pop_ratio_A36_F_to_R11_T", "pop_ratio_A37_F_to_R11_T", 
                           "pop_ratio_A38_F_to_R11_T", "pop_ratio_A39_F_to_R11_T", "pop_ratio_A40_F_to_R11_T", 
                           "pop_ratio_A41_F_to_R11_T", "pop_ratio_A42_F_to_R11_T", "pop_ratio_A43_F_to_R11_T", 
                           "pop_ratio_A44_F_to_R11_T", "pop_ratio_A45_F_to_R11_T", "pop_ratio_A46_F_to_R11_T", 
                           "pop_ratio_A47_F_to_R11_T")


assign(pop_ratio_A_F_to_R1_T, (get(Area[2])$pop$female/R01$pop$total)[,20:24])
for(iw in 2:4)
{
  assign(pop_ratio_A_F_to_R2_T[iw-1], (get(Area[(iw+1)])$pop$female/R02$pop$total)[,20:24])
}         
assign(pop_ratio_A_F_to_R3_T, (get(Area[6])$pop$female/R03$pop$total)[,20:24])
for(iw in 6:10)
{
  assign(pop_ratio_A_F_to_R4_T[iw-5], (get(Area[(iw+1)])$pop$female/R04$pop$total)[,20:24])
}
assign(pop_ratio_A_F_to_R5_T, (get(Area[12])$pop$female/R05$pop$total)[,20:24])
assign(pop_ratio_A_F_to_R6_T, (get(Area[13])$pop$female/R06$pop$total)[,20:24])
assign(pop_ratio_A_F_to_R7_T, (get(Area[14])$pop$female/R07$pop$total)[,20:24])
assign(pop_ratio_A_F_to_R8_T, (get(Area[15])$pop$female/R08$pop$total)[,20:24])
assign(pop_ratio_A_F_to_R9_T, (get(Area[16])$pop$female/R09$pop$total)[,20:24])
for(iw in 16:34)
{
  assign(pop_ratio_A_F_to_R10_T[iw-15], (get(Area[(iw+1)])$pop$female/R10$pop$total)[,20:24])
}
for(iw in 35:47)
{
  assign(pop_ratio_A_F_to_R11_T[iw-34], (get(Area[(iw+1)])$pop$female/R11$pop$total)[,20:24])
}


pop_ratio_A_M_to_R1_T = "pop_ratio_A1_M_to_R1_T"
pop_ratio_A_M_to_R2_T = c("pop_ratio_A2_M_to_R2_T", "pop_ratio_A3_M_to_R2_T", "pop_ratio_A4_M_to_R2_T")
pop_ratio_A_M_to_R3_T = "pop_ratio_A5_M_to_R3_T"
pop_ratio_A_M_to_R4_T = c("pop_ratio_A6_M_to_R4_T", "pop_ratio_A7_M_to_R4_T", "pop_ratio_A8_M_to_R4_T",
                          "pop_ratio_A9_M_to_R4_T", "pop_ratio_A10_M_to_R4_T")
pop_ratio_A_M_to_R5_T = "pop_ratio_A11_M_to_R5_T"
pop_ratio_A_M_to_R6_T = "pop_ratio_A12_M_to_R6_T"
pop_ratio_A_M_to_R7_T = "pop_ratio_A13_M_to_R7_T"
pop_ratio_A_M_to_R8_T = "pop_ratio_A14_M_to_R8_T"
pop_ratio_A_M_to_R9_T = "pop_ratio_A15_M_to_R9_T"
pop_ratio_A_M_to_R10_T = c("pop_ratio_A16_M_to_R10_T", "pop_ratio_A17_M_to_R10_T", "pop_ratio_A18_M_to_R10_T", 
                           "pop_ratio_A19_M_to_R10_T", "pop_ratio_A20_M_to_R10_T", "pop_ratio_A21_M_to_R10_T", 
                           "pop_ratio_A22_M_to_R10_T", "pop_ratio_A23_M_to_R10_T", "pop_ratio_A24_M_to_R10_T", 
                           "pop_ratio_A25_M_to_R10_T", "pop_ratio_A26_M_to_R10_T", "pop_ratio_A27_M_to_R10_T", 
                           "pop_ratio_A28_M_to_R10_T", "pop_ratio_A29_M_to_R10_T", "pop_ratio_A30_M_to_R10_T", 
                           "pop_ratio_A31_M_to_R10_T", "pop_ratio_A32_M_to_R10_T", "pop_ratio_A33_M_to_R10_T", 
                           "pop_ratio_A34_M_to_R10_T")
pop_ratio_A_M_to_R11_T = c("pop_ratio_A35_M_to_R11_T", "pop_ratio_A36_M_to_R11_T", "pop_ratio_A37_M_to_R11_T", 
                           "pop_ratio_A38_M_to_R11_T", "pop_ratio_A39_M_to_R11_T", "pop_ratio_A40_M_to_R11_T", 
                           "pop_ratio_A41_M_to_R11_T", "pop_ratio_A42_M_to_R11_T", "pop_ratio_A43_M_to_R11_T", 
                           "pop_ratio_A44_M_to_R11_T", "pop_ratio_A45_M_to_R11_T", "pop_ratio_A46_M_to_R11_T", 
                           "pop_ratio_A47_M_to_R11_T")


assign(pop_ratio_A_M_to_R1_T, (get(Area[2])$pop$male/R01$pop$total)[,20:24])
for(iw in 2:4)
{
  assign(pop_ratio_A_M_to_R2_T[iw-1], (get(Area[(iw+1)])$pop$male/R02$pop$total)[,20:24])
}         
assign(pop_ratio_A_M_to_R3_T, (get(Area[6])$pop$male/R03$pop$total)[,20:24])
for(iw in 6:10)
{
  assign(pop_ratio_A_M_to_R4_T[iw-5], (get(Area[(iw+1)])$pop$male/R04$pop$total)[,20:24])
}
assign(pop_ratio_A_M_to_R5_T, (get(Area[12])$pop$male/R05$pop$total)[,20:24])
assign(pop_ratio_A_M_to_R6_T, (get(Area[13])$pop$male/R06$pop$total)[,20:24])
assign(pop_ratio_A_M_to_R7_T, (get(Area[14])$pop$male/R07$pop$total)[,20:24])
assign(pop_ratio_A_M_to_R8_T, (get(Area[15])$pop$male/R08$pop$total)[,20:24])
assign(pop_ratio_A_M_to_R9_T, (get(Area[16])$pop$male/R09$pop$total)[,20:24])
for(iw in 16:34)
{
  assign(pop_ratio_A_M_to_R10_T[iw-15], (get(Area[(iw+1)])$pop$male/R10$pop$total)[,20:24])
}
for(iw in 35:47)
{
  assign(pop_ratio_A_M_to_R11_T[iw-34], (get(Area[(iw+1)])$pop$male/R11$pop$total)[,20:24])
}

# Row 5

pop_ratio_A_F_to_R1_F = "pop_ratio_A1_F_to_R1_F"
pop_ratio_A_F_to_R2_F = c("pop_ratio_A2_F_to_R2_F", "pop_ratio_A3_F_to_R2_F", "pop_ratio_A4_F_to_R2_F")
pop_ratio_A_F_to_R3_F = "pop_ratio_A5_F_to_R3_F"
pop_ratio_A_F_to_R4_F = c("pop_ratio_A6_F_to_R4_F", "pop_ratio_A7_F_to_R4_F", "pop_ratio_A8_F_to_R4_F",
                          "pop_ratio_A9_F_to_R4_F", "pop_ratio_A10_F_to_R4_F")
pop_ratio_A_F_to_R5_F = "pop_ratio_A11_F_to_R5_F"
pop_ratio_A_F_to_R6_F = "pop_ratio_A12_F_to_R6_F"
pop_ratio_A_F_to_R7_F = "pop_ratio_A13_F_to_R7_F"
pop_ratio_A_F_to_R8_F = "pop_ratio_A14_F_to_R8_F"
pop_ratio_A_F_to_R9_F = "pop_ratio_A15_F_to_R9_F"
pop_ratio_A_F_to_R10_F = c("pop_ratio_A16_F_to_R10_F", "pop_ratio_A17_F_to_R10_F", "pop_ratio_A18_F_to_R10_F", 
                           "pop_ratio_A19_F_to_R10_F", "pop_ratio_A20_F_to_R10_F", "pop_ratio_A21_F_to_R10_F", 
                           "pop_ratio_A22_F_to_R10_F", "pop_ratio_A23_F_to_R10_F", "pop_ratio_A24_F_to_R10_F", 
                           "pop_ratio_A25_F_to_R10_F", "pop_ratio_A26_F_to_R10_F", "pop_ratio_A27_F_to_R10_F", 
                           "pop_ratio_A28_F_to_R10_F", "pop_ratio_A29_F_to_R10_F", "pop_ratio_A30_F_to_R10_F", 
                           "pop_ratio_A31_F_to_R10_F", "pop_ratio_A32_F_to_R10_F", "pop_ratio_A33_F_to_R10_F", 
                           "pop_ratio_A34_F_to_R10_F")
pop_ratio_A_F_to_R11_F = c("pop_ratio_A35_F_to_R11_F", "pop_ratio_A36_F_to_R11_F", "pop_ratio_A37_F_to_R11_F", 
                           "pop_ratio_A38_F_to_R11_F", "pop_ratio_A39_F_to_R11_F", "pop_ratio_A40_F_to_R11_F", 
                           "pop_ratio_A41_F_to_R11_F", "pop_ratio_A42_F_to_R11_F", "pop_ratio_A43_F_to_R11_F", 
                           "pop_ratio_A44_F_to_R11_F", "pop_ratio_A45_F_to_R11_F", "pop_ratio_A46_F_to_R11_F", 
                           "pop_ratio_A47_F_to_R11_F")


assign(pop_ratio_A_F_to_R1_F, (get(Area[2])$pop$female/R01$pop$female)[,20:24])
for(iw in 2:4)
{
  assign(pop_ratio_A_F_to_R2_F[iw-1], (get(Area[(iw+1)])$pop$female/R02$pop$female)[,20:24])
}         
assign(pop_ratio_A_F_to_R3_F, (get(Area[6])$pop$female/R03$pop$female)[,20:24])
for(iw in 6:10)
{
  assign(pop_ratio_A_F_to_R4_F[iw-5], (get(Area[(iw+1)])$pop$female/R04$pop$female)[,20:24])
}
assign(pop_ratio_A_F_to_R5_F, (get(Area[12])$pop$female/R05$pop$female)[,20:24])
assign(pop_ratio_A_F_to_R6_F, (get(Area[13])$pop$female/R06$pop$female)[,20:24])
assign(pop_ratio_A_F_to_R7_F, (get(Area[14])$pop$female/R07$pop$female)[,20:24])
assign(pop_ratio_A_F_to_R8_F, (get(Area[15])$pop$female/R08$pop$female)[,20:24])
assign(pop_ratio_A_F_to_R9_F, (get(Area[16])$pop$female/R09$pop$female)[,20:24])
for(iw in 16:34)
{
  assign(pop_ratio_A_F_to_R10_F[iw-15], (get(Area[(iw+1)])$pop$female/R10$pop$female)[,20:24])
}
for(iw in 35:47)
{
  assign(pop_ratio_A_F_to_R11_F[iw-34], (get(Area[(iw+1)])$pop$female/R11$pop$female)[,20:24])
}

# Row 6

pop_ratio_A_M_to_R1_M = "pop_ratio_A1_M_to_R1_M"
pop_ratio_A_M_to_R2_M = c("pop_ratio_A2_M_to_R2_M", "pop_ratio_A3_M_to_R2_M", "pop_ratio_A4_M_to_R2_M")
pop_ratio_A_M_to_R3_M = "pop_ratio_A5_M_to_R3_M"
pop_ratio_A_M_to_R4_M = c("pop_ratio_A6_M_to_R4_M", "pop_ratio_A7_M_to_R4_M", "pop_ratio_A8_M_to_R4_M",
                          "pop_ratio_A9_M_to_R4_M", "pop_ratio_A10_M_to_R4_M")
pop_ratio_A_M_to_R5_M = "pop_ratio_A11_M_to_R5_M"
pop_ratio_A_M_to_R6_M = "pop_ratio_A12_M_to_R6_M"
pop_ratio_A_M_to_R7_M = "pop_ratio_A13_M_to_R7_M"
pop_ratio_A_M_to_R8_M = "pop_ratio_A14_M_to_R8_M"
pop_ratio_A_M_to_R9_M = "pop_ratio_A15_M_to_R9_M"
pop_ratio_A_M_to_R10_M = c("pop_ratio_A16_M_to_R10_M", "pop_ratio_A17_M_to_R10_M", "pop_ratio_A18_M_to_R10_M", 
                           "pop_ratio_A19_M_to_R10_M", "pop_ratio_A20_M_to_R10_M", "pop_ratio_A21_M_to_R10_M", 
                           "pop_ratio_A22_M_to_R10_M", "pop_ratio_A23_M_to_R10_M", "pop_ratio_A24_M_to_R10_M", 
                           "pop_ratio_A25_M_to_R10_M", "pop_ratio_A26_M_to_R10_M", "pop_ratio_A27_M_to_R10_M", 
                           "pop_ratio_A28_M_to_R10_M", "pop_ratio_A29_M_to_R10_M", "pop_ratio_A30_M_to_R10_M", 
                           "pop_ratio_A31_M_to_R10_M", "pop_ratio_A32_M_to_R10_M", "pop_ratio_A33_M_to_R10_M", 
                           "pop_ratio_A34_M_to_R10_M")
pop_ratio_A_M_to_R11_M = c("pop_ratio_A35_M_to_R11_M", "pop_ratio_A36_M_to_R11_M", "pop_ratio_A37_M_to_R11_M", 
                           "pop_ratio_A38_M_to_R11_M", "pop_ratio_A39_M_to_R11_M", "pop_ratio_A40_M_to_R11_M", 
                           "pop_ratio_A41_M_to_R11_M", "pop_ratio_A42_M_to_R11_M", "pop_ratio_A43_M_to_R11_M", 
                           "pop_ratio_A44_M_to_R11_M", "pop_ratio_A45_M_to_R11_M", "pop_ratio_A46_M_to_R11_M", 
                           "pop_ratio_A47_M_to_R11_M")


assign(pop_ratio_A_M_to_R1_M, (get(Area[2])$pop$male/R01$pop$male)[,20:24])
for(iw in 2:4)
{
  assign(pop_ratio_A_M_to_R2_M[iw-1], (get(Area[(iw+1)])$pop$male/R02$pop$male)[,20:24])
}         
assign(pop_ratio_A_M_to_R3_M, (get(Area[6])$pop$male/R03$pop$male)[,20:24])
for(iw in 6:10)
{
  assign(pop_ratio_A_M_to_R4_M[iw-5], (get(Area[(iw+1)])$pop$male/R04$pop$male)[,20:24])
}
assign(pop_ratio_A_M_to_R5_M, (get(Area[12])$pop$male/R05$pop$male)[,20:24])
assign(pop_ratio_A_M_to_R6_M, (get(Area[13])$pop$male/R06$pop$male)[,20:24])
assign(pop_ratio_A_M_to_R7_M, (get(Area[14])$pop$male/R07$pop$male)[,20:24])
assign(pop_ratio_A_M_to_R8_M, (get(Area[15])$pop$male/R08$pop$male)[,20:24])
assign(pop_ratio_A_M_to_R9_M, (get(Area[16])$pop$male/R09$pop$male)[,20:24])
for(iw in 16:34)
{
  assign(pop_ratio_A_M_to_R10_M[iw-15], (get(Area[(iw+1)])$pop$male/R10$pop$male)[,20:24])
}
for(iw in 35:47)
{
  assign(pop_ratio_A_M_to_R11_M[iw-34], (get(Area[(iw+1)])$pop$male/R11$pop$male)[,20:24])
}

## Population ratios: Area/Area_total (top level)

# Row 7

pop_ratio_A_F_to_A_T = c("pop_ratio_A1_F_to_A1_T", "pop_ratio_A2_F_to_A2_T", "pop_ratio_A3_F_to_A3_T",
                         "pop_ratio_A4_F_to_A4_T", "pop_ratio_A5_F_to_A5_T", "pop_ratio_A6_F_to_A6_T",
                         "pop_ratio_A7_F_to_A7_T", "pop_ratio_A8_F_to_A8_T", "pop_ratio_A9_F_to_A9_T",
                         "pop_ratio_A10_F_to_A10_T", "pop_ratio_A11_F_to_A11_T", "pop_ratio_A12_F_to_A12_T",
                         "pop_ratio_A13_F_to_A13_T", "pop_ratio_A14_F_to_A14_T",
                         "pop_ratio_A15_F_to_A15_T", "pop_ratio_A16_F_to_A16_T", "pop_ratio_A17_F_to_A17_T",
                         "pop_ratio_A18_F_to_A18_T", "pop_ratio_A19_F_to_A19_T", "pop_ratio_A20_F_to_A20_T",
                         "pop_ratio_A21_F_to_A21_T", "pop_ratio_A22_F_to_A22_T", "pop_ratio_A23_F_to_A23_T",
                         "pop_ratio_A24_F_to_A24_T", "pop_ratio_A25_F_to_A25_T", "pop_ratio_A26_F_to_A26_T",
                         "pop_ratio_A27_F_to_A27_T", "pop_ratio_A28_F_to_A28_T", "pop_ratio_A29_F_to_A29_T",
                         "pop_ratio_A30_F_to_A30_T", "pop_ratio_A31_F_to_A31_T", "pop_ratio_A32_F_to_A32_T",
                         "pop_ratio_A33_F_to_A33_T", "pop_ratio_A34_F_to_A34_T", "pop_ratio_A35_F_to_A35_T",
                         "pop_ratio_A36_F_to_A36_T", "pop_ratio_A37_F_to_A37_T", "pop_ratio_A38_F_to_A38_T",
                         "pop_ratio_A39_F_to_A39_T", "pop_ratio_A40_F_to_A40_T", "pop_ratio_A41_F_to_A41_T",
                         "pop_ratio_A42_F_to_A42_T", "pop_ratio_A43_F_to_A43_T", "pop_ratio_A44_F_to_A44_T",
                         "pop_ratio_A45_F_to_A45_T", "pop_ratio_A46_F_to_A46_T", "pop_ratio_A47_F_to_A47_T")

pop_ratio_A_M_to_A_T = c("pop_ratio_A1_M_to_A1_T", "pop_ratio_A2_M_to_A2_T", "pop_ratio_A3_M_to_A3_T",
                         "pop_ratio_A4_M_to_A4_T", "pop_ratio_A5_M_to_A5_T", "pop_ratio_A6_M_to_A6_T",
                         "pop_ratio_A7_M_to_A7_T", "pop_ratio_A8_M_to_A8_T", "pop_ratio_A9_M_to_A9_T",
                         "pop_ratio_A10_M_to_A10_T", "pop_ratio_A11_M_to_A11_T", "pop_ratio_A12_M_to_A12_T",
                         "pop_ratio_A13_M_to_A13_T", "pop_ratio_A14_M_to_A14_T",
                         "pop_ratio_A15_M_to_A15_T", "pop_ratio_A16_M_to_A16_T", "pop_ratio_A17_M_to_A17_T",
                         "pop_ratio_A18_M_to_A18_T", "pop_ratio_A19_M_to_A19_T", "pop_ratio_A20_M_to_A20_T",
                         "pop_ratio_A21_M_to_A21_T", "pop_ratio_A22_M_to_A22_T", "pop_ratio_A23_M_to_A23_T",
                         "pop_ratio_A24_M_to_A24_T", "pop_ratio_A25_M_to_A25_T", "pop_ratio_A26_M_to_A26_T",
                         "pop_ratio_A27_M_to_A27_T", "pop_ratio_A28_M_to_A28_T", "pop_ratio_A29_M_to_A29_T",
                         "pop_ratio_A30_M_to_A30_T", "pop_ratio_A31_M_to_A31_T", "pop_ratio_A32_M_to_A32_T",
                         "pop_ratio_A33_M_to_A33_T", "pop_ratio_A34_M_to_A34_T", "pop_ratio_A35_M_to_A35_T",
                         "pop_ratio_A36_M_to_A36_T", "pop_ratio_A37_M_to_A37_T", "pop_ratio_A38_M_to_A38_T",
                         "pop_ratio_A39_M_to_A39_T", "pop_ratio_A40_M_to_A40_T", "pop_ratio_A41_M_to_A41_T",
                         "pop_ratio_A42_M_to_A42_T", "pop_ratio_A43_M_to_A43_T", "pop_ratio_A44_M_to_A44_T",
                         "pop_ratio_A45_M_to_A45_T", "pop_ratio_A46_M_to_A46_T", "pop_ratio_A47_M_to_A47_T")

for(iw in 2:48)
{
  assign(pop_ratio_A_F_to_A_T[iw-1], (get(Area[iw])$pop$female/get(Area[iw])$pop$total)[,20:24])
  assign(pop_ratio_A_M_to_A_T[iw-1], (get(Area[iw])$pop$male/get(Area[iw])$pop$total)[,20:24])          
}

















