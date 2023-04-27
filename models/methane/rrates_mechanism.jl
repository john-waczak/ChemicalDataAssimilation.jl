my_ones = ones(nrow(df_params))


df_mech_rate_coeffs = DataFrame()
# copy in the times
df_mech_rate_coeffs.t = df_params.t
df_mech_rate_coeffs.w_ap = df_params.w_ap
df_mech_rate_coeffs.k_1 = @. 5.6e-34*df_params.N2*(df_params.temperature/300)^-2.6*df_params.O2*my_ones
df_mech_rate_coeffs.k_2 = @. 6.0e-34*df_params.O2*(df_params.temperature/300)^-2.6*df_params.O2*my_ones
df_mech_rate_coeffs.k_3 = @. 8.0e-12*exp(-2060/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_4 = @. df_rrate_coeffs.KMT01*my_ones
df_mech_rate_coeffs.k_5 = @. 5.5e-12*exp(188/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_6 = @. df_rrate_coeffs.KMT02*my_ones
df_mech_rate_coeffs.k_7 = @. 3.2e-11*exp(67/df_params.temperature)*df_params.O2*my_ones
df_mech_rate_coeffs.k_8 = @. 2.0e-11*exp(130/df_params.temperature)*df_params.N2*my_ones
df_mech_rate_coeffs.k_9 = @. 1.4e-12*exp(-1310/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_10 = @. 1.4e-13*exp(-2470/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_11 = @. 3.3e-39*exp(530/df_params.temperature)*df_params.O2*my_ones
df_mech_rate_coeffs.k_12 = @. 1.8e-11*exp(110/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_13 = @. 4.50e-14*exp(-1260/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_14 = @. df_rrate_coeffs.KMT03*my_ones
df_mech_rate_coeffs.k_15 = @. 2.14e-10*df_params.H2O*my_ones
df_mech_rate_coeffs.k_16 = @. 1.70e-12*exp(-940/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_17 = @. 7.7e-12*exp(-2100/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_18 = @. df_rrate_coeffs.KMT05*my_ones
df_mech_rate_coeffs.k_19 = @. 2.9e-12*exp(-160/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_20 = @. 2.03e-16*(df_params.temperature/300)^4.57*exp(693/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_21 = @. 4.8e-11*exp(250/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_22 = @. 2.20e-13*df_rrate_coeffs.KMT06*exp(600/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_23 = @. 1.90e-33*df_params.M*df_rrate_coeffs.KMT06*exp(980/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_24 = @. df_rrate_coeffs.KMT07*my_ones
df_mech_rate_coeffs.k_25 = @. df_rrate_coeffs.KMT08*my_ones
df_mech_rate_coeffs.k_26 = @. 2.0e-11*my_ones
df_mech_rate_coeffs.k_27 = @. 3.45e-12*exp(270/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_28 = @. df_rrate_coeffs.KMT09*my_ones
df_mech_rate_coeffs.k_29 = @. 3.2e-13*exp(690/df_params.temperature)*1.0*my_ones
df_mech_rate_coeffs.k_30 = @. 4.0e-12*my_ones
df_mech_rate_coeffs.k_31 = @. 2.5e-12*exp(260/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_32 = @. df_rrate_coeffs.KMT11*my_ones
df_mech_rate_coeffs.k_33 = @. 4.0e-32*exp(-1000/df_params.temperature)*df_params.M*my_ones
df_mech_rate_coeffs.k_34 = @. df_rrate_coeffs.KMT12*my_ones
df_mech_rate_coeffs.k_35 = @. 1.3e-12*exp(-330/df_params.temperature)*df_params.O2*my_ones
df_mech_rate_coeffs.k_36 = @. 6.00e-06*my_ones
df_mech_rate_coeffs.k_37 = @. 4.00e-04*my_ones
df_mech_rate_coeffs.k_38 = @. 1.20e-15*df_params.H2O*my_ones
df_mech_rate_coeffs.k_39 = @. df_photolysis.J_1*my_ones
df_mech_rate_coeffs.k_40 = @. df_photolysis.J_2*my_ones
df_mech_rate_coeffs.k_41 = @. df_photolysis.J_3*my_ones
df_mech_rate_coeffs.k_42 = @. df_photolysis.J_4*my_ones
df_mech_rate_coeffs.k_43 = @. df_photolysis.J_5*my_ones
df_mech_rate_coeffs.k_44 = @. df_photolysis.J_6*my_ones
df_mech_rate_coeffs.k_45 = @. df_photolysis.J_7*my_ones
df_mech_rate_coeffs.k_46 = @. df_photolysis.J_8*my_ones
df_mech_rate_coeffs.k_47 = @. df_rrate_coeffs.KMT04*my_ones
df_mech_rate_coeffs.k_48 = @. df_rrate_coeffs.KMT10*my_ones
df_mech_rate_coeffs.k_49 = @. 6.6e-12*exp(-1240/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_50 = @. 1.85e-12*exp(-1690/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_51 = @. 2.3e-12*exp(360/df_params.temperature)*0.001*my_ones
df_mech_rate_coeffs.k_52 = @. 2.3e-12*exp(360/df_params.temperature)*0.999*my_ones
df_mech_rate_coeffs.k_53 = @. df_rrate_coeffs.KMT13*my_ones
df_mech_rate_coeffs.k_54 = @. 1.2e-12*my_ones
df_mech_rate_coeffs.k_55 = @. 2*df_rrate_coeffs.KCH3O2*RO2ᵢ*7.18*exp(-885/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_56 = @. 2*df_rrate_coeffs.KCH3O2*RO2ᵢ*0.5*(1-7.18*exp(-885/df_params.temperature))*my_ones
df_mech_rate_coeffs.k_57 = @. 2*df_rrate_coeffs.KCH3O2*RO2ᵢ*0.5*(1-7.18*exp(-885/df_params.temperature))*my_ones
df_mech_rate_coeffs.k_58 = @. df_photolysis.J_41*my_ones
df_mech_rate_coeffs.k_59 = @. 5.3e-12*exp(190/df_params.temperature)*0.6*my_ones
df_mech_rate_coeffs.k_60 = @. 5.3e-12*exp(190/df_params.temperature)*0.4*my_ones
df_mech_rate_coeffs.k_61 = @. df_photolysis.J_11*my_ones
df_mech_rate_coeffs.k_62 = @. df_photolysis.J_12*my_ones
df_mech_rate_coeffs.k_63 = @. 5.5e-16*my_ones
df_mech_rate_coeffs.k_64 = @. 5.4e-12*exp(135/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_65 = @. df_photolysis.J_51*my_ones
df_mech_rate_coeffs.k_66 = @. 4.0e-13*exp(-845/df_params.temperature)*my_ones
df_mech_rate_coeffs.k_67 = @. 7.2e-14*exp(-1080/df_params.temperature)*df_params.O2*my_ones
df_mech_rate_coeffs.k_68 = @. df_rrate_coeffs.KMT14*my_ones
df_mech_rate_coeffs.k_69 = @. 2.85e-12*exp(-345/df_params.temperature)*my_ones



# save the file
CSV.write("./models/methane/rrate_coeffs_mech.csv", df_mech_rate_coeffs)
