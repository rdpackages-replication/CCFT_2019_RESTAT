clear all
qui do rdbwselect.ado
qui do rdrobust.ado

******* Setup
local n = 1000
set obs `n'
local p = 1
local q = 2
scalar c = 0

local hlist msepob cerpob mserd cerrd
local klist tri uni
local vlist nn hc1 hc2 hc3
local mlist 1

foreach model of local mlist {

disp "model: " `model'

	scalar cz_r0 = .49062457
	scalar cz_r1 = 1.0600066-.45019195
	scalar cz_r2 = 5.7429879-5.5137502
	scalar cz_r3 = 17.142579-20.605129
	scalar cz_r4 = 19.752974-13.317815
	scalar cz_r5 = 7.4750031-10.955869
	scalar cz_l0 = .49062457
	scalar cz_l1 = 1.0600066
	scalar cz_l2 = 5.7429879
	scalar cz_l3 = 17.142579
	scalar cz_l4 = 19.752974
	scalar cz_l5 = 7.4750031
	scalar sigma_y = 0.1295 
	scalar sigma_z = 0.13537

if (`model'==1) {
	scalar c_r0 = 0.52
	scalar c_r1 = 0.84   
	scalar c_r2 = -3 
	scalar c_r3 = 7.99 
	scalar c_r4 = -9.01 
	scalar c_r5 = 3.56
	scalar c_l0 = 0.48 
	scalar c_l1 = 1.27 
	scalar c_l2 = 7.18 
	scalar c_l3 = 20.21 
	scalar c_l4 = 21.54 
	scalar c_l5 = 7.33
	scalar c_rz = 0
	scalar c_lz = 0
	scalar sigma_yz = 0
}
if (`model'==2) {
	scalar c_r0 = .3853995
	scalar c_r1 = .6278497
	scalar c_r2 = -2.846475
	scalar c_r3 = 8.427833
	scalar c_r4 = -10.23893
	scalar c_r5 = 4.31672
	scalar c_l0 = .3650975
	scalar c_l1 = .9634399
	scalar c_l2 = 5.473362
	scalar c_l3 = 15.28258 
	scalar c_l4 = 15.8688
	scalar c_l5 = 5.140786
	scalar c_rz = .2824018 
	scalar c_lz = .2212406
	scalar sigma_yz = 0.2692*sigma_y*sigma_z
}
if (`model'==3) {
	scalar c_r0 = .3853995
	scalar c_r1 = .6278497
	scalar c_r2 = -2.846475
	scalar c_r3 = 8.427833
	scalar c_r4 = -10.23893
	scalar c_r5 = 4.31672
	scalar c_l0 = .3650975
	scalar c_l1 = .9634399
	scalar c_l2 = 5.473362
	scalar c_l3 = 15.28258 
	scalar c_l4 = 15.8688
	scalar c_l5 = 5.140786
	scalar c_rz = .2824018 
	scalar c_lz = .2212406
	scalar sigma_yz = 0
}
if (`model'==4) {
	scalar c_r0 = .3853995
	scalar c_r1 = .6278497
	scalar c_r2 = -2.846475
	scalar c_r3 = 8.427833
	scalar c_r4 = -10.23893
	scalar c_r5 = 4.31672
	scalar c_l0 = .3650975
	scalar c_l1 = .9634399
	scalar c_l2 = 5.473362
	scalar c_l3 = 15.28258 
	scalar c_l4 = 15.8688
	scalar c_l5 = 5.140786
	scalar c_rz = .2824018 
	scalar c_lz = .2212406
	scalar sigma_yz = 2*0.2692*sigma_y*sigma_z
}

local tau_pob = (c_r0 - c_l0) + (c_rz*cz_r0 - c_lz*cz_l0)
		
		foreach bwselect of local hlist {
			foreach kernel   of local klist {
				foreach vce      of local vlist {
					postfile sim_`model'_`bwselect'_`kernel'_`vce'_`n' h h_cov b b_cov tau_pob tau_cl tau_cl_cov tau_bc tau_bc_cov tau_bc_r1 tau_bc_r1_cov ///
													                                         	se_cl  se_cl_cov  se_rb  se_rb_cov  se_rb_r1  se_rb_r1_cov ///
					 using rdsim_`model'_`bwselect'_`kernel'_`vce'_`n', replace
				}
			}
		}
														
							
forvalues i = 1/1000 {
		disp `i'
		tempvar y y_l y_r x z z_l z_r u_y u_z 
		qui gen `x'   = 2*rbeta(2,4)-1 
		matrix C = (sigma_y^2, sigma_yz \ sigma_yz , sigma_z^2)
		drawnorm `u_y' `u_z', n(`n') cov(C)
		qui gen `z_l' = cz_l0 + cz_l1*`x' + cz_l2*`x'^2 + cz_l3*`x'^3 + cz_l4*`x'^4 + cz_l5*`x'^5 + `u_z'
		qui gen `z_r' = cz_r0 + cz_r1*`x' + cz_r2*`x'^2 + cz_r3*`x'^3 + cz_r4*`x'^4 + cz_r5*`x'^5 + `u_z'
		qui gen `y_l' = c_l0 + c_l1*`x' + c_l2*`x'^2 + c_l3*`x'^3 + c_l4*`x'^4 + c_l5*`x'^5 + c_lz*`z_l' + `u_y'
		qui gen `y_r' = c_r0 + c_r1*`x' + c_r2*`x'^2 + c_r3*`x'^3 + c_r4*`x'^4 + c_r5*`x'^5 + c_rz*`z_r' + `u_y'
		qui gen `y'   = .
		qui replace `y' = `y_l' if `x'<c
		qui replace `y' = `y_r' if `x'>=c
		qui gen `z'   = .
		qui replace `z' = `z_l' if `x'<c
		qui replace `z' = `z_r' if `x'>=c

		foreach bwselect of local hlist {
			foreach kernel   of local klist {
				foreach vce      of local vlist {
	
	if (`model'==1) {
		local h_pob_mse     = 0.1655316
		local b_pob_mse     = 0.2511249
		local h_pob_mse_cov = 0.1655316
		local b_pob_mse_cov = 0.2511249
		if (`n'==1000) {
			local h_pob_mse     = 0.1441037 
			local b_pob_mse     = 0.2274497
			local h_pob_mse_cov = 0.1441037
			local b_pob_mse_cov = 0.2274497
		}
		
		if ("`kernel'"=="uni"){
			local h_pob_mse     =  0.1301084
			local b_pob_mse     =  0.2224957
			local h_pob_mse_cov =  0.1301084 
			local b_pob_mse_cov =  0.2224957
			if (`n'==1000) {
				local h_pob_mse     = 0.1132659 
				local b_pob_mse     = 0.2015196
				local h_pob_mse_cov = 0.1132659 
				local b_pob_mse_cov = 0.2015196
			}
		}
	
	
	}
	if (`model'==2) {
		local h_pob_mse     = 0.1794462 
		local b_pob_mse     = 0.2638802
		local h_pob_mse_cov = 0.1809309 
		local b_pob_mse_cov = 0.2668585
		if (`n'==1000) {
			local h_pob_mse     = 0.1562170 
			local b_pob_mse     = 0.2390026
			local h_pob_mse_cov = 0.1575095 
			local b_pob_mse_cov = 0.2417001
		}
	
		if ("`kernel'"=="uni"){
			local h_pob_mse     = 0.1410453 
			local b_pob_mse     = 0.2337969
			local h_pob_mse_cov = 0.1422123 
			local b_pob_mse_cov = 0.2364356
			if (`n'==1000) {
				local h_pob_mse     = 0.1227870 
				local b_pob_mse     = 0.2117554
				local h_pob_mse_cov = 0.1238030 
				local b_pob_mse_cov = 0.2141453
			}
		}
		
	}
	if (`model'==3) {
		local h_pob_mse     = 0.1794462 
		local b_pob_mse     = 0.2638802
		local h_pob_mse_cov = 0.1771628 
		local b_pob_mse_cov = 0.2628768
		if (`n'==1000) {
			 local h_pob_mse     = 0.1562170 
			 local b_pob_mse     = 0.2390026
			 local h_pob_mse_cov = 0.1542292 
			 local b_pob_mse_cov = 0.2380937
		}
		
		if ("`kernel'"=="uni"){
			local h_pob_mse     = 0.1410453 
			local b_pob_mse     = 0.2337969
			local h_pob_mse_cov = 0.1392505 
			local b_pob_mse_cov = 0.2329078
			if (`n'==1000) {
				local h_pob_mse     = 0.1227870 
				local b_pob_mse     = 0.2117554
				local h_pob_mse_cov = 0.1212246 
				local b_pob_mse_cov = 0.2109501
			}
		}
	
	}
	if (`model'==4) {
		local h_pob_mse     = 0.1794462 
		local b_pob_mse     = 0.2638802
		local h_pob_mse_cov = 0.1844091 
		local b_pob_mse_cov = 0.2705128
		if (`n'==1000) {
			local h_pob_mse     = 0.1562170 
			local b_pob_mse     = 0.2390026
			local h_pob_mse_cov = 0.1605374 
			local b_pob_mse_cov = 0.2450098
		}
		
		if ("`kernel'"=="uni"){
			local h_pob_mse     = 0.1410453 
			local b_pob_mse     = 0.2337969
			local h_pob_mse_cov = 0.1449461 
			local b_pob_mse_cov = 0.2396733
			if (`n'==1000) {
				local h_pob_mse     = 0.1227870 
				local b_pob_mse     = 0.2117554
				local h_pob_mse_cov = 0.1261829 
				local b_pob_mse_cov = 0.2170778
			}
		}
	}
	
	if ("`bwselect'"=="msepob") {
		local h     = `h_pob_mse'
		local b     = `b_pob_mse'
		local h_cov = `h_pob_mse_cov'
		local b_cov = `b_pob_mse_cov'
	}
	
	if ("`bwselect'"=="cerpob") {
		scalar cer_h = `n'^(-(`p'/((3+`p')*(3+2*`p'))))
		scalar cer_b = `n'^(-(`q'/((3+`q')*(3+2*`q'))))
		
		local h_pob_cer = `h_pob_mse'*cer_h
		local b_pob_cer = `b_pob_mse'*cer_b
		local h_pob_cer_cov = `h_pob_mse_cov'*cer_h
		local b_pob_cer_cov = `b_pob_mse_cov'*cer_b
	
		local h     = `h_pob_cer'
		local b     = `b_pob_cer'
		local h_cov = `h_pob_cer_cov'
		local b_cov = `b_pob_cer_cov'
	}		
	
	if ("`bwselect'"=="mserd" | "`bwselect'"=="cerrd") {
			qui rdbwselect `y' `x', bwselect(`bwselect') vce(`vce') kernel(`kernel') scaleregul(1)
			local h = e(h_`bwselect')
			local b = e(b_`bwselect')
			qui rdbwselect `y' `x', bwselect(`bwselect') vce(`vce') kernel(`kernel') covs(`z') scaleregul(1)
			local h_cov = e(h_`bwselect')
			local b_cov = e(b_`bwselect')
	}
	
			qui rdrobust `y' `x', h(`h') b(`b') vce(`vce') kernel(`kernel')
			local tau_cl = e(tau_cl)
			local tau_bc = e(tau_bc)
			local se_cl  = e(se_tau_cl)
			local se_rb  = e(se_tau_rb)
			
			qui rdrobust `y' `x', h(`h') b(`h') vce(`vce') kernel(`kernel')
			local tau_cl_r1 = e(tau_cl)
			local tau_bc_r1 = e(tau_bc)
			local se_cl_r1  = e(se_tau_cl)
			local se_rb_r1  = e(se_tau_rb)
			
			qui rdrobust `y' `x', h(`h_cov') b(`b_cov') vce(`vce') kernel(`kernel') covs(`z')
			local tau_cl_cov = e(tau_cl)
			local tau_bc_cov = e(tau_bc)
			local se_cl_cov  = e(se_tau_cl)
			local se_rb_cov  = e(se_tau_rb)
			
			qui rdrobust `y' `x', h(`h_cov') b(`h_cov') vce(`vce') kernel(`kernel') covs(`z')  
			local tau_cl_r1_cov = e(tau_cl)
			local tau_bc_r1_cov = e(tau_bc)
			local se_cl_r1_cov  = e(se_tau_cl)
			local se_rb_r1_cov  = e(se_tau_rb)

			post sim_`model'_`bwselect'_`kernel'_`vce'_`n' (`h')  (`h_cov')  (`b')  (`b_cov') (`tau_pob')  (`tau_cl')  (`tau_cl_cov')  (`tau_bc')  (`tau_bc_cov')  (`tau_bc_r1')  (`tau_bc_r1_cov') ///
																											(`se_cl')   (`se_cl_cov')   (`se_rb')   (`se_rb_cov')   (`se_rb_r1')   (`se_rb_r1_cov')
			}
		}
	}
}

		foreach bwselect of local hlist {
			foreach kernel   of local klist {
				foreach vce      of local vlist {
					postclose sim_`model'_`bwselect'_`kernel'_`vce'_`n'
				}
			}
		}			

		foreach bwselect of local hlist {
			foreach kernel   of local klist {
				foreach vce      of local vlist {
					use rdsim_`model'_`bwselect'_`kernel'_`vce'_`n', clear
					scalar quant=-invnormal(abs((1-(95/100))/2))
					gen ec_cl    = abs((tau_cl-tau_pob)/se_cl)<=quant 
					gen ec_rb    = abs((tau_bc-tau_pob)/se_rb)<=quant
					gen ec_rb_r1 = abs((tau_bc_r1-tau_pob)/se_rb_r1)<=quant
					gen il_cl    = 2*quant*se_cl
					gen il_rb    = 2*quant*se_rb
					gen il_rb_r1 = 2*quant*se_rb_r1
					gen ec_cl_cov    = abs((tau_cl_cov-tau_pob)/se_cl_cov)<=quant 
					gen ec_rb_cov    = abs((tau_bc_cov-tau_pob)/se_rb_cov)<=quant
					gen ec_rb_r1_cov = abs((tau_bc_r1_cov-tau_pob)/se_rb_r1_cov)<=quant
					gen il_cl_cov    = 2*quant*se_cl_cov
					gen il_rb_cov    = 2*quant*se_rb_cov
					gen il_rb_r1_cov = 2*quant*se_rb_r1_cov
					saveold rdsim_`model'_`bwselect'_`kernel'_`vce'_`n', replace
				}
			}
		}
}

*cap log close


