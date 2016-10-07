################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/daycent/Calcpet.c \
../src/daycent/Cmpnfrac.c \
../src/daycent/Daylen.c \
../src/daycent/Diffusiv.c \
../src/daycent/Ferr.c \
../src/daycent/Floclr.c \
../src/daycent/Floclr_double.c \
../src/daycent/Floclr_double_in.c \
../src/daycent/Floclr_double_out.c \
../src/daycent/Flow.c \
../src/daycent/Flow_double.c \
../src/daycent/Flow_double_in.c \
../src/daycent/Flow_double_out.c \
../src/daycent/Flowd_double.c \
../src/daycent/Flowd_double_in.c \
../src/daycent/Flowd_double_out.c \
../src/daycent/Flowup.c \
../src/daycent/Flowup_double.c \
../src/daycent/Flowup_double_in.c \
../src/daycent/Flowup_double_out.c \
../src/daycent/Getdiff.c \
../src/daycent/H2oflux.c \
../src/daycent/HYDRUS_Init.c \
../src/daycent/HYDRUS_soilwater.c \
../src/daycent/Hwdrain.c \
../src/daycent/Initlyrs.c \
../src/daycent/Initsw.c \
../src/daycent/Leachdly.c \
../src/daycent/Nitrify.c \
../src/daycent/Petrad.c \
../src/daycent/Pevapdly.c \
../src/daycent/Pi_funcs.c \
../src/daycent/Potbse.c \
../src/daycent/Potbst.c \
../src/daycent/Pteevap.c \
../src/daycent/Rainflux.c \
../src/daycent/Setamov.c \
../src/daycent/Setasmos.c \
../src/daycent/Setlyrs.c \
../src/daycent/Showlyrs.c \
../src/daycent/Snowcent.c \
../src/daycent/Soiltemp.c \
../src/daycent/Svapor.c \
../src/daycent/Tanfunc.c \
../src/daycent/Therm.c \
../src/daycent/Trwtavg.c \
../src/daycent/Watrate.c \
../src/daycent/Watrbal.c \
../src/daycent/Watreqn.c \
../src/daycent/Watrflow_Hydrus.c \
../src/daycent/Watrlit.c \
../src/daycent/Watrstcr.c \
../src/daycent/Wfps.c \
../src/daycent/Wrtbio.c \
../src/daycent/Wrtsoiln.c \
../src/daycent/Wrtstemp.c \
../src/daycent/Wrtswc.c \
../src/daycent/Wrtwfps.c \
../src/daycent/balanceN.c \
../src/daycent/calcdefac.c \
../src/daycent/denitrify.c \
../src/daycent/flowd.c \
../src/daycent/fracbslos.c \
../src/daycent/getsoilprop.c \
../src/daycent/initdaily.c \
../src/daycent/initsite_tg.c \
../src/daycent/litstcrev.c \
../src/daycent/methane.c \
../src/daycent/nox_pulse.c \
../src/daycent/showminrl.c \
../src/daycent/snowmodel.c \
../src/daycent/soiltransp.c \
../src/daycent/swpotentl.c \
../src/daycent/tgmodel.c \
../src/daycent/updateN.c \
../src/daycent/wrtco2.c \
../src/daycent/wrtdeadc.c \
../src/daycent/wrtlivec.c \
../src/daycent/wrtmresp.c \
../src/daycent/wrtsoilc.c \
../src/daycent/wrtsysc.c \
../src/daycent/wrtwflux.c \
../src/daycent/wrtyearsum.c 

F_SRCS += \
../src/daycent/Adjlig.f \
../src/daycent/Agdrat.f \
../src/daycent/Anerob.f \
../src/daycent/Annacc.f \
../src/daycent/Bgdrat.f \
../src/daycent/Calciv.f \
../src/daycent/Candec.f \
../src/daycent/Catanf.f \
../src/daycent/Ckdata.f \
../src/daycent/Cmplig.f \
../src/daycent/Co2eff.f \
../src/daycent/Crop.f \
../src/daycent/Cropin.f \
../src/daycent/Csa_main.f \
../src/daycent/Csched.f \
../src/daycent/Cultin.f \
../src/daycent/Cultiv.f \
../src/daycent/Cutrtn.f \
../src/daycent/Cycle.f \
../src/daycent/Declig.f \
../src/daycent/Decomp.f \
../src/daycent/Dedrem.f \
../src/daycent/Droot.f \
../src/daycent/Dshoot.f \
../src/daycent/Eachyr.f \
../src/daycent/Erosn.f \
../src/daycent/Esched.f \
../src/daycent/Extend.f \
../src/daycent/Falstd.f \
../src/daycent/Fertin.f \
../src/daycent/Firein.f \
../src/daycent/Firrtn.f \
../src/daycent/Fixin.f \
../src/daycent/Fracis.f \
../src/daycent/Frem.f \
../src/daycent/Frespr.f \
../src/daycent/Fsfunc.f \
../src/daycent/Getwth.f \
../src/daycent/Gpdf.f \
../src/daycent/Grazin.f \
../src/daycent/Grem.f \
../src/daycent/Grochk.f \
../src/daycent/Growth.f \
../src/daycent/Harvin.f \
../src/daycent/Harvst.f \
../src/daycent/Inprac.f \
../src/daycent/Irrgin.f \
../src/daycent/Irrigt.f \
../src/daycent/Killrt.f \
../src/daycent/Lacalc.f \
../src/daycent/Laprod.f \
../src/daycent/Line.f \
../src/daycent/Litburn.f \
../src/daycent/Litdec.f \
../src/daycent/Livrem.f \
../src/daycent/Message.f \
../src/daycent/Mnracc.f \
../src/daycent/Nutrlm.f \
../src/daycent/Omadin.f \
../src/daycent/Ozonein.f \
../src/daycent/Parse.f \
../src/daycent/Partit.f \
../src/daycent/Potcrp.f \
../src/daycent/Potfor.f \
../src/daycent/Potprod.f \
../src/daycent/Pprdwc.f \
../src/daycent/Prcgrw.f \
../src/daycent/Predec.f \
../src/daycent/Prelim.f \
../src/daycent/Pschem.f \
../src/daycent/Ramp.f \
../src/daycent/Readblk.f \
../src/daycent/Respir.f \
../src/daycent/Restrp.f \
../src/daycent/Rtimp.f \
../src/daycent/Savarp.f \
../src/daycent/Schedl.f \
../src/daycent/Simsom.f \
../src/daycent/Sitein.f \
../src/daycent/Sitein_grid.f \
../src/daycent/Soilos.f \
../src/daycent/Somdec.f \
../src/daycent/Sumcar.f \
../src/daycent/Surftemp.f \
../src/daycent/Treein.f \
../src/daycent/Trees.f \
../src/daycent/Tremin.f \
../src/daycent/Wdeath.f \
../src/daycent/Woodec.f \
../src/daycent/Wrtbin.f \
../src/daycent/Wthini.f \
../src/daycent/adjustpar.f \
../src/daycent/cropDynC.f \
../src/daycent/csa_detiv.f \
../src/daycent/dailymoist_hydrus.f \
../src/daycent/default.f \
../src/daycent/fltce.f \
../src/daycent/fozone_forest.f \
../src/daycent/froota.f \
../src/daycent/grazrst.f \
../src/daycent/initialize.f \
../src/daycent/leafa.f \
../src/daycent/phshift.f \
../src/daycent/shwave.f \
../src/daycent/tcalc.f \
../src/daycent/treeDynC.f \
../src/daycent/treegrow.f 

OBJS += \
./src/daycent/Adjlig.o \
./src/daycent/Agdrat.o \
./src/daycent/Anerob.o \
./src/daycent/Annacc.o \
./src/daycent/Bgdrat.o \
./src/daycent/Calciv.o \
./src/daycent/Calcpet.o \
./src/daycent/Candec.o \
./src/daycent/Catanf.o \
./src/daycent/Ckdata.o \
./src/daycent/Cmplig.o \
./src/daycent/Cmpnfrac.o \
./src/daycent/Co2eff.o \
./src/daycent/Crop.o \
./src/daycent/Cropin.o \
./src/daycent/Csa_main.o \
./src/daycent/Csched.o \
./src/daycent/Cultin.o \
./src/daycent/Cultiv.o \
./src/daycent/Cutrtn.o \
./src/daycent/Cycle.o \
./src/daycent/Daylen.o \
./src/daycent/Declig.o \
./src/daycent/Decomp.o \
./src/daycent/Dedrem.o \
./src/daycent/Diffusiv.o \
./src/daycent/Droot.o \
./src/daycent/Dshoot.o \
./src/daycent/Eachyr.o \
./src/daycent/Erosn.o \
./src/daycent/Esched.o \
./src/daycent/Extend.o \
./src/daycent/Falstd.o \
./src/daycent/Ferr.o \
./src/daycent/Fertin.o \
./src/daycent/Firein.o \
./src/daycent/Firrtn.o \
./src/daycent/Fixin.o \
./src/daycent/Floclr.o \
./src/daycent/Floclr_double.o \
./src/daycent/Floclr_double_in.o \
./src/daycent/Floclr_double_out.o \
./src/daycent/Flow.o \
./src/daycent/Flow_double.o \
./src/daycent/Flow_double_in.o \
./src/daycent/Flow_double_out.o \
./src/daycent/Flowd_double.o \
./src/daycent/Flowd_double_in.o \
./src/daycent/Flowd_double_out.o \
./src/daycent/Flowup.o \
./src/daycent/Flowup_double.o \
./src/daycent/Flowup_double_in.o \
./src/daycent/Flowup_double_out.o \
./src/daycent/Fracis.o \
./src/daycent/Frem.o \
./src/daycent/Frespr.o \
./src/daycent/Fsfunc.o \
./src/daycent/Getdiff.o \
./src/daycent/Getwth.o \
./src/daycent/Gpdf.o \
./src/daycent/Grazin.o \
./src/daycent/Grem.o \
./src/daycent/Grochk.o \
./src/daycent/Growth.o \
./src/daycent/H2oflux.o \
./src/daycent/HYDRUS_Init.o \
./src/daycent/HYDRUS_soilwater.o \
./src/daycent/Harvin.o \
./src/daycent/Harvst.o \
./src/daycent/Hwdrain.o \
./src/daycent/Initlyrs.o \
./src/daycent/Initsw.o \
./src/daycent/Inprac.o \
./src/daycent/Irrgin.o \
./src/daycent/Irrigt.o \
./src/daycent/Killrt.o \
./src/daycent/Lacalc.o \
./src/daycent/Laprod.o \
./src/daycent/Leachdly.o \
./src/daycent/Line.o \
./src/daycent/Litburn.o \
./src/daycent/Litdec.o \
./src/daycent/Livrem.o \
./src/daycent/Message.o \
./src/daycent/Mnracc.o \
./src/daycent/Nitrify.o \
./src/daycent/Nutrlm.o \
./src/daycent/Omadin.o \
./src/daycent/Ozonein.o \
./src/daycent/Parse.o \
./src/daycent/Partit.o \
./src/daycent/Petrad.o \
./src/daycent/Pevapdly.o \
./src/daycent/Pi_funcs.o \
./src/daycent/Potbse.o \
./src/daycent/Potbst.o \
./src/daycent/Potcrp.o \
./src/daycent/Potfor.o \
./src/daycent/Potprod.o \
./src/daycent/Pprdwc.o \
./src/daycent/Prcgrw.o \
./src/daycent/Predec.o \
./src/daycent/Prelim.o \
./src/daycent/Pschem.o \
./src/daycent/Pteevap.o \
./src/daycent/Rainflux.o \
./src/daycent/Ramp.o \
./src/daycent/Readblk.o \
./src/daycent/Respir.o \
./src/daycent/Restrp.o \
./src/daycent/Rtimp.o \
./src/daycent/Savarp.o \
./src/daycent/Schedl.o \
./src/daycent/Setamov.o \
./src/daycent/Setasmos.o \
./src/daycent/Setlyrs.o \
./src/daycent/Showlyrs.o \
./src/daycent/Simsom.o \
./src/daycent/Sitein.o \
./src/daycent/Sitein_grid.o \
./src/daycent/Snowcent.o \
./src/daycent/Soilos.o \
./src/daycent/Soiltemp.o \
./src/daycent/Somdec.o \
./src/daycent/Sumcar.o \
./src/daycent/Surftemp.o \
./src/daycent/Svapor.o \
./src/daycent/Tanfunc.o \
./src/daycent/Therm.o \
./src/daycent/Treein.o \
./src/daycent/Trees.o \
./src/daycent/Tremin.o \
./src/daycent/Trwtavg.o \
./src/daycent/Watrate.o \
./src/daycent/Watrbal.o \
./src/daycent/Watreqn.o \
./src/daycent/Watrflow_Hydrus.o \
./src/daycent/Watrlit.o \
./src/daycent/Watrstcr.o \
./src/daycent/Wdeath.o \
./src/daycent/Wfps.o \
./src/daycent/Woodec.o \
./src/daycent/Wrtbin.o \
./src/daycent/Wrtbio.o \
./src/daycent/Wrtsoiln.o \
./src/daycent/Wrtstemp.o \
./src/daycent/Wrtswc.o \
./src/daycent/Wrtwfps.o \
./src/daycent/Wthini.o \
./src/daycent/adjustpar.o \
./src/daycent/balanceN.o \
./src/daycent/calcdefac.o \
./src/daycent/cropDynC.o \
./src/daycent/csa_detiv.o \
./src/daycent/dailymoist_hydrus.o \
./src/daycent/default.o \
./src/daycent/denitrify.o \
./src/daycent/flowd.o \
./src/daycent/fltce.o \
./src/daycent/fozone_forest.o \
./src/daycent/fracbslos.o \
./src/daycent/froota.o \
./src/daycent/getsoilprop.o \
./src/daycent/grazrst.o \
./src/daycent/initdaily.o \
./src/daycent/initialize.o \
./src/daycent/initsite_tg.o \
./src/daycent/leafa.o \
./src/daycent/litstcrev.o \
./src/daycent/methane.o \
./src/daycent/nox_pulse.o \
./src/daycent/phshift.o \
./src/daycent/showminrl.o \
./src/daycent/shwave.o \
./src/daycent/snowmodel.o \
./src/daycent/soiltransp.o \
./src/daycent/swpotentl.o \
./src/daycent/tcalc.o \
./src/daycent/tgmodel.o \
./src/daycent/treeDynC.o \
./src/daycent/treegrow.o \
./src/daycent/updateN.o \
./src/daycent/wrtco2.o \
./src/daycent/wrtdeadc.o \
./src/daycent/wrtlivec.o \
./src/daycent/wrtmresp.o \
./src/daycent/wrtsoilc.o \
./src/daycent/wrtsysc.o \
./src/daycent/wrtwflux.o \
./src/daycent/wrtyearsum.o 

C_DEPS += \
./src/daycent/Calcpet.d \
./src/daycent/Cmpnfrac.d \
./src/daycent/Daylen.d \
./src/daycent/Diffusiv.d \
./src/daycent/Ferr.d \
./src/daycent/Floclr.d \
./src/daycent/Floclr_double.d \
./src/daycent/Floclr_double_in.d \
./src/daycent/Floclr_double_out.d \
./src/daycent/Flow.d \
./src/daycent/Flow_double.d \
./src/daycent/Flow_double_in.d \
./src/daycent/Flow_double_out.d \
./src/daycent/Flowd_double.d \
./src/daycent/Flowd_double_in.d \
./src/daycent/Flowd_double_out.d \
./src/daycent/Flowup.d \
./src/daycent/Flowup_double.d \
./src/daycent/Flowup_double_in.d \
./src/daycent/Flowup_double_out.d \
./src/daycent/Getdiff.d \
./src/daycent/H2oflux.d \
./src/daycent/HYDRUS_Init.d \
./src/daycent/HYDRUS_soilwater.d \
./src/daycent/Hwdrain.d \
./src/daycent/Initlyrs.d \
./src/daycent/Initsw.d \
./src/daycent/Leachdly.d \
./src/daycent/Nitrify.d \
./src/daycent/Petrad.d \
./src/daycent/Pevapdly.d \
./src/daycent/Pi_funcs.d \
./src/daycent/Potbse.d \
./src/daycent/Potbst.d \
./src/daycent/Pteevap.d \
./src/daycent/Rainflux.d \
./src/daycent/Setamov.d \
./src/daycent/Setasmos.d \
./src/daycent/Setlyrs.d \
./src/daycent/Showlyrs.d \
./src/daycent/Snowcent.d \
./src/daycent/Soiltemp.d \
./src/daycent/Svapor.d \
./src/daycent/Tanfunc.d \
./src/daycent/Therm.d \
./src/daycent/Trwtavg.d \
./src/daycent/Watrate.d \
./src/daycent/Watrbal.d \
./src/daycent/Watreqn.d \
./src/daycent/Watrflow_Hydrus.d \
./src/daycent/Watrlit.d \
./src/daycent/Watrstcr.d \
./src/daycent/Wfps.d \
./src/daycent/Wrtbio.d \
./src/daycent/Wrtsoiln.d \
./src/daycent/Wrtstemp.d \
./src/daycent/Wrtswc.d \
./src/daycent/Wrtwfps.d \
./src/daycent/balanceN.d \
./src/daycent/calcdefac.d \
./src/daycent/denitrify.d \
./src/daycent/flowd.d \
./src/daycent/fracbslos.d \
./src/daycent/getsoilprop.d \
./src/daycent/initdaily.d \
./src/daycent/initsite_tg.d \
./src/daycent/litstcrev.d \
./src/daycent/methane.d \
./src/daycent/nox_pulse.d \
./src/daycent/showminrl.d \
./src/daycent/snowmodel.d \
./src/daycent/soiltransp.d \
./src/daycent/swpotentl.d \
./src/daycent/tgmodel.d \
./src/daycent/updateN.d \
./src/daycent/wrtco2.d \
./src/daycent/wrtdeadc.d \
./src/daycent/wrtlivec.d \
./src/daycent/wrtmresp.d \
./src/daycent/wrtsoilc.d \
./src/daycent/wrtsysc.d \
./src/daycent/wrtwflux.d \
./src/daycent/wrtyearsum.d 


# Each subdirectory must supply rules for building sources it contributes
src/daycent/%.o: ../src/daycent/%.f
	@echo 'Building file: $<'
	@echo 'Invoking: GNU Fortran Compiler'
	$(FC) -fno-underscoring -O0 -g -Wall -c -fmessage-length=0 -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/daycent/Adjlig.o: ../src/daycent/Adjlig.f

src/daycent/Agdrat.o: ../src/daycent/Agdrat.f

src/daycent/Anerob.o: ../src/daycent/Anerob.f

src/daycent/Annacc.o: ../src/daycent/Annacc.f

src/daycent/Bgdrat.o: ../src/daycent/Bgdrat.f

src/daycent/Calciv.o: ../src/daycent/Calciv.f

src/daycent/%.o: ../src/daycent/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	$(CC) -O0 -g3 -Wall -c -fmessage-length=0 -fPIC -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/daycent/Candec.o: ../src/daycent/Candec.f

src/daycent/Catanf.o: ../src/daycent/Catanf.f

src/daycent/Ckdata.o: ../src/daycent/Ckdata.f

src/daycent/Cmplig.o: ../src/daycent/Cmplig.f

src/daycent/Co2eff.o: ../src/daycent/Co2eff.f

src/daycent/Crop.o: ../src/daycent/Crop.f

src/daycent/Cropin.o: ../src/daycent/Cropin.f

src/daycent/Csa_main.o: ../src/daycent/Csa_main.f

src/daycent/Csched.o: ../src/daycent/Csched.f

src/daycent/Cultin.o: ../src/daycent/Cultin.f

src/daycent/Cultiv.o: ../src/daycent/Cultiv.f

src/daycent/Cutrtn.o: ../src/daycent/Cutrtn.f

src/daycent/Cycle.o: ../src/daycent/Cycle.f

src/daycent/Declig.o: ../src/daycent/Declig.f

src/daycent/Decomp.o: ../src/daycent/Decomp.f

src/daycent/Dedrem.o: ../src/daycent/Dedrem.f

src/daycent/Droot.o: ../src/daycent/Droot.f

src/daycent/Dshoot.o: ../src/daycent/Dshoot.f

src/daycent/Eachyr.o: ../src/daycent/Eachyr.f

src/daycent/Erosn.o: ../src/daycent/Erosn.f

src/daycent/Esched.o: ../src/daycent/Esched.f

src/daycent/Extend.o: ../src/daycent/Extend.f

src/daycent/Falstd.o: ../src/daycent/Falstd.f

src/daycent/Fertin.o: ../src/daycent/Fertin.f

src/daycent/Firein.o: ../src/daycent/Firein.f

src/daycent/Firrtn.o: ../src/daycent/Firrtn.f

src/daycent/Fixin.o: ../src/daycent/Fixin.f

src/daycent/Fracis.o: ../src/daycent/Fracis.f

src/daycent/Frem.o: ../src/daycent/Frem.f

src/daycent/Frespr.o: ../src/daycent/Frespr.f

src/daycent/Fsfunc.o: ../src/daycent/Fsfunc.f

src/daycent/Getwth.o: ../src/daycent/Getwth.f

src/daycent/Gpdf.o: ../src/daycent/Gpdf.f

src/daycent/Grazin.o: ../src/daycent/Grazin.f

src/daycent/Grem.o: ../src/daycent/Grem.f

src/daycent/Grochk.o: ../src/daycent/Grochk.f

src/daycent/Growth.o: ../src/daycent/Growth.f

src/daycent/Harvin.o: ../src/daycent/Harvin.f

src/daycent/Harvst.o: ../src/daycent/Harvst.f

src/daycent/Inprac.o: ../src/daycent/Inprac.f

src/daycent/Irrgin.o: ../src/daycent/Irrgin.f

src/daycent/Irrigt.o: ../src/daycent/Irrigt.f

src/daycent/Killrt.o: ../src/daycent/Killrt.f

src/daycent/Lacalc.o: ../src/daycent/Lacalc.f

src/daycent/Laprod.o: ../src/daycent/Laprod.f

src/daycent/Line.o: ../src/daycent/Line.f

src/daycent/Litburn.o: ../src/daycent/Litburn.f

src/daycent/Litdec.o: ../src/daycent/Litdec.f

src/daycent/Livrem.o: ../src/daycent/Livrem.f

src/daycent/Message.o: ../src/daycent/Message.f

src/daycent/Mnracc.o: ../src/daycent/Mnracc.f

src/daycent/Nutrlm.o: ../src/daycent/Nutrlm.f

src/daycent/Omadin.o: ../src/daycent/Omadin.f

src/daycent/Ozonein.o: ../src/daycent/Ozonein.f

src/daycent/Parse.o: ../src/daycent/Parse.f

src/daycent/Partit.o: ../src/daycent/Partit.f

src/daycent/Potcrp.o: ../src/daycent/Potcrp.f

src/daycent/Potfor.o: ../src/daycent/Potfor.f

src/daycent/Potprod.o: ../src/daycent/Potprod.f

src/daycent/Pprdwc.o: ../src/daycent/Pprdwc.f

src/daycent/Prcgrw.o: ../src/daycent/Prcgrw.f

src/daycent/Predec.o: ../src/daycent/Predec.f

src/daycent/Prelim.o: ../src/daycent/Prelim.f

src/daycent/Pschem.o: ../src/daycent/Pschem.f

src/daycent/Ramp.o: ../src/daycent/Ramp.f

src/daycent/Readblk.o: ../src/daycent/Readblk.f

src/daycent/Respir.o: ../src/daycent/Respir.f

src/daycent/Restrp.o: ../src/daycent/Restrp.f

src/daycent/Rtimp.o: ../src/daycent/Rtimp.f

src/daycent/Savarp.o: ../src/daycent/Savarp.f

src/daycent/Schedl.o: ../src/daycent/Schedl.f

src/daycent/Simsom.o: ../src/daycent/Simsom.f

src/daycent/Sitein.o: ../src/daycent/Sitein.f

src/daycent/Sitein_grid.o: ../src/daycent/Sitein_grid.f

src/daycent/Soilos.o: ../src/daycent/Soilos.f

src/daycent/Somdec.o: ../src/daycent/Somdec.f

src/daycent/Sumcar.o: ../src/daycent/Sumcar.f

src/daycent/Surftemp.o: ../src/daycent/Surftemp.f

src/daycent/Treein.o: ../src/daycent/Treein.f

src/daycent/Trees.o: ../src/daycent/Trees.f

src/daycent/Tremin.o: ../src/daycent/Tremin.f

src/daycent/Wdeath.o: ../src/daycent/Wdeath.f

src/daycent/Woodec.o: ../src/daycent/Woodec.f

src/daycent/Wrtbin.o: ../src/daycent/Wrtbin.f

src/daycent/Wthini.o: ../src/daycent/Wthini.f

src/daycent/adjustpar.o: ../src/daycent/adjustpar.f

src/daycent/cropDynC.o: ../src/daycent/cropDynC.f

src/daycent/csa_detiv.o: ../src/daycent/csa_detiv.f

src/daycent/dailymoist_hydrus.o: ../src/daycent/dailymoist_hydrus.f

src/daycent/default.o: ../src/daycent/default.f

src/daycent/fltce.o: ../src/daycent/fltce.f

src/daycent/fozone_forest.o: ../src/daycent/fozone_forest.f

src/daycent/froota.o: ../src/daycent/froota.f

src/daycent/grazrst.o: ../src/daycent/grazrst.f

src/daycent/initialize.o: ../src/daycent/initialize.f

src/daycent/leafa.o: ../src/daycent/leafa.f

src/daycent/phshift.o: ../src/daycent/phshift.f

src/daycent/shwave.o: ../src/daycent/shwave.f

src/daycent/tcalc.o: ../src/daycent/tcalc.f

src/daycent/treeDynC.o: ../src/daycent/treeDynC.f

src/daycent/treegrow.o: ../src/daycent/treegrow.f


