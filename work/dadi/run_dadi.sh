#!/usr/bin/env bash

export OPTIM=100

#for pop_pair in cc_ei cc_lo cc_cm ei_lo ei_cm cm_lo
#do
infer_dm() {
	pop_pair=$1
	# Infer models
#	dadi-cli InferDM --fs All.$pop_pair.folded.downproj.fs --model no_mig --cpus 32 --output-prefix All.$pop_pair.folded.downproj.no_mig --lbounds 1e-3 1e-3 0 --ubounds 2 2 10 --global-optimization --optimizations $OPTIM
#	dadi-cli InferDM --fs All.$pop_pair.folded.downproj.fs --model split_mig --cpus 32 --output-prefix All.$pop_pair.folded.downproj.split_mig --lbounds 1e-3 1e-3 0 0 --ubounds 2 2 10 5 --global-optimization --optimizations $OPTIM
#	dadi-cli InferDM --fs All.$pop_pair.folded.downproj.fs --model sec_contact_sym_mig --cpus 32 --output-prefix All.$pop_pair.folded.downproj.sec_contact_sym_mig --lbounds 1e-3 1e-3 0 0 0 --ubounds 2 2 5 10 10 --global-optimization --optimizations $OPTIM
#	dadi-cli InferDM --fs All.$pop_pair.folded.downproj.fs --model asym_mig --cpus 32 --output-prefix All.$pop_pair.folded.downproj.asym_mig --lbounds 1e-3 1e-3 0 0 0 --ubounds 2 2 10 5 5 --global-optimization --optimizations $OPTIM
#	dadi-cli InferDM --fs All.$pop_pair.folded.downproj.fs --model sec_contact_asym_mig --cpus 32 --output-prefix All.$pop_pair.folded.downproj.sec_contact_asym_mig --lbounds 1e-3 1e-3 0 0 0 0 --ubounds 2 2 5 5 10 10 --global-optimization --optimizations $OPTIM

	# Get best fit runs
	dadi-cli BestFit --input-prefix All.$pop_pair.folded.downproj.no_mig.InferDM --lbounds 1e-3 1e-3 0 --ubounds 2 2 10
	dadi-cli BestFit --input-prefix All.$pop_pair.folded.downproj.split_mig.InferDM --lbounds 1e-3 1e-3 0 0 --ubounds 2 2 10 5
	dadi-cli BestFit --input-prefix All.$pop_pair.folded.downproj.sec_contact_sym_mig.InferDM --lbounds 1e-3 1e-3 0 0 0 --ubounds 2 2 5 10 10
	dadi-cli BestFit --input-prefix All.$pop_pair.folded.downproj.asym_mig.InferDM --lbounds 1e-3 1e-3 0 0 0 --ubounds 2 2 10 5 5
	dadi-cli BestFit --input-prefix All.$pop_pair.folded.downproj.sec_contact_asym_mig.InferDM --lbounds 1e-3 1e-3 0 0 0 0 --ubounds 2 2 5 5 10 10

	# Plot SFS vs model
	vmin=1e-4
	if [ $pop_pair = 'cc_ei' ]
	then
		vmin=1e-4
	else
		vmin=1e-2
	fi
	dadi-cli Plot --fs All.$pop_pair.folded.downproj.fs --model no_mig --demo-popt All.$pop_pair.folded.downproj.no_mig.InferDM.bestfits --output All.$pop_pair.folded.downproj.no_mig.pdf --vmin $vmin
	dadi-cli Plot --fs All.$pop_pair.folded.downproj.fs --model split_mig --demo-popt All.$pop_pair.folded.downproj.split_mig.InferDM.bestfits --output All.$pop_pair.folded.downproj.split_mig.pdf --vmin $vmin
	dadi-cli Plot --fs All.$pop_pair.folded.downproj.fs --model sec_contact_sym_mig --demo-popt All.$pop_pair.folded.downproj.sec_contact_sym_mig.InferDM.bestfits --output All.$pop_pair.folded.downproj.sec_contact_sym_mig.pdf --vmin $vmin
	dadi-cli Plot --fs All.$pop_pair.folded.downproj.fs --model asym_mig --demo-popt All.$pop_pair.folded.downproj.asym_mig.InferDM.bestfits --output All.$pop_pair.folded.downproj.asym_mig.pdf --vmin $vmin
	dadi-cli Plot --fs All.$pop_pair.folded.downproj.fs --model sec_contact_asym_mig --demo-popt All.$pop_pair.folded.downproj.sec_contact_asym_mig.InferDM.bestfits --output All.$pop_pair.folded.downproj.sec_contact_asym_mig.pdf --vmin $vmin

	# Calculate confidence intervals
	dadi-cli StatDM --fs All.$pop_pair.folded.downproj.fs --model no_mig --demo-popt All.$pop_pair.folded.downproj.no_mig.InferDM.bestfits --bootstrapping-dir "${pop_pair}_bootstrap/" --output All.$pop_pair.folded.downproj.no_mig.godambe.ci
	dadi-cli StatDM --fs All.$pop_pair.folded.downproj.fs --model split_mig --demo-popt All.$pop_pair.folded.downproj.split_mig.InferDM.bestfits --bootstrapping-dir "${pop_pair}_bootstrap/" --output All.$pop_pair.folded.downproj.split_mig.godambe.ci
	dadi-cli StatDM --fs All.$pop_pair.folded.downproj.fs --model sec_contact_sym_mig --demo-popt All.$pop_pair.folded.downproj.sec_contact_sym_mig.InferDM.bestfits --bootstrapping-dir "${pop_pair}_bootstrap/" --output All.$pop_pair.folded.downproj.sec_contact_sym_mig.godambe.ci
	dadi-cli StatDM --fs All.$pop_pair.folded.downproj.fs --model asym_mig --demo-popt All.$pop_pair.folded.downproj.asym_mig.InferDM.bestfits --bootstrapping-dir "${pop_pair}_bootstrap/" --output All.$pop_pair.folded.downproj.asym_mig.godambe.ci
	dadi-cli StatDM --fs All.$pop_pair.folded.downproj.fs --model sec_contact_asym_mig --demo-popt All.$pop_pair.folded.downproj.sec_contact_asym_mig.InferDM.bestfits --bootstrapping-dir "${pop_pair}_bootstrap/" --output All.$pop_pair.folded.downproj.sec_contact_asym_mig.godambe.ci
}
export -f infer_dm
parallel -j 6 infer_dm {} ::: cc_ei cc_lo cc_cm ei_lo ei_cm cm_lo
#done
