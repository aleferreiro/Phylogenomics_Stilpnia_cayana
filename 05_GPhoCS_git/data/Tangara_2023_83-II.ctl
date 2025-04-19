GENERAL-INFO-START

	seq-file Tangara_2023_83-II.gphocs
	trace-file Tangara_2023_83-II_TopoRAxML_run1.mcmc
	num-loci 13044  


	burn-in 0
	mcmc-iterations 5000000
	mcmc-sample-skip 100
	iterations-per-log  100
	logs-per-line       100

	#****** Priors **************************
	tau-theta-alpha 1.0
	tau-theta-beta 10000
	tau-theta-print 10000

	mig-rate-alpha 0.001
	mig-rate-beta 0.00001
	mig-rate-print 1.0

	find-finetunes TRUE
	#find-finetunes-num-steps	10000

	locus-mut-rate CONST
GENERAL-INFO-END 
##############################################
CURRENT-POPS-START

	POP-START
		name Cayana
		#samples UFG5654 d UFG5672 d UFG5750 d UFG5762 d UFG5661 d UFG5664 d UFG4007 d UFG4125 d
		samples UFG5654 d UFG5672 d UFG5661 d  UFG4007 d 
	POP-END

	POP-START
		name Marajo
		samples UFG4222 d UFG4228 d UFG4233 d
	POP-END

	#POP-START
	#	name Caatinga
	#	samples UFG5635 d UFG5640 d UFG5648 d UFG5649 d UFG5669 d UFG5698 d UFG5765 d UFG5670 d
	#POP-END

	#POP-START
	#	name EastCerrado
	#	samples Piri253 d UFG2997 d UFG3150 d UFG5131 d UFG5636 d UFG5639 d UFG5646 d UFG5913 d
	#POP-END

	#POP-START
	#	name WestCerrado
	#	samples UFG5642 d UFG5691 d UFG5694 d UFG5756 d UFG5757 d UFG5758 d UFG5759 d UFG5760 d
	#POP-END



POP-START
		name Caatinga
		samples UFG5635 d UFG5640 d UFG5648 d UFG5649 d
	POP-END

	POP-START
		name EastCerrado
		samples Piri253 d UFG2997 d UFG3150 d UFG5131 d 
	POP-END

	POP-START
		name WestCerrado
		samples UFG5642 d UFG5691 d UFG5694 d UFG5756 d 
	POP-END


CURRENT-POPS-END
##############################################
ANCESTRAL-POPS-START

	POP-START
		name CAA-EC
		children Caatinga EastCerrado
	POP-END

	POP-START
		name CAA-EC-WC
		children CAA-EC WestCerrado
	POP-END

	POP-START
		name CAY-MJ
		children Cayana Marajo
	POP-END

	POP-START
		name root
		children  CAY-MJ CAA-EC-WC
	POP-END

ANCESTRAL-POPS-END
##############################################
MIG-BANDS-START
## Cayana
	BAND-START
		source Cayana
		target Marajo
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Marajo 
		target Cayana
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Cayana
		target Caatinga
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Caatinga 
		target Cayana
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Cayana
		target EastCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source EastCerrado 
		target Cayana
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Cayana
		target WestCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source WestCerrado 
		target Cayana
		mig-rate-print 1.0
	BAND-END
# Marajo
	BAND-START
		source Marajo
		target Caatinga
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Caatinga 
		target Marajo
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Marajo
		target EastCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source EastCerrado 
		target Marajo
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source Marajo
		target WestCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source WestCerrado 
		target Marajo
		mig-rate-print 1.0
	BAND-END
# Caatinga
	BAND-START
		source Caatinga 
		target EastCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source EastCerrado
		target Caatinga
		mig-rate-print 1.0
	BAND-END
	BAND-START
		source Caatinga 
		target WestCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source WestCerrado
		target Caatinga
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source EastCerrado 
		target WestCerrado
		mig-rate-print 1.0
	BAND-END

	BAND-START
		source WestCerrado
		target EastCerrado
		mig-rate-print 1.0
	BAND-END
# Ancestral pops 1-3
	#BAND-START
	#	source CAA-EC 
	#	target CAY-MJ
	#	mig-rate-print 1.0
	#BAND-END

	#BAND-START
	#	source CAY-MJ
	#	target CAA-EC
	#	mig-rate-print 1.0
	#BAND-END

# Ancestral pops 2-3
	#BAND-START
	#	source CAA-EC-WC 
	#	target CAY-MJ
	#	mig-rate-print 1.0
	#BAND-END

	#BAND-START
	#	source CAY-MJ
	#	target CAA-EC-WC
	#	mig-rate-print 1.0
	#BAND-END

MIG-BANDS-END

