
# Ouput options

OPT   +=  -DSTELLARAGE
#OPT   +=  -DMETALS
OPT   +=  -DOUTPUTPOTENTIAL
OPT   +=  -DOUTPUTACCELERATION
OPT   +=  -DOUTPUTCHANGEOFENTROPY
OPT   +=  -DOUTPUTTIMESTEP
OPT   +=  -DOUTPUTCOOLRATE
OPT   +=  -DMAGNETIC_FIELD
#OPT   +=  -DOUTPUT_DIV_CURL
OPT   +=  -DOUTPUT_VORTICITY

# simulation options

#OPT   +=  -DLONGIDS
OPT   +=  -DCOOLING	
OPT   +=  -DSFR   #only for runs with P-Gadget2
OPT   +=  -DWINDS # this enable the output delaytime 

#code options (gadget of gizmo)
OPT += -DCODE_GIZMO

#merger or isolated
#OPT += -DMERGER

# Reading/writting options
#OPT   +=  -DPRINT_HEADER
#OPT   +=  -DPRINT_ASCII_ONE_TYPE
#OPT   +=  -DPRINT_ASCII_ALL_TYPES
#OPT   +=  -DPRINT_ALL_ASCII   


#writting options Satellite/host after split
#OPT += -DPRINT_DATA

# Additional routines

OPT   +=  -DINPUT_OUTPUT_GADGET
OPT   +=  -DTRANSLATIONS_ROTATION
OPT   +=  -DANALYSIS
OPT   +=  -DBINS
OPT   += -DGADGET2_FORMAT

#option for surface density and dispersion velocity computation
#OPT += -DSURF_DENS_PROFILE

#option for softening potential function
#OPT += -DSOFT_POTENTIAL

#option for temperature/epicyclic frecuency analysis (needs SURF_DENS_PROFILE active)
#OPT += -DTOOMRE

#option for doing Star formation analysis
#OPT += -DSF_ANALYSIS

#option for mass profile (needs SURF_DENS_PROFILE active.. NO)
#OPT += -DMASS_PROFILE

#option for integrated sfr estimation at every snap
OPT += -DINTEGRATED_SFR

#option for gas sfr estimation at every snap
OPT += -DGAS_SFR


CC=gcc
UBICATIONGSL=/$HOME/local
#CFLAGS=-g -I. -I$(UBICATIONGSL)/include $(OPT)
CFLAGS= -g -I. -O0  -Wall -I$(UBICATIONGSL)/include $(OPT) 
LFLAGS= -lm -L$(UBICATIONGSL)/lib -lgsl -lgslcblas

clean:
	rm -rf *.o* *~ *.x

%.x:%.o
	$(CC) $^ $(LFLAGS) -o $@
#	cp $@ /home/lauraarboleda/local/bin
#cp $@ $@.x 



