time_division=6      #time split per hour
tolerance=0.1
#radius_GC = 4.25
scale = 1  #in the real size
radius_GC=1*scale
radius_MFR = radius_GC*0.5  #initial radius, it becomes larger when glomeruli  synapse with other cells
                            # no synapse observed by P15 at glomeruli

scale_radiation_range=1.0
scale_radiation_intensity=0.5
radiation1=radius_MFR*2.0*scale_radiation_range
radiation2=radius_MFR*5.0*scale_radiation_range
radiation3=radius_MFR*8.0*scale_radiation_range
rate1=1/2
rate2=1/8
rate3=1/15


intensity_radiation1=(0.01/time_division)*scale_radiation_intensity
intensity_radiation2=(0.005/time_division)*scale_radiation_intensity
intensity_radiation3=(0.001/time_division)*scale_radiation_intensity

height_PCL=radius_GC*42
area_length=height_PCL
area_width=height_PCL

Mig_Timing_Variation=10

simulation_time=24*20*time_division  #20days each 24 hours



