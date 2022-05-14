3D Reconstruction of the Cerebellar Input Layer<br/>
By Ikhwan Jeon, 2kaijeon@gmail.com<br/>
Supervisor: Dr. Taegon Kim<br/>

Required Libraries:
Numpy, Matplotlibe, VPython, Seaborn

Parameter Configuration:
-Physical Settings: Num_childs, Num_parents, Layer_expansion<br/>
-For debugging: Simulation_on, GC_implementation, vpython, time_sleep<br/>
-Data Save : Save_name_Synapses, Save_name_cell_locations<br/>
-To be removed: Blink, show_intensity

Description:<br/>
-Simulation<br/>
1. Layer Expansion<br/>
2. MFs<br/>
3. GCs -migration stop either if goes below current IGL or probablilistic choice from molecule gauge<br/>
4. Molecule absorption and migration end<br/>
5. Synapse Formation<br/>
-Data analysis<br/>
1. Molecular activity - Connectivity Relationship<br/>
2. Spatial Distribution of Cell body

Simulation Timeline<br/>
1. Configuration<br/>
2. Loop<br/>
For all time steps:<br/>
- ~P7: MFs move around<br/>
- P7~: <br/>
 a. MF_activity_coordination<br/>
 b. Layer_expansion<br/>
 c. Generation of GCs<br/>
 d. Migration<br/>
 e. Synapse formation<br/>
3. Data Extraction
 - Cell objects: GCs & MFs

To be added:<br/>
- Continuous variation of molalcular activity<br/>
- Equalizers of MF Activity Variation<br/>
- Visualization of Synapse formation?<br/>

Analysis of extracted data
1. Connectivity
2. Cell Distribution
3. Modularity