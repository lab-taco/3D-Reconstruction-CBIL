3D Reconstruction of the Cerebellar Input Layer<br/>
By Ikhwan Jeon, 2kaijeon@gmail.com<br/>
Supervisor: Dr. Taegon Kim<br/>

Poster: Ikhwan Jeon, Taegon Kim, ”Developmental modeling of modularized cerebellar input layer neural networks”, Poster presented at 25th Annual Meeting of the Korean Society for Brain and Neural Sciences. 2022 May 19; Incheon, Korea
[PDF](https://drive.google.com/file/d/1-6HaHTQEUoZtUSzVzadXosr22OKI98Md/view?usp=sharing)

Required Libraries:
Numpy, Matplotlib, VPython, Seaborn

Parameter Configuration:
-Physical Settings: Num_childs, Num_parents, Layer_expansion<br/>
-For debugging: Simulation_on, GC_implementation, vpython, time_sleep<br/>
-Data Save : Save_name_Synapses, Save_name_cell_locations<br/>

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

Analysis of extracted data
1. Connectivity
2. Cell Distribution
3. Modularity