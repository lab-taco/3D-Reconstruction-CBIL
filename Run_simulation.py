from Simulation.Parameters import *
from Simulation.Init_Setting.Set_SpaceTime import *
from Simulation.Init_Setting.Growth_Control import *

from Simulation.Cell_Gen import *
from Simulation.Cell_Body_Controls import *

from Simulation.Molecule_Controls import *
from Simulation.form_synapse import *

from Simulation.Data_Extraction import *
from vpython import *
import numpy as np
import math
import matplotlib.pyplot as plt 
import time



#from hanging_threads import start_monitoring
#start_monitoring(seconds_frozen=10, test_interval=100)

#collision_check_duration=6
Superposing=50
TWO_STEPS=True
if TWO_STEPS: step='2step'
else: step='1step'

def main(Num_childs, Num_parents, Simulation_on, Blink, \
            GC_implementation, Layer_expansion, FLAG_vpython, time_sleep, \
                Save_name_Synapses, Save_name_cell_locations, Save_name_Cell_Objects):

    start_time = time.time()
    
    print('-------Parameter status-----------------------------\n'\
        , 'Num_childs:', Num_childs, 'Num_parents:', Num_parents, 'Simulation:',Simulation_on\
        , 'Blink:',Blink,'\nGC_implementation:',GC_implementation, 'Layer_expansion:',Layer_expansion\
        , 'vpython:',FLAG_vpython, 'time_sleep:', time_sleep\
        , 'Superposing:',Superposing, 'TWO_STEPS:',TWO_STEPS)
    print('DATA_PATH:', DATA_PATH, '\n-------------------------------------------------------')
    #Save_name_Synapses, Save_name_cell_locations='save_tmp', 'save_tmp'
    Simulation_on=True
    #GC_implementation=False
    #vpython=Falses
    time_sleep=False
    Analysis_on=not Simulation_on

    #Calculate number of cells
    max_MF=Num_parents*Num_childs + Num_parents
    MF_total=max_MF
    Max_GC=3*max_MF
    
    """
    Layer division, time display, Sample MFs ----------------------------------------------------
    """     
    if Simulation_on:
        time_display=class_time_display(area_length, height_PCL, time_division=time_division, vpython=FLAG_vpython)        
        layer_division(area_length, height_PCL, area_width, labeling=False, vpython=FLAG_vpython) #against PCL, the ceiling
        if Layer_expansion: 
            IGL_expansion=moving_layer(area_length, height_PCL, area_width, vpython=FLAG_vpython) #against WM, the bottom
        print('Simulation starts...')
        MFs, Color_map_MF=sample_MFs_clusters(Num_parents, Num_childs, max_MF, MF_Mig_Timing_Variation, \
                                area_length,area_width, height_PCL, vpython=FLAG_vpython)
        
        #num_cells_per_MF_clusters(MFs)
    
    
    #while True: rate(30)
    """
    Growth controls of Granule cell & IGL depth.
    Turn on graph to be True if you want to see the # of GCs.
    """    
    depth_IGL_table=exponential_growth(graph=Analysis_on)    
    #num_GCs=logistic_growth('GC', Max_GC,graph=Analysis_on)
    num_GCs=linear_growth(Max_GC,graph=Analysis_on)
    if Analysis_on: sys.exit()
    
    """
    MF Blinking test
    if Blink: Blink(MFs)
    
    Intensity mapping?
    show_intensity=False
    if show_intensity:
        intensity_map=intensity_mapping(MFs)
        #data_save('intensity_map', intensity_map)
        print('int_map shape:',intensity_map.shape)
    """
    
    """
    Main loop ----------------------------------------------------
    """
    GCs=[]    
    #Num_variation_GC_color = 100 if Max_GC>100 else Max_GC    
    Num_variation_GC_color = GC_Mig_Timing_Variation
    Color_map_GC = colormap(Num_variation_GC_color, 'GC')
    current_num_GCs=0  #current num cell
    Simul_Start, Postnatal_days, P7_passed, P14_passed, P20_passed=False, False, False, False, False
    FLAG_SYNAPSTIC_FORMATION=True 
    #End_Time=21*24*time_division
    #End_Time=20*24*time_division
    
    MF_activity_pattern=curvs_generations(MF_Mig_Timing_Variation,time_division, SP=Superposing,draw=False)
    if len(Color_map_MF)!=len(MF_activity_pattern): raise Exception('coloring & mf activity curve mismatch')
    while(Simulation_on and not P20_passed): #main loop
        #gives time final GCs to finish migration and form syanpse

        if time_sleep:
            sleep(0.00001)
            
        #MF grow up toward PCL till P7 at which 30% of MFs invade PCL & GCs start migrate down to IGL
        if not P7_passed: 
            MF_growup(MFs)
        elif not P14_passed: 
            MFs_repulsed(MFs) 

        if not Postnatal_days: # Embryonic stage
            if not Simul_Start: 
                print("Embryonic stage")
                Simul_Start=True
            num_MFs_above_PCL=0
            for mf in MFs:
                if mf.body.pos.y>=height_PCL:
                    num_MFs_above_PCL+=1
            #Start of postnatal days P0 at which 16% of MFs invade PCL
            if num_MFs_above_PCL>=int(0.16*MF_total):                 
                Postnatal_days=True
                time_display.display.text=('P%d'%time_display.day)
        
        elif Postnatal_days: # Postnatal stage            
            #Postnatal Timing
            if not P7_passed and time_display.counter==7*24*time_division: P7_passed=True # 7days * 24 hours
            if not P14_passed and time_display.counter==14*24*time_division: # 14days * 24 hours
                P14_passed=True   #at the end of P14, take all MFs into collision control process
                MFs_migration_complete(MFs)
                all_cells=GCs + MFs
                all_superposition_check(all_cells)
                                    

            #IGL expansion
            #if time_display.counter<len(depth_IGL_table):
            
            #current_depth=depth_IGL_table[time_display.counter-1]
            current_depth=depth_IGL_table[time_display.counter]
            if Layer_expansion: #expand IGL according to current depth by table
                IGL_expansion.expand(current_depth)
                #if P20_passed: IGL_expansion.expansion_end()
            #Calculation of MF Molecular Activity
            #MF_activities=MF_activity_coordinator(time_display.counter, time_division)
                    

            #Generation of GCs
            if GC_implementation and not P20_passed:
                current_num_GCs=int(num_GCs[(time_display.counter)])
                while len(GCs)<current_num_GCs: 
                    cell_position=random_sample_2d(area_length,area_width,height_PCL,'GC')
                    Ind_GC_Color = int(len(GCs)/Max_GC*Num_variation_GC_color)
                    cell_color=Color_map_GC[Ind_GC_Color]
                    GCs.append(Cells('GC', init_position=cell_position, color_=cell_color, vpython=FLAG_vpython))
            # Migration & Synapse formation
            if not P14_passed: target_cells=GCs
            else: target_cells=GCs + MFs  #keep update target_cells for collision control of MFs
            #GC_migrates(target_cells, current_depth, time_display.counter, collision_check_duration)
            Cell_migratedown(target_cells, current_depth, time_display.counter)
            #Molecule_absorption(GCs, MFs, collision_check_cells=target_cells, vpython=FLAG_vpython)
            Stop_or_not(GCs, target_cells, FLAG_vpython)
            if current_num_GCs>0:
                if GCs[0].flag_arrived_final_destination:
                    Sum_mf_activities=MF_activities_summation(MFs, Color_map_MF, MF_activity_pattern, time_display)
                    if Sum_mf_activities<0: 
                        print('Sum_mf_activities < 0 at time:', time_display.counter)                    
                    else:
                        Form_Synapse(GCs, MFs, Sum_mf_activities, Two_steps=TWO_STEPS)
                        #Initial_Contacts(GCs, MFs)

            time_display.time_count() # count time steps for postnatal days
            if time_display.counter==20*24*time_division: P20_passed=True # 20days * 24 hours 
    print('Simulation end')
    
    #while True:
    #    sleep(10)
    #    pass
    #sys.exit()

    #Simulation ends
    final_collision_check(GCs+MFs)

    elapsed_time = time.time() - start_time
    print('After processes end, total elapsed time:', elapsed_time)
    
    
    if not Save_name_Cell_Objects=='No save':
        from datetime import date, datetime    
        save_time=datetime.now().strftime("%d-%m-%Y-%H%M")        
        data_label='SP{}_{}_{}{}'.format(str(Superposing), step, str(Num_parents),'Parents')

        #data_save(data_name, data, PATH, dir_name)
        data_save('GC_Objects'+data_label, GCs, DATA_PATH, save_time)
        data_save('MF_Objects'+data_label, MFs, DATA_PATH, save_time)

        data_save('GC_Colormap'+data_label, Color_map_GC, DATA_PATH, save_time)
        data_save('MF_Colormap'+data_label, Color_map_MF, DATA_PATH, save_time)
    
        print('Cell Object data saved')

    print('Len GC:', len(GCs), 'Len GC Color', len(Color_map_GC))
    print('Len MF:', len(MFs), 'Len MF Color', len(Color_map_MF))

    num_synapse(GCs, MFs)
    len_synapse(GCs, MFs)
    #sys.exit()  #Ver. 220514
    #early=[]
    #late=[]
    #for cell in GCs:
    #    if cell.body.color==color.red:
    #        early.append(cell)
    #    elif cell.body.color==color.green:
    #        late.append(cell)
#
    #print('# GCs of  early:',len(early),\
    #             'which is %.2f percent of all cells'%(len(early)/len(all_cells)*100))
    #print('# all_cells:',len(all_cells))
    #after_migration_coloring(early, all_cells)
#
#
    #if not GC_implementation or not Simulation_on:
    #    print('GC_implementation off, simulation off')
    #    sys.exit()

    #Save_name_Synapses='No save'
    #Save_name_Synapses='tolerance'
    Save_name_Synapses='Graph_2023May'
    edges_for_all_MFs=[]
    if not Save_name_Synapses=='No save':
        for index, mf in enumerate(MFs):
            synapses=mf.synapse_partners
            edges=[] #element structure: [index_mf, index_parter_gcs]
            for syn_partner in synapses:
                edges.append([index, syn_partner])
                #edges.append([index, syn_partner, [GCs[syn_partner].color] ])
            #edges.append(mf.color)
            edges_for_all_MFs.append(edges)
    
        #edges_for_all_MFs.append(['Num_GCs', Max_GC])
        #edges_for_all_MFs= edges_for_all_MFs.astype('object')
        edges_for_all_MFs= np.asarray(edges_for_all_MFs, dtype='object')
        #data_save('Volume_filling2', edges_for_all_MFs)

        from datetime import date, datetime    
        now=datetime.now().strftime("%d-%m-%Y-%H%M")
        data_name=Save_name_Synapses+'_{}_{}_{}'.format(str(Num_parents),'MFs', 'edges')+now
        data_save(data_name, edges_for_all_MFs, DATA_PATH, save_time)
        print('Synapse graph data saved')

    print('Simulation Data Save finished')
    """
    Have to change w.r.t new Color_map_GC type
    early=[]
    mid=[]
    late=[]
    for cell in GCs:
        if cell.body.color==color.red:
            early.append(cell)
        elif cell.body.color==color.green:
            late.append(cell)
        else: mid.append(cell)

    print('# GCs of  early:',len(early),\
                 'which is %.2f percent of all GCs'%(len(early)/len(GCs)*100))
    print('# all_cells:',len(all_cells))
    """
    #after_migration_coloring(early, all_cells)
    
    #sleep(100000)
    while True: rate(30)
    sys.exit()
    if not GC_implementation or not Simulation_on:
        print('GC_implementation off, simulation off')
        sys.exit()
    if not Save_name_cell_locations=='No save':
        Save_name_cell_locations='Tomography_linears2'
        #data_name=Save_name_cell_locations+'_{}_{}_{}'.format(str(Num_parents) ,'MFs' ,'normalGCs')
        #data_extraction_3d([c for c in GCs if c not in early], data_name)
        data_name=Save_name_cell_locations+'_{}_{}_{}'.format(str(Num_parents) ,'MFs' ,'earlyGCs')+now
        data_extraction_3d(early, data_name)
        data_name=Save_name_cell_locations+'_{}_{}_{}'.format(str(Num_parents) ,'MFs' ,'midGCs')+now
        data_extraction_3d([c for c in mid], data_name)
        data_name=Save_name_cell_locations+'_{}_{}_{}'.format(str(Num_parents) ,'MFs' ,'lateGCs')+now
        data_extraction_3d([c for c in late], data_name)                                    
        data_name=Save_name_cell_locations+'_{}_{}_{}'.format(str(Num_parents) ,'MFs' ,'MFs')+now
        data_extraction_3d([c for c in MFs], data_name)

        print('Data extraction finished')

    

import argparse
if __name__ == '__main__':

    #start_time = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--Num_childs', type=int, dest='Num_childs',
                        default=DEFEAULT_Num_childs,
                        help='Number of child cells at each Mossy fiber cluster.')
    parser.add_argument('--Num_parents', type=int, dest='Num_parents',
                        default=DEFEAULT_Num_parents,
                        help='Number of parents cells(# of clusters) for Mossy fiber clusters.')

    parser.add_argument('--Simulation_on', type=int, dest='Simulation_on',
                        #boolen for argparse is tricky, So used 0,1 instead
                        default=1,
                        help='if simulation off, migration does not start.\
                              if simulation on with vpython off, migration start but visualization does not')
    parser.add_argument('--Blink', type=int, dest='Blink',
                        default=0,
                        help='Blinking represents the activation pattern of Mossy fibers.')

    parser.add_argument('--GC_implementation', type=int, dest='GC_implementation',
                        default=1,
                        help='if GC_implementation is off, only MFs show up.')
    parser.add_argument('--Layer_expansion', type=int, dest='Layer_expansion',
                        default=1,
                        help='Visualization of Layers and layer_expansion.')

    parser.add_argument('--FLAG_vpython', type=int, dest='FLAG_vpython',
                        default=1,
                        help='if vpython off, only calculation start without visualization by vpython.')

    parser.add_argument('--time_sleep', type=int, dest='time_sleep',
                        default=0,
                        help='To see simulations slowly, insert time delay in the loop.')

    parser.add_argument('--Save_name_Synapses', type=str, dest='Save_name_Synapses',
                        default='save_tmp',
                        help='Extract and Save the connection of Synapses made through a simulation.\
                                input "No save" for saving nothing')

    parser.add_argument('--Save_name_cell_locations', type=str, dest='Save_name_cell_locations',
                        default='save_tmp',
                        help='Save_name_cell_locations.\
                             input "No save" for saving nothing')

    parser.add_argument('--Save_name_Cells', type=str, dest='Save_name_Cell_Objects',
                        default='save_cell_objects',
                        help='Save_name_cell_locations.\
                             input "No save" for saving nothing')    
    
    args = parser.parse_args()
    
    main(args.Num_childs, args.Num_parents, args.Simulation_on, args.Blink, \
         args.GC_implementation, args.Layer_expansion, args.FLAG_vpython, args.time_sleep, \
         args.Save_name_Synapses, args.Save_name_cell_locations, args.Save_name_Cell_Objects)
    


    