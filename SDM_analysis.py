import math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


for SIMGA in [2]:
    for Kay in [5]:
        FILENAME='OUT_SED_1000_SDM_Sigma_{}_K_{}.all_ind_info'.format(SIMGA,Kay)


        DATA_FILE=open(FILENAME,'r')

        DATA_DC={}
        LABELS=DATA_FILE.readline().strip().split()

        for line in DATA_FILE:
            line=line.strip().split()
            for numb in range(0,len(line)):
                line[numb]=float(line[numb])
            
            if line[0] not in DATA_DC:
                DATA_DC[line[0]]=[]
            else:
                DATA_DC[line[0]].append(line)

        print(len(DATA_DC))


        ###################### plot for population through time
        GENS=[]
        POP_SIZE=[]

        for k in range(2,len(DATA_DC)+1):
            GENS.append(k)
            POP_SIZE.append(len(DATA_DC[k]))


        plt.figure(figsize=(9, 3))
        plt.plot(GENS,POP_SIZE)

        plt.savefig('Popsize_plot_{}_{}.pdf'.format(str(Kay),str(SIMGA)))



        ########import and plot original env map

        MAP_FILE=open('test_SDM_vector_small_transition_nrow_40_ncol_60.txt','r')
        ORIGINAL_MAP=[]
        for line in MAP_FILE:
            line=line.strip()
            if line=="NA":
                line=0.0
            if line!='NA' and float(line)<=0:
                line=0.0
            else:
                line=float(line)
            ORIGINAL_MAP.append(line)

        ORIGINAL_MAP=np.array(ORIGINAL_MAP)
        ORIGINAL_MAP=np.transpose(ORIGINAL_MAP.reshape(40,60))


        normal = matplotlib.colors.Normalize(np.amin(ORIGINAL_MAP),np.amax(ORIGINAL_MAP))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.pcolormesh(ORIGINAL_MAP,cmap=colormap,norm=normal)
        plt.colorbar().ax.tick_params(labelsize=90) 
        plt.savefig('Enviromental_Map_plot.pdf')
        plt.close()

        ########Calculate and plot the occurance of individuals per square





        BOXES=[[] for x in range(0,40*60)]
        BOXES_FIT=[[] for x in range(0,40*60)]
        BOXES_COMP=[[] for x in range(0,40*60)]
        BOXES_OPTIM=[[] for x in range(0,40*60)]


        for gen in range(50,150):
            print(gen)
            BoxCounter=0
            for Xbox in range(0,40):
                for Ybox in range(0,60):
                    for ind in DATA_DC[gen]:
                        if ( (ind[2]>=Xbox) & (ind[2]<=Xbox+1) & (ind[3]>=Ybox) &  (ind[3]<Ybox+1) ):
                            BOXES[BoxCounter].append(1.0)
                            BOXES_FIT[BoxCounter].append(ind[4])
                            BOXES_COMP[BoxCounter].append(ind[7])
                            BOXES_OPTIM[BoxCounter].append(ind[6]) 
                            
                    BoxCounter+=1

        BOXES=[sum(x)/100 for x in BOXES]
        BOXES_FIT=[np.mean(x) for x in BOXES_FIT]
        BOXES_COMP=[np.mean(x) for x in BOXES_COMP]
        BOXES_OPTIM=[np.mean(x) for x in BOXES_OPTIM]


        BOXES=np.array(BOXES)
        BOXES=np.transpose(BOXES.reshape(40,60))
        print(BOXES)

        BOXES_FIT=np.transpose(np.array(BOXES_FIT).reshape(40,60))
        BOXES_COMP=np.transpose(np.array(BOXES_COMP).reshape(40,60))
        BOXES_OPTIM=np.transpose(np.array(BOXES_OPTIM).reshape(40,60))
        

        print(BOXES_FIT)
        print(BOXES_COMP)




        normal = matplotlib.colors.Normalize(np.amin(BOXES),np.amax(BOXES))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.pcolormesh(BOXES,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 
        plt.savefig('Distribution_plot_{}_{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()

        normal = matplotlib.colors.Normalize(np.amin(BOXES_FIT),np.amax(BOXES_FIT))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.pcolormesh(BOXES_FIT,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 
        plt.savefig('Fitness_Distribution_plot_{}_{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()

        normal = matplotlib.colors.Normalize(np.amin(BOXES_COMP),np.amax(BOXES_COMP))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.pcolormesh(BOXES_COMP,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 
        plt.savefig('Competition_Distribution_plot_{}_{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()

        normal = matplotlib.colors.Normalize(np.amin(BOXES_OPTIM),np.amax(BOXES_OPTIM))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.pcolormesh(BOXES_OPTIM,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 
        plt.savefig('Optimum_Distribution_plot_{}_{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()




        ##################### output matrices
        OUT_MATRIX=open('Occurance_Matrix_{}_{}'.format(str(Kay),str(SIMGA)),'w')
        OUT_FIT_MAT=open('Fitness_Matrix_{}_{}'.format(str(Kay),str(SIMGA)),'w')
        OUT_COMP_MAT=open('Competition_Matrix_{}_{}'.format(str(Kay),str(SIMGA)),'w')
        OUT_OPT_MAT=open('Optimum_Matrix_{}_{}'.format(str(Kay),str(SIMGA)),'w')


        BOXES=list(BOXES)
        BOXES= [x for sublist in BOXES for x in sublist ] 

        BOXES_FIT=list(BOXES_FIT)
        BOXES_FIT= [x for sublist in BOXES_FIT for x in sublist ] 

        BOXES_COMP=list(BOXES_COMP)
        BOXES_COMP= [x for sublist in BOXES_COMP for x in sublist ] 

        BOXES_OPTIM=list(BOXES_OPTIM)
        BOXES_OPTIM= [x for sublist in BOXES_OPTIM for x in sublist ] 



        for k in range(0,len(BOXES)):
            OUT_MATRIX.write(str(BOXES[k])+'\n')
            OUT_FIT_MAT.write(str(BOXES_FIT[k])+'\n')
            OUT_COMP_MAT.write(str(BOXES_COMP[k])+'\n') 
            OUT_OPT_MAT.write(str(BOXES_OPTIM[k])+'\n') 
            
