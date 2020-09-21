import math
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from skbio import DistanceMatrix
from skbio.stats.distance import mantel

###### import original evn map


MAP_FILE=open('test_SDM_vector_small_transition_nrow_40_ncol_60.txt','r')
ORIGINAL_MAP=[]
for line in MAP_FILE:
    line=line.strip()
    if line=="NA":
        line=0.0
    else:
        line=float(line)
        if line<=0:
            line=0
    ORIGINAL_MAP.append(line)

ORIGINAL_MAP=np.array(ORIGINAL_MAP)
ORIGINAL_MAP=np.transpose(ORIGINAL_MAP.reshape(40,60))



###### import generated occurance map

for Kay in [1,2,3,5,8]:
    for SIMGA in [2]:

            
        OC_MAP_FILE=open('Occurance_Matrix_{}_{}'.format(Kay,SIMGA),'r')



        OC_MAP=[]
        for line in OC_MAP_FILE:
            line=line.strip()
            if line=="NA":
                line=0.0
            else:
                line=float(line)
                if line<=0:
                    line=0
                
            OC_MAP.append(line)

        OC_MAP=np.array(OC_MAP)
        OC_MAP=OC_MAP.reshape(60,40)
        OC_MAP_normal=(OC_MAP-np.amin(OC_MAP))/(np.amax(OC_MAP)-np.amin(OC_MAP))
        OC_MAP_avg=(OC_MAP/150)/OC_MAP.size
        print(OC_MAP_avg)


        COMP_MAP_FILE=open('Competition_Matrix_{}_{}'.format(Kay,SIMGA),'r')

        COMP_MAP=[]
        for line in COMP_MAP_FILE:
            line=line.strip()
            if line=="NA":
                line=0.0
            else:
                line=float(line)
                if line<=0:
                    line=0
                
            COMP_MAP.append(line)
        
        COMP_MAP=np.array(COMP_MAP)
        COMP_MAP=COMP_MAP.reshape(60,40)
        
        FITN_MAP_FILE=open('Fitness_Matrix_{}_{}'.format(Kay,SIMGA),'r')

        FITN_MAP=[]
        for line in FITN_MAP_FILE:
            line=line.strip()
            if line=="NA":
                line=0.0
            else:
                line=float(line)
                if line<=0:
                    line=0
                
            FITN_MAP.append(line)
        
        FITN_MAP=np.array(FITN_MAP)
        FITN_MAP=FITN_MAP.reshape(60,40)

        #######
        ## Calculating differences 


        OC_MAP_normal_list=[ x for y in list(OC_MAP_normal) for x in y]
        OC_MAP_avg_list=[ x for y in list(OC_MAP_avg) for x in y]
        ORIGINAL_MAP_list=[ x for y in list(ORIGINAL_MAP) for x in y]
        DIFF_MAP=[]
        DIFF_MAP2=[]
        for k in range(0,len(ORIGINAL_MAP_list)):
            DIFF_MAP.append(OC_MAP_normal_list[k]-ORIGINAL_MAP_list[k])
            DIFF_MAP2.append(OC_MAP_avg_list[k]-ORIGINAL_MAP_list[k])



        DIFF_MAP=(np.asarray(DIFF_MAP)).reshape(60,40)
        DIFF_MAP2=(np.asarray(DIFF_MAP2)).reshape(60,40)
        CORREL=np.corrcoef(ORIGINAL_MAP.reshape(1,ORIGINAL_MAP.size),OC_MAP.reshape(1,OC_MAP.size),rowvar=True)






        ##############
        ##### Plotting 


        normal = matplotlib.colors.Normalize(np.amin(DIFF_MAP),np.amax(DIFF_MAP))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.title('Difference between Occurance Map - SD Map',size=50)
        plt.text(-7,-5,'Mean differense between occurance matrix (normalised) and SD Map : {}'.format(np.mean(abs(DIFF_MAP))),fontsize=65)
        plt.text(-7,-7,'Correlation between Occurance map and Original map : {}'.format(CORREL[0,1]),fontsize=65)
        plt.pcolormesh(DIFF_MAP,cmap=colormap,norm=normal)
        plt.colorbar().ax.tick_params(labelsize=90) 

        plt.savefig('Differences_plot_Normalised_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()





        normal = matplotlib.colors.Normalize(np.amin(OC_MAP_normal),np.amax(OC_MAP_normal))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.title('Normalised Occur',size=50)
        plt.pcolormesh(OC_MAP_normal,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 

        plt.savefig('Normalised_Occur_plot_Normalised_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()




        normal = matplotlib.colors.Normalize(np.amin(OC_MAP),np.amax(OC_MAP))
        colormap = plt.get_cmap("plasma")
        plt.figure(figsize=(60, 40))
        plt.title('Difference between Occurance Map - SD Map',size=50)
        plt.pcolormesh(OC_MAP,cmap=colormap)
        plt.colorbar().ax.tick_params(labelsize=90) 

        plt.savefig('Occurance_Map_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()
        
        
        
        
        ##### Scatter Plots

        
        ORIGINAL_MAP_list=[ x for y in list(ORIGINAL_MAP) for x in y]
        OC_MAP_list=[ x for y in list(OC_MAP) for x in y]
        FITN_MAP_list=[ x for y in list(FITN_MAP) for x in y]
        COMP_MAP_list=[ x for y in list(COMP_MAP) for x in y]
        DIFF_MAP_list=[ x for y in list(DIFF_MAP) for x in y]
        
        
        
        
        plt.figure(figsize=(60, 40))
        plt.title('Correlation of Differences map and Original Map',size=50)
        plt.xlabel("Difference per Grid",size=50)
        plt.ylabel("SDM Value per Grid",size=50)
        plt.tick_params(axis='both', which='major', labelsize=25)
        plt.tick_params(axis='both', which='minor', labelsize=16)
        plt.scatter(DIFF_MAP_list, ORIGINAL_MAP_list, color='r')
        plt.savefig('Correlation_Between_Differences_OriginalMap_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()

        plt.figure(figsize=(60, 40))
        plt.title('Correlation of Differences map and Competition',size=50)
        plt.xlabel("Difference per Grid",size=50)
        plt.ylabel("Competition per Grid",size=50)
        plt.tick_params(axis='both', which='major', labelsize=25)
        plt.tick_params(axis='both', which='minor', labelsize=16)
        plt.scatter(DIFF_MAP_list, COMP_MAP_list, color='r')
        plt.savefig('Correlation_Between_Differences_Competition_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()
        
        plt.figure(figsize=(60, 40))
        plt.title('Correlation of Differences map and Fitness',size=50)
        plt.xlabel("Difference per Grid",size=50)
        plt.ylabel("Fitness per Grid",size=50)
        plt.tick_params(axis='both', which='major', labelsize=25)
        plt.tick_params(axis='both', which='minor', labelsize=16)
        plt.scatter(DIFF_MAP_list, FITN_MAP_list, color='r')
        plt.savefig('Correlation_Between_Differences_Fitnes_K{}_Sig{}.pdf'.format(str(Kay),str(SIMGA)))
        plt.close()





        ####Printing out stats

        print('Actual average population per generation : ',np.sum(OC_MAP)/100)
        print('Expected population : ',np.sum(ORIGINAL_MAP)*Kay)




        

        print(CORREL)
        print('Mean differense between occurance matrix (normalised) and SD Map : ', np.mean(abs(DIFF_MAP)))
        print('Mean differense between occurance matrix (normalised) and SD Map : ', np.mean(abs(DIFF_MAP2)))
