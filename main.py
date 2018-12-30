from lib import get_file,setup,getCOGpval, getCOGfunction, getinfoNCBI, get_best_hit
from matplotlib_venn import venn3,venn3_circles
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


files = ["HMM_precision.txt","HMM_coverage.txt","DIAMOND_precision.txt","DIAMOND_coverage.txt","Blastp_thermococcis_precision.txt","Blastp_thermococcales_precision.txt"]
filename = ["HMM one-to-one","HMM coverage","DIAMOND one-to-one","DIAMOND coverage","Blastp NR Thermococcis","Blastp NR Thermococcales"]


if __name__ == '__main__':
    while True:

        # We store the user answers
        answer1,answer2 = setup()

        # Get functional annotations from HMM method
        HMM_function = getCOGfunction(get_file(files[answer1[0]]))

        # Same for DIAMOND
        DIAMOND_function = getCOGfunction(get_file(files[answer1[1]]))

        # Get p_values from HMM method
        HMM_p_vals = getCOGpval(get_file(files[answer1[0]]))

        # Same for DIAMOND
        DIAMOND_p_vals = getCOGpval(get_file(files[answer1[1]]))

        # Finally, get accession numbers from NR
        NR_accession = getinfoNCBI(get_file(files[answer2]))


        # P-values and functional annotations are pulled with their respective accession number
        HMMset = [i[0] for i in HMM_function]
        DIAMONDset =[i[0] for i in DIAMOND_function]
        NRset = [i[0] for i in NR_accession]

        # We set our data in order to remove duplicates
        set1 = set(HMMset)
        set2 = set(DIAMONDset)
        set3 = set(NRset)


        # We pull intersecting values between HMM and DIAMOND methods
        HMM_inter_NR = set1.intersection(set3)
        DIAMOND_inter_NR = set2.intersection(set3)
        inter_accession = set(HMM_inter_NR).intersection(set(DIAMOND_inter_NR))

        # Then we reorganize our data in order to create a dataframe
        PDIA = sorted({Id: pval for Id, pval in DIAMOND_p_vals if Id in inter_accession})
        ids = [i for i in PDIA]
        HMM_p = [pval for id, pval in sorted(HMM_p_vals) if id in ids]
        HMM_f = [function for id, function in sorted(HMM_function) if id in ids]
        DIAMOND_p = [pval for id, pval in sorted(DIAMOND_p_vals) if id in ids]
        DIAMOND_f = [function for id, function in sorted(DIAMOND_function) if id in ids]



        # Some values are lacking, we fill the list with 0
        while len(HMM_p)<207:
            HMM_p.append(0)

        # We create a dataframe to store functional annotations and p-values
        d = {'Accession': ids, 'HMM function': HMM_f, 'HMM p values': HMM_p,'DIAMOND function': DIAMOND_f, 'DIAMOND p values': DIAMOND_p}
        df = pd.DataFrame(data=d)
        print(df)

        # Then, we save the dataframe as a csv file
        table_name = "Functional_annotation_"+filename[answer1[0]]+"_"+filename[answer1[0]]+"_"+filename[answer2]+".csv"
        df.to_csv(table_name,sep=',', header=True,columns = ['Accession', 'HMM function', 'HMM p values','DIAMOND function','DIAMOND p values'])
        print("Number of orthologs found in NR Database : ", len(set3))
        print("Number of orthologs found in COG Database using HMM method : ", len(set1))
        print("Number of orthologs found in COG Database using DIAMOND method : ", len(set2))
        print("Number of intersecting orthologs : ",len(inter_accession))


        # Finally, we create a Venn diagram to get a visualization of intersecting orthologs
        plt.title(filename[answer1[0]] + " vs " + filename[answer1[1]] + " vs " + filename[answer2])
        ven = venn3((set1, set2, set3), (filename[answer1[0]], filename[answer1[1]], filename[answer2]))
        venn = venn3_circles((set1,set2,set3))
        plt.show()
        get_best_hit(HMM_p_vals, HMM_function)
        print("#################################################")

















