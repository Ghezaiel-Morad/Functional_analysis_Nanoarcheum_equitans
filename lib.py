from matplotlib_venn import venn3,venn3_circles
import pandas as pd
import matplotlib.pyplot as plt

def setup():
    print("### Hello wolrd, this program allows you to compare genome alignement results from the COG Database with a conventional one")
    print("-> Here, we will align Nanoarcheum equitans genome with sequences belonging to the NCBI NR databases")
    print("-> Blasts were runned using conventional parameters")
    print("Enjoy ! ")
    for i in range(4):
        print(" ")
    print(">>>> Setup comparaison : ")
    Coverage = ["c", "C", "Coverage"]
    Precision = ["p", "P", "Precision"]
    Thermococcales = ["T","Thermococcales"]
    Thermococcis = ["t","Thermococcis"]
    answer1 = input("Select COG method : (C = Coverage / P = precision ) : ")
    answer2 = input("Select NCBI organism specific database : (Thermococcales = T / Thermococcis = t) : ")
    COG_mode = 0
    NR_mode = 0
    if answer1 in Coverage:
        COG_mode = (1,3)
    if answer1 in Precision:
        COG_mode = (0,2)
    if answer1 not in Coverage:
        if answer1 not in Precision:
            print("\nInvalid COG mode entry, please retry\n")
            setup()
    if answer2 in Thermococcales:
        NR_mode = 5
    if answer2 in Thermococcis:
        NR_mode = 4
    if answer2 not in Thermococcales:
        if answer2 not in Thermococcis:
            print("\nInvalid NR mode entry, please retry\n")
            setup()
    print("\nResults:\n")
    return COG_mode, NR_mode

def get_file(file):
    path = file
    fo= open(path,"r")
    line = fo.readlines()
    # Close opend file
    fo.close()
    return line


def getCOGpval(line):
    pvals=[]
    for i in range(len(line)):
        Line = line[i].split("\t")
        for j in Line:
            if "e" in j and len(j)<9:
                p_val = (Line[0],j)
                pvals.append(p_val)

    return pvals

def getCOGfunction(line):
    function=[]
    for i in range(len(line)):
        Line = line[i].split("\t")
        func = (Line[0],Line[-1])
        function.append(func)
    return function

def getinfoNCBI(line):
    info = []
    for i in range(len(line)):
        Line = line[i].split(",")
        for j in Line:
            if "e" in j:
                inf = (Line[0],j)
                info.append(inf)
        if len(info)==0:
                inf = (Line[0],"NS")
                info.append(inf)

    return info


def get_best_hit(pval, function):
        best_hit_p_val = sorted([i[0] for i in pval if float(i[1])<1.0e-200])
        best_hit_function = [f for f in function if f[0] in best_hit_p_val]
        for i in best_hit_function:
            print(i)

