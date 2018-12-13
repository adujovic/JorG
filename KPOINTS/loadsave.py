import re
import numpy as np

def load_POSCAR(inputName):
    """
        Reading POSCAR file. Extensive testing required.
                                                        """

    # Returned data:
    directions    = []  # crystal directions

    with open(inputName,"r+") as inFile:
        for i,raw in enumerate(inFile.readlines()):
            line = re.sub("^\s*","",raw)  # remove all blank characters from begining of the line
            line = re.sub("\s+"," ",line) # replace all blank characters in a row to a single space
            if i == 0:
                comment = line[:-1]
            elif i == 1:
                scale = np.float64(line[:-1]) # scaling factor
                if(scale < 0.0):
                    scale = 1.0
            elif i in range(2,5):
                try:
                    directions.append(scale*np.fromstring(line,sep=" ")) # crystal directions
                except:
                    print("Error reading file %s in line %d:\nCan't convert crystal directions."%(inputName,i))
                    exit(-2)
                if(len(directions[-1]) != 3):   
                    print("Error reading file %s in line %d:\nCrystal directions has %d != 3 dimensions!"%(inputName,i,len(directions[-1])))
                    exit(-3)
            elif i == 5:        
                return directions

    return None
