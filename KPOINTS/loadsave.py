import re
import numpy as np

def load_POSCAR(inputName):
    """
        Reading POSCAR file. Extensive testing required.
                                                        """
    # Returned data:
    directions    = []  # crystal directions

    with open(inputName,"r+") as inFile:
        rawText = inFile.readlines()[:5]
        for i in range(5):
            rawText[i] = re.sub("^\s+","",rawText[i])  # remove all blank characters from the begining of the rawText[i]
            rawText[i] = re.sub("\s+"," ",rawText[i])  # replace all blank characters in a row to a single space
            rawText[i] = re.sub("\s+$","",rawText[i])  # remove all blank characters from the end of the rawText[i]

        scale   = np.float64(rawText[1])
        for i in range(2,5):
            directions.append(scale*np.fromstring(rawText[i],sep=" ")) # crystal directions
        return directions

    return None
