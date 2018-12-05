class color:
    """
        From:
        https://stackoverflow.com/questions/8924173/how-do-i-print-bold-text-in-python
    """    
    BF = '\033[1m'
    IT = '\033[3m'
    UN = '\033[4m'
    BLINK = '\033[5m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    END = '\033[0m'

def print_atom(atom,end='\n',vector=color.END,element=color.BF+color.YELLOW):
    print(element+str(atom[0])+color.END+' [ '+vector+'% 2.5f % 2.5f % 2.5f'%(*atom[1],)+color.END+' ]',end=end)
