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
    output  = element + str(atom[0]) + color.END
    output += ' [ ' + vector
    output += '{:= 10.5f} {:= 10.5f} {:= 10.5f}'.format(*atom[1],)
    output += color.END + ' ] '
    print(output,end=end)


def print_case(caseID,atom,atomID,wyckoffPosition,distance,end='\n',caseStyle=color.IT+color.CYAN,numberStyle=color.IT,distanceStyle=color.END):
    output  = "Case "
    output += caseStyle+"{:3d}".format(caseID)+color.END
    output += " atom No. "
    output += numberStyle+"{:3d}".format(atomID)+color.END
    output += distanceStyle+" @ {:1s} & {:6.2f} Å | ".format(wyckoffPosition,distance)+color.END
    print(output,end='')
    print_atom(atom,end=end)

import numpy as np
def print_crystal(directions,cell,linewidth=42,atomNames=None):
    if(linewidth < 42):
        print("Too low linewidth!")
    
    report = {}
    if atomNames is None:
        for atom in cell:
            if atom[0] in report.keys():
                report[atom[0]] += 1
            else:
                report[atom[0]] = 1
    else:
        for atom in cell:
            if atomNames[atom[0]] in report.keys():
                report[atomNames[atom[0]]] += 1
            else:
                report[atomNames[atom[0]]] = 1

    print("+"+linewidth*'-'+"+")
    print("|"+color.BF+"Composition:".center(linewidth)+color.END+"|")
    sumOfAtoms = 0
    for record in report:
        name   = "{:4s}:   ".format(record)
        number = "#{:4d}".format(report[record])
        label  = color.BF+color.YELLOW+name+color.END+color.IT+number+color.END
        offset = len(label)-len(name)-len(number)
        print('|'+label.center(linewidth+offset)+'|')
        sumOfAtoms += report[record]

    name   = "cell:   "
    number = "#{:4d}".format(sumOfAtoms)
    label  = color.BF+color.DARKCYAN+name+color.END+color.IT+number+color.END
    offset = len(label)-len(name)-len(number)
    print('|'+label.center(linewidth+offset)+'|')
    

    print("|"+linewidth*' '+"|")
    print("|"+color.BF+'Crystal axes:'.center(linewidth)+color.END+"|")

    names = [color.BF+color.RED  +'a⃗'+color.END,
             color.BF+color.GREEN+'b⃗'+color.END,
             color.BF+color.BLUE +'c⃗'+color.END]
    for n,d in zip(names,directions):
        data = "{:s} = [ {:= 10.5f} {:= 10.5f} {:= 10.5f} ]".format(n,*d,)
        print('|'+data.center(linewidth+len(n)-1)+"|")

    print("|"+linewidth*' '+"|")
    print("|"+color.BF+'Crystal directions:'.center(linewidth)+color.END+"|")

    names = [color.RED+'a'+color.END, color.GREEN+'b'+color.END, color.BLUE+'c'+color.END]
    for n,d in zip(names,directions):
        data = "{:s} = {:= 8.5f} Å".format(n,np.linalg.norm(d))
        print('|'+data.center(linewidth+len(n)-1)+"|")

    print("|"+linewidth*' '+"|")
    print("|"+color.BF+'Crystal angles:'.center(linewidth)+color.END+"|")

    names = [color.RED+'α'+color.END, color.GREEN+'β'+color.END, color.BLUE+'γ'+color.END]
    vectors = [directions[0],directions[1],directions[2],directions[0],directions[1]]
    for i,n in enumerate(names):
        dotProduct = np.dot(vectors[i+1],vectors[i+2])
        dotProduct /= np.linalg.norm(vectors[i+1])
        dotProduct /= np.linalg.norm(vectors[i+2])
        angle = np.arccos(np.clip(dotProduct,-1,1))/np.pi
        data = "{:s} = {:= 5.3f}π = {:3d}°".format(n,angle,int(180.0*angle))
        print('|'+data.center(linewidth+len(n)-1)+"|")

    print("+"+linewidth*'-'+"+")



def print_moments(moments,linewidth=42):
    if(linewidth < 42):
        print("Too low linewidth!")

    colors = [color.BLUE, color.DARKCYAN, color.CYAN]
    print("+"+linewidth*'-'+"+")
    print("|"+"Magnetic moments read:".center(linewidth)+"|")
    print('| ',end='')
    sumOfChars = 1
    for i,moment in enumerate(moments):
        momentStr = "{:=+3.1f} ".format(moment)
        if sumOfChars > linewidth-len(momentStr)-1:
            print(' '*(linewidth-sumOfChars),end='')
            print('|')
            print('| ',end='')
            sumOfChars = 1
        print(colors[i%3]+momentStr+color.END,end='')
        sumOfChars += len(momentStr) 
    print(' '*(linewidth-sumOfChars),end='')
    print('|')
    print("+"+linewidth*'-'+"+")


def print_label(text,linewidth=42):
    if len(text) > linewidth - 2:
        label = text[:linewidth-5]+"..."
    else:
        label = text.center(linewidth-2)
    label = color.BF+color.YELLOW+label+color.END
    print("+"+linewidth*'-'+"+")
    print("| "+label+" |")
    print("+"+linewidth*'-'+"+")
