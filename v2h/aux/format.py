class color:
    """
        From:
        https://stackoverflow.com/questions/8924173/how-do-i-print-bold-text-in-python
    """    
    BF    = '\033[1m'
    IT    = '\033[3m'
    UN    = '\033[4m'
    BLINK = '\033[5m'
    END   = '\033[0m'
    INV   = '\033[7m'
    HID   = '\033[8m'
    DEFAULT       = "\033[39m"
    BLACK         = "\033[30m"
    DARKRED       = "\033[31m"
    DARKGREEN     = "\033[32m"
    DARKYELLOW    = "\033[33m"
    DARKBLUE      = "\033[34m"
    DARKMAGENTA   = "\033[35m"
    DARKCYAN      = "\033[36m"
    GRAY          = "\033[37m"
    DARKGRAY      = "\033[90m"
    RED           = "\033[91m"
    GREEN         = "\033[92m"
    YELLOW        = "\033[93m"
    BLUE          = "\033[94m"
    MAGENTA       = "\033[95m"
    CYAN          = "\033[96m"
    WHITE         = "\033[97m"

from sys import stdout
def print_atom(atom,end='\n',vector=color.END,elementStyle=color.BF+color.YELLOW,stream=stdout):
    output  = elementStyle + str(atom[0]) + color.END
    output += ' [ ' + vector
    output += '{:= 10.5f} {:= 10.5f} {:= 10.5f}'.format(*atom[1],)
    output += color.END + ' ] '
    stream.write(output+end)

def print_case(caseID, atom, atomID,
               wyckoffPosition, distance,
               caseStyle=color.IT+color.CYAN,
               numberStyle=color.IT,
               distanceStyle=color.END, end='\n',stream=stdout):
    output  = "Case "
    output += caseStyle+"{:3d}".format(caseID)+color.END
    output += " atom No. "
    output += numberStyle+"{:3d}".format(atomID)+color.END
    output += distanceStyle+" @ {:1s} & {:6.2f} Å | ".format(wyckoffPosition,distance)+color.END
    stream.write(output)
    print_atom(atom,end=end)

import numpy as np
def print_crystal(directions, cell,
                  linewidth=44,atomNames=None,
                  labelStyle=color.BF,
                  elementStyle=color.BF+color.YELLOW,
                  numberStyle=color.IT
                  ,stream=stdout):
    if linewidth < 44:
        stream.write("Too low linewidth!"+'\n')
    
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

    stream.write("+"+linewidth*'-'+"+"+'\n')
    stream.write("|"+labelStyle+"Composition:".center(linewidth)+color.END+"|"+'\n')
    sumOfAtoms = 0
    delimiter = ' | '
    for record in report:
        name   = "{:4s}:   ".format(record)
        number = "#{:4d}".format(report[record])
        label  = delimiter+elementStyle+name+color.END+numberStyle+number+color.END
        offset = len(label)-len(name)-len(number)-len(delimiter)
        stream.write('|'+label.center(linewidth+offset)+'|'+'\n')
        sumOfAtoms += report[record]

    delimiter = '+| '
    name   = "cell:   "
    number = "#{:4d}".format(sumOfAtoms)
    label  = delimiter+elementStyle+name+color.END+numberStyle+number+color.END
    offset = len(label)-len(name)-len(number)-len(delimiter)
    stream.write("|"+((len(name)+len(number)+len(delimiter)+1)*'-').center(linewidth)+"|"+'\n')
    stream.write('|'+label.center(linewidth+offset)+'|'+'\n')
    

    stream.write("|"+linewidth*' '+"|"+'\n')
    stream.write("|"+color.BF+'Crystal axes:'.center(linewidth)+color.END+"|"+'\n')

    names = [color.BF+color.RED  +'a⃗'+color.END,
             color.BF+color.GREEN+'b⃗'+color.END,
             color.BF+color.BLUE +'c⃗'+color.END]
    for n,d in zip(names,directions):
        data = "{:s} = [ {:= 10.5f} {:= 10.5f} {:= 10.5f} ]".format(n,*d,)
        stream.write('|'+data.center(linewidth+len(n)-1)+"|"+'\n')

    stream.write("|"+linewidth*' '+"|"+'\n')
    stream.write("|"+labelStyle+'Crystal directions:'.center(linewidth)+color.END+"|"+'\n')

    names = [color.RED+'a'+color.END, color.GREEN+'b'+color.END, color.BLUE+'c'+color.END]
    for n,d in zip(names,directions):
        data = "{:s} = {:= 8.5f} Å".format(n,np.linalg.norm(d))
        stream.write('|'+data.center(linewidth+len(n)-1)+"|"+'\n')

    stream.write("|"+linewidth*' '+"|"+'\n')
    stream.write("|"+labelStyle+'Crystal angles:'.center(linewidth)+color.END+"|"+'\n')

    names = [color.RED+'α'+color.END, color.GREEN+'β'+color.END, color.BLUE+'γ'+color.END]
    vectors = [directions[0],directions[1],directions[2],directions[0],directions[1]]
    for i,n in enumerate(names):
        dotProduct = np.dot(vectors[i+1],vectors[i+2])
        dotProduct /= np.linalg.norm(vectors[i+1])
        dotProduct /= np.linalg.norm(vectors[i+2])
        angle = np.arccos(np.clip(dotProduct,-1,1))/np.pi
        data = "{:s} = {:= 5.3f}π = {:3d}°".format(n,angle,int(180.0*angle))
        stream.write('|'+data.center(linewidth+len(n)-1)+"|"+'\n')

    stream.write("+"+linewidth*'-'+"+"+'\n')



def print_moments(moments,cell=None,atomNames=None,linewidth=44,
                  colors=[color.DARKYELLOW,
                          color.DARKGREEN,
                          color.DARKBLUE,
                          color.DARKMAGENTA],
                  elementStyle=color.BF+color.GRAY,
                  labelStyle=color.DARKGRAY,stream=stdout):
    if linewidth < 44:
        stream.write("Too low linewidth!"+'\n')

    stream.write("+"+linewidth*'-'+"+"+'\n')
    stream.write("|"+"Magnetic moments read:".center(linewidth)+"|"+'\n')
    stream.write('| ')
    sumOfChars = 1
    for i,moment in enumerate(moments):
        elementStr = ''
        numberStr  = ''
        if cell is not None and atomNames is not None:
            elementStr = "{:2s}".format(atomNames[cell[i][0]])
            numberStr  = "({:3d}): ".format(i+1)
        momentStr = "{:=+3.1f} ".format(moment)
        length = len(elementStr)+len(numberStr)+len(momentStr)
        if sumOfChars > linewidth-length-1:
            stream.write(' '*(linewidth-sumOfChars))
            stream.write('|'+'\n')
            stream.write('| ')
            sumOfChars = 1
        stream.write(elementStyle+elementStr+color.END+labelStyle
                +numberStr+colors[i%len(colors)]+momentStr+color.END)
        sumOfChars += length
    stream.write(' '*(linewidth-sumOfChars))
    stream.write('|'+'\n')
    stream.write("+"+linewidth*'-'+"+"+'\n')


def print_label(text,linewidth=44,labelStyle=color.BF+color.YELLOW,stream=stdout):
    if len(text) > linewidth - 2:
        label = text[:linewidth-5]+"..."
    else:
        label = text.center(linewidth-2)
    label = labelStyle+label+color.END
    stream.write("+"+linewidth*'-'+"+"+'\n')
    stream.write("| "+label+" |"+'\n')
    stream.write("+"+linewidth*'-'+"+"+'\n')
