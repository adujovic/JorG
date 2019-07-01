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
from copy import copy
class standard:
    values = {'linewidth' : 88,
              'end'       : '\n',
              'stream'    : stdout}

    @staticmethod
    def fix(**kwargs):
        result = copy(standard.values)
        result.update(kwargs)
        return result

def safe_update(primary,secondary):
    secondary.update(primary)
    return secondary

def print_vector(vector,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'vectorStyle' : color.DARKCYAN})
    txt = ' [ ' + kwargs['vectorStyle']
    for x in vector:
        txt += ' {:=10.5f}'.format(x)
    txt += color.END + ' ] '
    output  = "|"
    output += txt.center(kwargs['linewidth']+len(kwargs['vectorStyle']+color.END))
    output += "|"+kwargs['end']
    kwargs['stream'].write(output)

def print_atom(atom,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'vectorStyle'  : color.END,
                                 'elementStyle' : color.BF+color.DARKYELLOW,
                                 'center'       : False})
    output  = kwargs['elementStyle'] + str(atom[0]) + color.END
    output += ' [ ' + kwargs['vectorStyle']
    output += '{:= 10.5f} {:= 10.5f} {:= 10.5f}'.format(*atom[1],)
    output += color.END + ' ] '
    if kwargs['center']:
        offset = len(kwargs['elementStyle'] + color.END + kwargs['vectorStyle']+color.END)
        output = output.center(kwargs['linewidth']+offset)
        output = "|" + output
        output += "|"
    kwargs['stream'].write(output+kwargs['end'])

def print_case(atom,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'wyckoffPosition' : " ",
                                 'distance'        : -1.0,
                                 'caseStyle'       : color.IT+color.CYAN,
                                 'numberStyle'     : color.IT,
                                 'caseID'          : 1,
                                 'atomID'          : 1,
                                 'distanceStyle'   : color.END})
    output  = "Case " + kwargs['caseStyle']+"{:3d}".format(kwargs['caseID'])+color.END
    output += " atom No. " + kwargs['numberStyle']+"{:3d}".format(kwargs['atomID'])+color.END
    output += kwargs['distanceStyle']
    output += " @ {:1s} & {:6.2f} Å | ".format(kwargs['wyckoffPosition'],kwargs['distance'])
    output += color.END
    kwargs['stream'].write(output)
    print_atom(atom,end=kwargs['end'])

def format_record(record,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'number' : 0})
    name   = "{:4s}:   ".format(record)
    num    = "#{:4d}".format(kwargs['number'])
    label  = kwargs['delimiter']    \
           + kwargs['elementStyle'] \
           + name                   \
           + color.END              \
           + kwargs['numberStyle']  \
           + num                    \
           + color.END
    offset = len(label)-len(name)-len(num)-len(kwargs['delimiter'])
    return label,offset

def print_composition(cell,**kwargs):
    report = {}
    for atom in cell:
        if atom[0] in report.keys():
            report[atom[0]] += 1
        else:
            report[atom[0]] = 1

    kwargs['stream'].write("+"+kwargs['linewidth']*'-'+"+"+'\n')
    kwargs['stream'].write("|"+kwargs['labelStyle']
                              +"Composition:".center(kwargs['linewidth'])
                              +color.END+"|"+'\n')
    sumOfAtoms = 0
    for record in report:
        label,offset = format_record(record,number=report[record],**kwargs)
        kwargs['stream'].write('|'+label.center(kwargs['linewidth']+offset)+'|'+'\n')
        sumOfAtoms += report[record]

    label,offset = format_record("cell: ",number=sumOfAtoms,**kwargs)
    kwargs['stream'].write("|"+((kwargs['linewidth']//4)*'-').center(kwargs['linewidth'])+"|"+'\n')
    kwargs['stream'].write('|'+label.center(kwargs['linewidth']+offset)+'|'+'\n')

def print_axes(directions,**kwargs):
    kwargs['stream'].write("|"+color.BF
                              +'Crystal axes:'.center(kwargs['linewidth'])
                              +color.END+"|"+'\n')
    names = color_names('a⃗','b⃗','c⃗')
    for n,d in zip(names,directions):
        data = "{:s} = [ {:= 10.5f} {:= 10.5f} {:= 10.5f} ]".format(n,*d,)
        kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(n)-1)+"|"+'\n')

def print_directions(directions,**kwargs):
   kwargs['stream'].write("|"+kwargs['labelStyle']
                             +'Crystal directions:'.center(kwargs['linewidth'])
                             +color.END+"|"+'\n')
   names = color_names('a','b','c')
   for n,d in zip(names,directions):
       data = "{:s} = {:= 8.5f} Å".format(n,np.linalg.norm(d))
       kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(n)-1)+"|"+'\n')

def print_angles(directions,**kwargs):
    kwargs['stream'].write("|"+kwargs['labelStyle']
                              +'Crystal angles:'.center(kwargs['linewidth'])
                              +color.END+"|"+'\n')

    names = color_names('α','β','γ')
    tmpCycle = cycle(np.flipud(directions))
    for n in names:
        vectorJ = next(tmpCycle)
        vectorI = next(tmpCycle)
        dotProduct  = np.dot(vectorI,vectorJ)
        dotProduct /= np.linalg.norm(vectorI)
        dotProduct /= np.linalg.norm(vectorJ)
        angle = np.arccos(np.clip(dotProduct,-1,1))/np.pi
        data = "{:s} = {:= 5.3f}π = {:3d}°".format(n,angle,int(np.around(180.0*angle)))
        kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(n)-1)+"|"+'\n')

def empty_line(**kwargs):
    kwargs['stream'].write("|"+kwargs['linewidth']*' '+"|"+'\n')

def line(**kwargs):
    kwargs['stream'].write("+"+kwargs['linewidth']*'-'+"+"+'\n')

from itertools import cycle
def color_names(*args):
    colors = cycle([color.RED,color.GREEN,color.BLUE])
    return [color.BF+next(colors)+str(arg)+color.END for arg in args]

import numpy as np
def print_crystal(cell, **kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'directions'   : np.eye(3),
                                 'labelStyle'   : color.BF,
                                 'elementStyle' : color.BF+color.DARKYELLOW,
                                 'numberStyle'  : color.IT,
                                 'delimiter'    : ' | '})

    print_composition(cell,**kwargs)
    empty_line(**kwargs)

    print_axes(**kwargs)
    empty_line(**kwargs)

    print_directions(**kwargs)
    empty_line(**kwargs)

    print_angles(**kwargs)
    line(**kwargs)

def check_line(sumOfChars,length,**kwargs):
    if sumOfChars > kwargs['linewidth']-length-1:
        kwargs['stream'].write(' '*(kwargs['linewidth']-sumOfChars))
        kwargs['stream'].write('|'+'\n')
        kwargs['stream'].write('| ')
        sumOfChars = 1
    return sumOfChars

def print_moments(moments,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'colors' : cycle([color.DARKGREEN,
                                                   color.DARKBLUE,
                                                   color.DARKMAGENTA]),
                                 'elementStyle' : color.BF+color.GRAY,
                                 'labelStyle'   : color.DARKGRAY})

    line(**kwargs)

    kwargs['stream'].write("|"+"Magnetic moments read:".center(kwargs['linewidth'])+"|"+'\n')
    kwargs['stream'].write('| ')
    sumOfChars = 1
    for i,moment in enumerate(moments):
        elementStr = ''
        numberStr  = ''
        if 'cell' in kwargs:
            elementStr = "{:2s}".format(kwargs['cell'][i][0])
            numberStr  = "({:3d}): ".format(i+1)
        momentStr = "{:=+3.1f} ".format(moment)
        length = len(elementStr)+len(numberStr)+len(momentStr)
        sumOfChars = check_line(sumOfChars,length,**kwargs)
        kwargs['stream'].write(kwargs['elementStyle']
                              +elementStr + color.END
                              +kwargs['labelStyle']
                              +numberStr
                              +next(kwargs['colors'])+momentStr
                              +color.END)
        sumOfChars += length
    kwargs['stream'].write(' '*(kwargs['linewidth']-sumOfChars))
    kwargs['stream'].write('|'+'\n')

    line(**kwargs)


def print_label(text,**kwargs):
    kwargs = standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'elementStyle' : color.BF+color.DARKGREEN,
                                 'labelStyle'   : color.BF+color.DARKYELLOW,
                                 'atoms'        : None,
                                 'vectorStyle'  : color.END})
    if len(text) > kwargs['linewidth'] - 2:
        label = text[:kwargs['linewidth']-5]+"..."
    else:
        label = text.center(kwargs['linewidth']-2)
    label = kwargs['labelStyle']+label+color.END
    line(**kwargs)
    kwargs['stream'].write("| "+label+" |"+'\n')

    if kwargs['atoms'] is not None:
        for atom in kwargs['atoms']:
            print_atom(atom,center=True,**kwargs)
    line(**kwargs)
