class Color:
    """
        From:
        https://stackoverflow.com/questions/8924173/how-do-i-print-bold-text-in-python
    """
    bold  = '\033[1m'
    italic= '\033[3m'
    underl= '\033[4m'
    blink = '\033[5m'
    end   = '\033[0m'
    inv   = '\033[7m'
    hid   = '\033[8m'
    default       = "\033[39m"
    black         = "\033[30m"
    darkred       = "\033[31m"
    darkgreen     = "\033[32m"
    darkyellow    = "\033[33m"
    darkblue      = "\033[34m"
    darkmagenta   = "\033[35m"
    darkcyan      = "\033[36m"
    gray          = "\033[37m"
    darkgray      = "\033[90m"
    red           = "\033[91m"
    green         = "\033[92m"
    yellow        = "\033[93m"
    blue          = "\033[94m"
    magenta       = "\033[95m"
    cyan          = "\033[96m"
    white         = "\033[97m"

from sys import stdout
from copy import copy
class Standard:
    values = {'linewidth' : 88,
              'end'       : '\n',
              'stream'    : stdout}

    @staticmethod
    def fix(**kwargs):
        result = copy(Standard.values)
        result.update(kwargs)
        return result

def safe_update(primary,secondary):
    secondary.update(primary)
    return secondary

def print_vector(vector,**kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'vectorStyle' : Color.darkcyan})
    txt = ' [ ' + kwargs['vectorStyle']
    try:
        txt += len(vector)*' {:=10.5f}'.format(*vector,)
    except IndexError:
        txt += str(vector)
    txt += Color.end + ' ] '
    output  = "|"
    output += txt.center(kwargs['linewidth']+len(kwargs['vectorStyle']+Color.end))
    output += "|"+kwargs['end']
    kwargs['stream'].write(output)

def print_atom(atom,**kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'vectorStyle'  : Color.end,
                                 'elementStyle' : Color.bold+Color.darkyellow,
                                 'center'       : False})
    output  = kwargs['elementStyle'] + str(atom[0]) + Color.end
    output += ' [ ' + kwargs['vectorStyle']
    output += '{:= 10.5f} {:= 10.5f} {:= 10.5f}'.format(*atom[1],)
    output += Color.end + ' ] '
    if kwargs['center']:
        offset = len(kwargs['elementStyle'] + Color.end + kwargs['vectorStyle']+Color.end)
        output = output.center(kwargs['linewidth']+offset)
        output = "|" + output
        output += "|"
    kwargs['stream'].write(output+kwargs['end'])

def print_case(atom,**kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'wyckoffPosition' : " ",
                                 'distance'        : -1.0,
                                 'caseStyle'       : Color.italic+Color.cyan,
                                 'numberStyle'     : Color.italic,
                                 'caseID'          : 1,
                                 'atomID'          : 1,
                                 'distanceStyle'   : Color.end})
    output  = "Case " + kwargs['caseStyle']+"{:3d}".format(kwargs['caseID'])+Color.end
    output += " atom No. " + kwargs['numberStyle']+"{:3d}".format(kwargs['atomID'])+Color.end
    output += kwargs['distanceStyle']
    output += " @ {:1s} & {:6.2f} Å | ".format(kwargs['wyckoffPosition'],kwargs['distance'])
    output += Color.end
    kwargs['stream'].write(output)
    print_atom(atom,end=kwargs['end'])

def format_record(record,**kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'number' : 0})
    name   = "{:4s}:   ".format(record)
    num    = "#{:4d}".format(kwargs['number'])
    label  = kwargs['delimiter']    \
           + kwargs['elementStyle'] \
           + name                   \
           + Color.end              \
           + kwargs['numberStyle']  \
           + num                    \
           + Color.end
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
                              +Color.end+"|"+'\n')
    sumOfAtoms = 0
    for record in report:
        label,offset = format_record(record,number=report[record],**kwargs)
        kwargs['stream'].write('|'+label.center(kwargs['linewidth']+offset)+'|'+'\n')
        sumOfAtoms += report[record]

    label,offset = format_record("cell: ",number=sumOfAtoms,**kwargs)
    kwargs['stream'].write("|"+((kwargs['linewidth']//4)*'-').center(kwargs['linewidth'])+"|"+'\n')
    kwargs['stream'].write('|'+label.center(kwargs['linewidth']+offset)+'|'+'\n')

def print_axes(directions,**kwargs):
    kwargs['stream'].write("|"+Color.bold
                              +'Crystal axes:'.center(kwargs['linewidth'])
                              +Color.end+"|"+'\n')
    names = color_names('a⃗','b⃗','c⃗')
    for n,d in zip(names,directions):
        data = "{:s} = [ {:= 10.5f} {:= 10.5f} {:= 10.5f} ]".format(n,*d,)
        kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(n)-1)+"|"+'\n')

def print_directions(directions,**kwargs):
    kwargs['stream'].write("|"+kwargs['labelStyle']
                              +'Crystal directions:'.center(kwargs['linewidth'])
                              +Color.end+"|"+'\n')
    names = color_names('a','b','c')
    for name,direction in zip(names,directions):
        data = "{:s} = {:= 8.5f} Å".format(name,np.linalg.norm(direction))
        kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(name)-1)+"|"+'\n')

def print_angles(directions,**kwargs):
    kwargs['stream'].write("|"+kwargs['labelStyle']
                              +'Crystal angles:'.center(kwargs['linewidth'])
                              +Color.end+"|"+'\n')

    names = color_names('α','β','γ')
    tmpCycle = cycle(np.flipud(directions))
    for name in names:
        vectorJ = next(tmpCycle)
        vectorI = next(tmpCycle)
        dotProduct  = np.dot(vectorI,vectorJ)
        dotProduct /= np.linalg.norm(vectorI)
        dotProduct /= np.linalg.norm(vectorJ)
        angle = np.arccos(np.clip(dotProduct,-1,1))/np.pi
        data = "{:s} = {:= 5.3f}π = {:3d}°".format(name,angle,int(np.around(180.0*angle)))
        kwargs['stream'].write('|'+data.center(kwargs['linewidth']+len(name)-1)+"|"+'\n')

def empty_line(**kwargs):
    kwargs['stream'].write("|"+kwargs['linewidth']*' '+"|"+'\n')

def line(**kwargs):
    kwargs['stream'].write("+"+kwargs['linewidth']*'-'+"+"+'\n')

from itertools import cycle
def color_names(*args):
    colors = cycle([Color.red,Color.green,Color.blue])
    return [Color.bold+next(colors)+str(arg)+Color.end for arg in args]

import numpy as np
def print_crystal(cell, **kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'directions'   : np.eye(3),
                                 'labelStyle'   : Color.bold,
                                 'elementStyle' : Color.bold+Color.darkyellow,
                                 'numberStyle'  : Color.italic,
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
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'colors' : cycle([Color.darkgreen,
                                                   Color.darkblue,
                                                   Color.darkmagenta]),
                                 'elementStyle' : Color.bold+Color.gray,
                                 'labelStyle'   : Color.darkgray})

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
                              +elementStr + Color.end
                              +kwargs['labelStyle']
                              +numberStr
                              +next(kwargs['colors'])+momentStr
                              +Color.end)
        sumOfChars += length
    kwargs['stream'].write(' '*(kwargs['linewidth']-sumOfChars))
    kwargs['stream'].write('|'+'\n')

    line(**kwargs)


def print_label(text,**kwargs):
    kwargs = Standard.fix(**kwargs)
    kwargs = safe_update(kwargs,{'elementStyle' : Color.bold+Color.darkgreen,
                                 'labelStyle'   : Color.bold+Color.darkyellow,
                                 'atoms'        : None,
                                 'vectorStyle'  : Color.end})
    if len(text) > kwargs['linewidth'] - 2:
        label = text[:kwargs['linewidth']-5]+"..."
    else:
        label = text.center(kwargs['linewidth']-2)
    label = kwargs['labelStyle']+label+Color.end
    line(**kwargs)
    kwargs['stream'].write("| "+label+" |"+'\n')

    if kwargs['atoms'] is not None:
        for atom in kwargs['atoms']:
            print_atom(atom,center=True,**kwargs)
    line(**kwargs)
