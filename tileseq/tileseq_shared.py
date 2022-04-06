
import argparse

def isWithin(cooCoord, cooRange, fIsInclusive=True):
    if fIsInclusive: return (cooCoord <= max(cooRange) and cooCoord >= min(cooRange))
    else:            return (cooCoord <  max(cooRange) and cooCoord >  min(cooRange))

def chrStartStop(s):
    return s.split(':')[0],int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1])
    
def chrTupStartStop(s):
    return s.split(':')[0],(int(s.split(':')[1].split('-')[0]),int(s.split(':')[1].split('-')[1]))

def coordConversion( typeFrom='[01]', typeTo='[00]' ):
    intvCorrectionFrom=[0,0]
    intvCorrectionFrom[0] = ( typeFrom[0]=='(' and 1 or
                          typeFrom[0]=='[' and 0 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 
    intvCorrectionFrom[1] = ( typeFrom[3]==')' and 0 or
                          typeFrom[3]==']' and 1 or
                          0 ) +\
                        ( typeFrom[1]=='0' and 0 or
                          typeFrom[1]=='1' and -1 or
                          0 ) 

    intvCorrectionTo=[0,0]
    intvCorrectionTo[0] = ( typeTo[0]=='(' and 1 or
                          typeTo[0]=='[' and 0 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 
    intvCorrectionTo[1] = ( typeTo[3]==')' and 0 or
                          typeTo[3]==']' and 1 or
                          0 ) +\
                        ( typeTo[1]=='0' and 0 or
                          typeTo[1]=='1' and -1 or
                          0 ) 

    return ( intvCorrectionFrom[0]-intvCorrectionTo[0], 
             intvCorrectionFrom[1]-intvCorrectionTo[1] )


def chrStartStopArg_11incl_to_00incl(s):
    try:
        chrom,start,stop=chrStartStop(s)
        return (chrom,
        	    coordConversion('[11]','[00]')[0]+start,
        	    coordConversion('[11]','[00]')[0]+stop )
    except:
        raise argparse.ArgumentTypeError('range must be chrom:start-stop')


import Bio.Seq
import Bio.SeqUtils 

mTransTbl = {}
for x in 'ACGT': 
    for y in 'ACGT':
        for z in 'ACGT':
            mTransTbl['%s%s%s'%(x,y,z)] = Bio.Seq.Seq( '%s%s%s'%(x,y,z) ).translate().__str__()

maa1to3 = {}
for aa in 'ACDEFGHIKLMNPQRSTVWY*': maa1to3[aa] = Bio.SeqUtils.seq3(aa)
