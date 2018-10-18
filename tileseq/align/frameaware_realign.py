import sys
from optparse import OptionParser
from collections import defaultdict

import copy
 
from jkutils import *

import math

import numpy as np

import pysam

import collections

import os
import sys

import frameaware_realign_core

def main():   
    opts=OptionParser()

    opts.add_option('','--faRef',dest='faRef')  # reference of subassembly casette, against which SA reads are aligned
    opts.add_option('','--refName',dest='refName')
    
    opts.add_option('','--inBam',dest='inBam')
    opts.add_option('','--outBam',dest='outBam')

    opts.add_option('','--corngCds',dest='corngCds')

    (o,args)=opts.parse_args()

    corzCds = ( int(o.corngCds.split(',')[0])-1, int(o.corngCds.split(',')[1])-1 )

    bamIn = pysam.Samfile( o.inBam, 'rb' )
    bamOut = pysam.Samfile( o.outBam, 'wb', template=bamIn )

    faRef = pysam.Fastafile( o.faRef )
    seqRef = faRef.fetch( o.refName ).upper().encode('utf-8')

    x=[]

    frameaware_realign_core.realignFrameaware(
        bamIn,
        bamOut,
        maxReadLen=600,
        seqRef=seqRef,
        cozCdsStart=corzCds[0],
        cozCdsEnd=corzCds[1],
        matchScore=1,
        mismatchScore=-1,
        gapopenScoreInframe=-2,
        gapopenScoreOutframe=-3,
        gapextendScore=-1,
        ldebug=x )

    # for i in x:
    #     print i[0],i[1],i[3],i[4]


if __name__ == '__main__':                
    main()