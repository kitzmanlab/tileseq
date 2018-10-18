# cython: c_string_type=unicode, c_string_encoding=utf8


cimport cython
 # from libc.stdlib cimport malloc, free
cimport numpy as np
import numpy as np

import os
import sys

from pysam.libcsamfile cimport AlignedSegment

cpdef cigToStr( cig ):
    cdef str x=''
    for ce in cig:
        if ce[0]==0:
            x+='%dM'%ce[1]
        elif ce[0]==1:
            x+='%dI'%ce[1]
        elif ce[0]==2:
            x+='%dD'%ce[1]
        elif ce[0]==3:
            x+='%dN'%ce[1]
        elif ce[0]==4:
            x+='%dS'%ce[1]
        elif ce[0]==5:
            x+='%dH'%ce[1]
        elif ce[0]==6:
            x+='%dP'%ce[1]
    return x

cpdef convertSeqToArray( 
    const unsigned char[:] seq, 
    int[::1] ar ):

    cdef int N = len(seq)
    for 0 <= i < N:
        if seq[i]=='A' or seq[i]=='a':
            ar[i]=0
        elif seq[i]=='C' or seq[i]=='c':
            ar[i]=1
        elif seq[i]=='G' or seq[i]=='g':
            ar[i]=2
        elif seq[i]=='T' or seq[i]=='t':
            ar[i]=3
        elif seq[i]=='N' or seq[i]=='n':
            ar[i]=4            


cpdef realignFrameaware(
    bamIn,
    bamOut,
    int maxReadLen,
    const unsigned char[:] seqRef,
    int cozCdsStart,
    int cozCdsEnd,
    int matchScore,
    int mismatchScore,
    int gapopenScoreInframe,
    int gapopenScoreOutframe,
    int gapextendScore,
    list ldebug
    ):

    cdef AlignedSegment rd
    cdef int hasDel
    cdef int cozRdStart,cozRdEnd

    cdef int[::1] arseqRef = -9999999*np.ones( len(seqRef), dtype=np.int32 )
    convertSeqToArray( seqRef, arseqRef )

    cdef int[::1] arseqRd = -9999999*np.ones( maxReadLen, dtype=np.int32 )

    cdef int maxSizeOnRef = maxReadLen + 100

    cdef int frameAtS2start

    cdef int[:,:,::1] scoreMtx = -9999999*np.ones( (maxReadLen,maxSizeOnRef,6), dtype=np.int32 )

    cdef int rlen

    for rd in bamIn:
        hasDel=0

        for cpart in rd.cigar:
            if cpart[0]==2:
                hasDel=1
                break

        if hasDel==0:
            # no deletion; passthru
            bamOut.write(rd)
        else:

            # a handful of reads seem to have deletions at their starts - deletion rather than softmasked
            # mapping position vs ref includes these (i.e., if 5D then the first mapped bases are not until pos+5)
            # but the .query does NOT include these bases
            # 
            #  -> simply skip these reads.
            if rd.cigar[0][0] == 2:
                continue

            cozRdStart=rd.pos
            cozRdEnd=cozRdStart+rd.alen-1

            # in case there is a really large deletion, (permanently) expand the score martix to accomodate.
            while cozRdEnd-cozRdStart+1 > maxSizeOnRef:
                maxSizeOnRef += 100
                scoreMtx = -9999999*np.ones( (maxReadLen,maxSizeOnRef,6), dtype=np.int32 )

            arseqRd[:]=-999999
            rlen = rd.qlen
            convertSeqToArray( rd.query.encode('utf-8'), arseqRd )

            frameAtS2start = (cozRdStart - cozCdsStart)%3

            #print 'rd name ',rd.qname
            #print 'rd.rdStart ', cozRdStart, ' frameAtS2start ', frameAtS2start
            #print 'rd.seq ',rd.query
            #print ' ofsCdsBounds:  ', cozCdsStart-cozRdStart,'  ', cozCdsEnd-cozRdStart
            #print 'rd.alen, passed s2len ',( rd.alen, arseqRef[ cozRdStart:cozRdEnd+1 ].shape[0] )
            #print ' ',(scoreMtx.shape[0],scoreMtx.shape[1],scoreMtx.shape[2])


            scoreMtx[:,:,:] = 0

            populateScoreMatrix_GlobalPenalizeOutOfFrame( 
                scoreMtx,
                arseqRd,
                arseqRef[ cozRdStart:cozRdEnd+1 ],
                rlen,
                rd.alen,
                cozCdsStart-cozRdStart,
                cozCdsEnd-cozRdStart,
                matchScore,
                mismatchScore,
                gapopenScoreOutframe,
                gapopenScoreInframe,
                gapopenScoreOutframe,
                gapextendScore)
    
            cig = tracebkToCigar_Global( 
                scoreMtx, 
                arseqRd, 
                arseqRef[ cozRdStart:cozRdEnd+1 ],
                rlen,
                rd.alen )

            # do we need to clip the read at all?
            cigReplace = []
            for ce in rd.cigar:
                if ce[0]>=4:
                    cigReplace.append(ce)
                else:
                    break
            cigReplace = cigReplace+cig
            for ce in rd.cigar[::-1]:
                if ce[0]>=4:
                    cigReplace.append(ce)
                else:
                    break

            # strip off the NM tag 
            # place the old cigar tag in 'X9'
            rd.tags = [ t for t in rd.tags if t[0]!='NM' ] + [('X9',cigToStr( rd.cigar ))]

            #ldebug.append( ( rd.query,
            #                 seqRef[cozRdStart:cozRdEnd+1],
            #                 np.array(scoreMtx.copy()),
            #                 rd.cigar,
            #                 cig ) )

            rd.cigar = cigReplace
            bamOut.write(rd)






#@cython.boundscheck(False)
#@cython.wraparound(False)
cdef populateScoreMatrix_GlobalPenalizeOutOfFrame(  
    int[:,:,::1] matrix, 
    int[::1] seq1, 
    int[::1] seq2, 
    int s1len, int s2len,
    int ofsCdsStartInS2,
    int ofsCdsEndInS2,
    int match, 
    int mismatch, 
    int s1Gapopen,
    int s2GapopenInframe, 
    int s2GapopenOutofframe,
    int gapextend):

    """

    ofsCdsStartIn2 is the 0-based offset of the cds start, relative to the start of 

    matrix should be at least (s1len+1, s2len+1, == 4 )
      dim 2 (last dimension),  0 --> M[i,j]
                               1 --> I[i,j]
                               2 --> J[i,j]
                               3 --> traceback M [ 1=D, 2=V, 3=H ]
                               4 --> traceback M [ 1=D, 2=V, 3=H ]
                               5 --> traceback M [ 1=D, 2=V, 3=H ]


    """

    cdef int i,j

    cdef int matchscore = 0

    cdef int upleftscore = 0
    cdef int leftscore = 0
    cdef int upscore = 0

    cdef int temp1, temp2, temp3, temp4, gapuse, direc

    cdef int tbM,tbI,tbJ

    # initialize M
    matrix[ 0, 0, 0 ] = 0
    matrix[ 0, 1:, 0 ] = -99999999
    matrix[ 1:, 0, 0 ] = -88888888

    # init I (vertical gaps in s1 === insertion relative to ref, s2)
    matrix[ 0, 0, 1 ] = -77777777
    matrix[ 0, 1:, 1] = -77777777
    matrix[ 1, 0, 1 ] = s1Gapopen
    for 2 <= i <= s1len:
        matrix[ i, 0, 1 ] = matrix[ i-1, 0, 1 ] + gapextend

    matrix[ 1:, 0, 3 ] = 2

    # init J (vertical gaps in s1 === insertion relative to ref, s2)
    matrix[ 0, 0, 2 ] = -66666666
    matrix[ 1:, 0, 2] = -66666666
    matrix[ 0, 1, 2 ] = s2GapopenOutofframe
    for 2 <= j <= s2len:
        matrix[ 0, j, 2 ] = matrix[ 0, j-1, 2 ] + gapextend

    matrix[ 0, 1:, 3 ] = 3

    # row by row, starting at row 1
    for i in range(1,s1len+1):
        # column by column, starting at column 1
        for j in range(1,s2len+1):

            if seq1[i-1]==seq2[j-1]:
                matchscore = match
            else:
                matchscore = mismatch

            # compute I[i,j]
            temp1 = matrix[ i - 1, j, 0 ] + s1Gapopen
            temp2 = matrix[ i - 1, j, 1 ] + gapextend
            matrix[ i, j, 1 ] = max( temp1, temp2 )
            if temp1 >= temp2:
                tbI = 1
            else:
                tbI = 2

            matrix[ i, j, 4 ] = tbI

            # compute J[i,j]

            # if chosen, the gap open base will be j+1 [1based] = j [0based]
            if ( ( (j-1) - ofsCdsStartInS2) % 3 == 0) and (j-1 >= ofsCdsStartInS2) and (j-1 <= ofsCdsEndInS2):
                gapuse = s2GapopenInframe
            else:
                gapuse = s2GapopenOutofframe

            temp3 = matrix[ i, j - 1, 0 ] + gapuse
            temp4 = matrix[ i, j - 1, 2 ] + gapextend
            matrix[ i, j, 2 ] = max( temp3, temp4 )

            if temp3 >= temp4:
                tbJ = 1
            else:
                tbJ = 3

            matrix[ i, j, 5 ] = tbJ


            upleftscore = matrix[ i-1, j-1, 0 ] + matchscore
            upscore = matrix[ i-1, j-1, 1 ] + matchscore
            leftscore = matrix[ i-1, j-1, 2 ] + matchscore

            #3 --> traceback  [  1=D, 2=V, 3=H ]

            bestscore = max( upscore, upleftscore, leftscore )
            matrix[ i, j, 0 ] = bestscore

            if upleftscore >= leftscore and upleftscore >= upscore:
                tbM = 1
            elif leftscore > upleftscore and leftscore >= upscore:
                tbM = 3
            elif upscore > upleftscore and upscore >= leftscore:
                tbM = 2
            else:
                assert 1==0, (i,j)

            matrix[ i, j, 3 ] = tbM

            """
            if matrix[ i, j, 1 ] >= matrix[ i, j, 0 ] and matrix[ i, j, 1 ] >= matrix[ i, j, 2 ]:
                direc = tbI
            elif matrix[ i, j, 0 ] >= matrix[ i, j, 1 ] and matrix[ i, j, 0 ] >= matrix[ i, j, 2 ]:
                direc = tbM
            elif matrix[ i, j, 2 ] >= matrix[ i, j, 0 ] and matrix[ i, j, 2 ] >= matrix[ i, j, 1 ]:
                direc = tbJ
            else:
                assert 1==0

            matrix[ i, j, 3 ] = direc"""

            """
            if abs(i-j)<5:
                print '\n\n i=%d j=%d'%(i,j), ' s_i,j= ',matchscore
                print '   s2gapopenhere=',gapuse
                print ' upleft=',upleftscore, '  left=',leftscore, ' upscore=',upscore
                print ' M[i,j]=%d    I[i,j]=%d   J[i,j]=%d'%( matrix[i,j,0], matrix[i,j,1], matrix[i,j,2] )
                print 'tb MIJ ',tbM,tbI,tbJ
            """



#@cython.boundscheck(True)
#@cython.wraparound(True)
cdef tracebkToCigar_Global( int[:,:,::1] matrix, 
                            int[::1] seq1, 
                            int[::1] seq2, 
                            int s1len, 
                            int s2len ):

    #print 'tracebkToCigar_Global'

    cdef list cig = []

    cdef int row, col

    row = s1len
    col = s2len
    
    cdef int state=1, stateNext=0
    cdef int numInCurStart=0

    # trace back to beginning
    while ( row>0 and col>0 ):

        #print row,col,list(matrix[row,col,:])

        if state == 1 :
            #print 'M'
            stateNext = matrix[ row, col, 3+state-1 ]

            row -= 1
            col -= 1

        elif state == 2:
            #print 'I'
            stateNext = matrix[ row, col, 3+state-1 ]

            row -= 1

        elif state == 3:
            #print 'D'
            stateNext = matrix[ row, col, 3+state-1 ]

            col -= 1

        numInCurStart += 1

        if stateNext != state:
            cig.append( ( state - 1, numInCurStart ) )
            #print ' -->  ',state-1,numInCurStart
            numInCurStart = 0 

        state = stateNext
    
    if numInCurStart != 0 :
        cig.append( (state - 1, numInCurStart) )

    cig=cig[::-1]

    return cig