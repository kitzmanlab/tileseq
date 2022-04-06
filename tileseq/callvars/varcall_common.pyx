cimport cython
cimport numpy as np
import numpy as np

from pysam.libcfaidx cimport FastaFile
from pysam.libcalignedsegment cimport AlignedSegment

# cython: c_string_type=str, c_string_encoding=ascii
cdef class Variant:
    
    cdef GappedAlign galign

    cdef readonly int coz_ref_start, coz_ref_end
    cdef readonly int coz_query_start, coz_query_end

    cdef readonly int gofs_start, gofs_end

    cdef readonly str gapped_query_seq
    cdef readonly str gapped_ref_seq

    cdef readonly int n_match, n_mm, n_ins, n_del, n_mm_tv, n_mm_ti

    cdef readonly int[:] gapped_query_qual
    cdef readonly int[2] flanking_quals  #-1 means there is no flanking base
    
    cpdef int get_min_bq_incl_flanks(self):

        # print(np.array(self.gapped_query_qual))
        # print(np.array(self.flanking_quals)) 

        cdef int themin=int(1e6)
        cdef int i 
        for i in range(0, len(self.gapped_query_qual)):
            if self.gapped_query_qual[i]>=0:
                themin = min(themin, self.gapped_query_qual[i])
        if self.flanking_quals[0]>=0:
            themin=min(themin, self.flanking_quals[0])
        if self.flanking_quals[1]>=0:
            themin=min(themin, self.flanking_quals[1])
        return themin

    cpdef float get_mean_bq_incl_flanks(self):
        cdef int agg=0, n=0
        cdef int i 
        for i in range(0, len(self.gapped_query_qual)):
            if self.gapped_query_qual[i]>=0:
                agg += self.gapped_query_qual[i]
                n+=1
        if self.flanking_quals[0]>=0:
            agg+=self.flanking_quals[0]
            n+=1
        if self.flanking_quals[1]>=0:
            agg+=self.flanking_quals[1]
            n+=1
        if n>0:
            return float(agg)/float(n)
        else:
            return -1

    def __str__(self):
        return "\n".join(
            ["%d match/%d MM/%d INS/%d DEL"%( self.n_match, self.n_mm, self.n_ins, self.n_del ),
             "ref:%d-%d"%(self.coz_ref_start, self.coz_ref_end),
             "query:%d-%d"%(self.coz_query_start, self.coz_query_end),
             "gapped:%d-%d"%(self.gofs_start, self.gofs_end),
              "REF  %s"%(self.gapped_ref_seq),
              "READ %s"%(self.gapped_query_seq) ] )


    def __init__( self,
                  GappedAlign galign,
                  int gofs_start,
                  int gofs_end ):
    
        self.gapped_query_seq=""
        self.gapped_ref_seq=""

        cdef int gi

        self.n_match = self.n_mm = self.n_ins = self.n_del = 0

        self.n_mm_tv = self.n_mm_ti = 0

        self.galign = galign

        self.gofs_start=gofs_start
        self.gofs_end=gofs_end

        self.gapped_query_qual = np.zeros( (self.gofs_end-self.gofs_start+1), dtype=np.int32 )

        self.coz_ref_start = galign.coz_ref_start + galign.gali_ugofsRef[ gofs_start ]
        self.coz_ref_end = galign.coz_ref_start + galign.gali_ugofsRef[ gofs_end ]

        cdef str bpr, bpq

        # get flanking qualities
        self.flanking_quals[0]=-1
        gi = self.gofs_start - 1
        while gi >= 0 and (galign.gali_ugofsRead[gi]>=0):
            if galign.gali_nongapRdRefMatch[ gi, 0 ] != 0:
                self.flanking_quals[ 0 ] = <int>( galign.rd.query_qualities[ galign.gali_ugofsRead[gi] ] )
                break
            gi-=1

        self.flanking_quals[1]=-1
        gi = self.gofs_end + 1
        while gi <= galign.gali_len-1 and (galign.gali_ugofsRead[gi]<len(galign.query_seq)):
            if galign.gali_nongapRdRefMatch[ gi, 0 ] != 0:
                self.flanking_quals[ 1 ] = <int>( galign.rd.query_qualities[ galign.gali_ugofsRead[gi] ] )
                break
            gi+=1

        for gi in range( self.gofs_start, self.gofs_end+1 ):
            if galign.gali_nongapRdRefMatch[ gi, 0 ] != 0:
                # nongap in read
                if galign.gali_nongapRdRefMatch[ gi, 1 ] != 0:
                    # nongap in ref and read
                    # is it a match?
                    if galign.gali_nongapRdRefMatch[ gi, 2 ] != 0:
                        self.n_match += 1
                        self.gapped_query_seq += (galign.query_seq[ galign.gali_ugofsRead[gi] ]).lower()
                        self.gapped_ref_seq += galign.ref_seq[ galign.gali_ugofsRef[gi] ]
                    else:
                        self.n_mm += 1

                        bpq = galign.query_seq[ galign.gali_ugofsRead[gi] ].upper() 
                        bpr = galign.ref_seq[ galign.gali_ugofsRef[gi] ].upper() 
                        
                        if bpq=='A':
                            if bpr=='G':
                                self.n_mm_ti += 1
                            elif (bpr=='C' or bpr=='T'):
                                self.n_mm_tv += 1
                        elif bpq=='C':
                            if bpr=='T':
                                self.n_mm_ti += 1
                            elif (bpr=='A' or bpr=='G'):
                                self.n_mm_tv += 1
                        elif bpq=='G':
                            if bpr=='A':
                                self.n_mm_ti += 1
                            elif (bpr=='C' or bpr=='T'):
                                self.n_mm_tv += 1
                        elif bpq=='T':
                            if bpr=='C':
                                self.n_mm_ti += 1
                            elif (bpr=='A' or bpr=='C'):
                                self.n_mm_tv += 1

                        self.gapped_query_seq += bpq
                        # self.gapped_query_seq += (galign.query_seq[ galign.gali_ugofsRead[gi] ]).upper()
                        self.gapped_ref_seq += galign.ref_seq[ galign.gali_ugofsRef[gi] ]

                    self.gapped_query_qual[ gi - self.gofs_start ] = <int>( galign.rd.query_qualities[ galign.gali_ugofsRead[gi] ] )
                else:
                    # nongap in query, gap in ref --> ins
                    self.n_ins += 1
                    self.gapped_query_seq += galign.query_seq[ galign.gali_ugofsRead[gi] ]
                    self.gapped_ref_seq += '-'

                    self.gapped_query_qual[ gi - self.gofs_start ] = <int>( galign.rd.query_qualities[ galign.gali_ugofsRead[gi] ] )
            else:
                # gap in read
                if galign.gali_nongapRdRefMatch[ gi, 1 ] != 0:
                    # nongap in ref --> del
                    self.gapped_query_seq += '-'
                    self.gapped_ref_seq += galign.ref_seq[ galign.gali_ugofsRef[gi] ]
                    self.n_del+=1

                    self.gapped_query_qual[ gi - self.gofs_start ] = -1
                else:
                    # this shouldn't happen
                    assert 1==2




cdef class GappedAlign:

    # coordinate of starting position within the reference
    # 0-based
    # if alignment starts w/ insertion then = reference position preceeding the insertion
    cdef readonly int coz_ref_start

    # 0-based
    # if alignment ends w/ insertion then = reference position preceeding insertion
    cdef readonly int coz_ref_end

    # gapped alignment
    #   [:,0] --> nongap in read
    #   [:,1] --> nongap in ref
    #   [:,2] --> match
    cdef readonly int[:,:] gali_nongapRdRefMatch 

    # ungapped sequence offsets, by gapped alignment position
    cdef readonly int[:] gali_ugofsRef
    cdef readonly int[:] gali_ugofsRead

    # length of gapped alignment
    cdef readonly int gali_len


    # cdef const unsigned char[:] query_seq 
    # cdef const unsigned char[:] ref_seq 
    cdef readonly str query_seq 
    cdef readonly str ref_seq 
    cdef AlignedSegment rd


    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __init__( self,
                  AlignedSegment rd,
                  FastaFile reffa ):
        
        self.rd=rd

        self.coz_ref_start = self.rd.reference_start
        self.coz_ref_end = self.rd.reference_end - 1 # pysam defines ref end as 1 past last aligned based

        self.ref_seq = reffa.fetch( self.rd.reference_name, self.rd.reference_start, self.rd.reference_end+1 ).upper() #.encode('ascii')
        self.query_seq = self.rd.query_sequence#.encode('ascii')

        cdef int ici=0, Nci=len(self.rd.cigartuples)
        cdef int ci0
        cdef int cozgPrev=0, cozgNext=0

        cdef int iung_ref=0, iung_query=0, ctr=0

        self.gali_len = 0
        for ici in range(0,Nci): #0 <= ici < Nci:
            self.gali_len += self.rd.cigartuples[ici][1]

        self.gali_nongapRdRefMatch=np.zeros( (self.gali_len, 3), dtype=np.int32 )
        self.gali_ugofsRef=np.zeros( (self.gali_len), dtype=np.int32 )
        self.gali_ugofsRead=np.zeros( (self.gali_len), dtype=np.int32 )

        for ici in range(0,Nci): # for 0 <= ici < Nci:
            ci0 = self.rd.cigartuples[ ici ][ 0 ]

            if ci0!=5:
                cozgPrev = cozgNext
                cozgNext = cozgPrev + self.rd.cigartuples[ ici ][ 1 ]

            # print('prev->next',cozgPrev,cozgNext,ci0)

            if ci0==0:
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 0 ] = 1 
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 1 ] = 1  

                for ctr in range(cozgPrev, cozgNext) : #cozgPrev <= ctr < cozgNext:
                    # print('CIGAR0 ctr,iung_query,iung_ref: ',ctr,iung_query,iung_ref)
                    self.gali_ugofsRef[ctr] = iung_ref
                    self.gali_ugofsRead[ctr] = iung_query
                    self.gali_nongapRdRefMatch[ ctr, 2 ] = <int>(self.query_seq[iung_query]==self.ref_seq[iung_ref])
                    iung_ref+=1
                    iung_query+=1

            elif ci0==1: # insertion to ref
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 0 ] = 1
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 1 ] = 0

                for ctr in range(cozgPrev, cozgNext) :  # cozgPrev <= ctr < cozgNext:
                    # print('CIGAR INS ctr,iung_query,iung_ref: ',ctr,iung_query,iung_ref)
                    self.gali_ugofsRef[ctr] = iung_ref
                    self.gali_ugofsRead[ctr] = iung_query
                    self.gali_nongapRdRefMatch[ ctr, 2 ] = 0
                    iung_query+=1

            elif ci0==2: # deletion from ref
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 0 ] = 0
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 1 ] = 1

                for ctr in range(cozgPrev, cozgNext) : #  cozgPrev <= ctr < cozgNext:
                    # print('CIGAR DEL ctr,iung_query,iung_ref: ',ctr,iung_query,iung_ref)
                    self.gali_ugofsRef[ctr] = iung_ref
                    self.gali_ugofsRead[ctr] = iung_query
                    self.gali_nongapRdRefMatch[ ctr, 2 ] = 0
                    iung_ref+=1

            elif ci0==4: # masked query
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 0 ] = -1
                self.gali_nongapRdRefMatch[ cozgPrev:cozgNext, 1 ] = 0

                for ctr in range(cozgPrev, cozgNext) : # for  cozgPrev <= ctr < cozgNext:
                    # print('CIGAR MASK ctr,iung_query,iung_ref: ',ctr,iung_query,iung_ref)
                    self.gali_ugofsRef[ctr] = iung_ref
                    self.gali_ugofsRead[ctr] = iung_query
                    self.gali_nongapRdRefMatch[ ctr, 2 ] = 0
                    iung_query+=1


    cpdef debug_dump( self ):
        return {'ugofsRef':self.gali_ugofsRef, 
                'ugofsRead':self.gali_ugofsRead,
                'nongapReadRefMm':self.gali_nongapRdRefMatch }


    # if we can't query this interval for some reason, return None
    # else return a list of VariantRecords
    cpdef list find_variants_aggbydist( self, 
                                   int coz_target_start,
                                   int coz_target_end,
                                   int maxdist=0 ):

        cdef list lvars=[]
        cdef tuple gintv = self.ref_to_gapped_interval( coz_target_start, coz_target_end )

        # if not gintv[0]:
            # return None

        cdef int gi = gintv[1]
        cdef int gj
        cdef int gimax = gintv[2]

        # where does current variant start?  (<0 --> not in a variant)
        cdef int gi_start_curvar = -1
        cdef int gi_end_curvar 
        
        while gi <= gimax:

            # not in a variant. does one start here?
            if (self.gali_nongapRdRefMatch[gi,0]!=1) | (self.gali_nongapRdRefMatch[gi,1]!=1) | (self.gali_nongapRdRefMatch[gi,2]!=1):

                gi_start_curvar = gi
                gi_end_curvar = gi

                gj = gi_end_curvar + 1
                while gj <= gi_end_curvar + 1 + maxdist:
                    if gj > gimax:
                        break
                    if (self.gali_nongapRdRefMatch[gj,0]!=1) | (self.gali_nongapRdRefMatch[gj,1]!=1) | (self.gali_nongapRdRefMatch[gj,2]!=1):
                        gi_end_curvar = gj
                
                    gj+=1

                # print('found variant gapped ofs range []':,gi_start_curvar,gi_end_curvar)

                v = Variant( self, gi_start_curvar, gi_end_curvar )

                lvars.append(v)

                gi = gi_end_curvar + 1  

            else:
                gi+=1

        return lvars


    # TODO - could optimize - the exon break calculation will be repeated again and again...

    # assumes only 1 exon overlapping target region
    #
    # strand = 1 or -1
    # if plus, then coz_exon_start <= coz_exon_end
    # if minus then coz_exon_start >= coz_exon_end
    # 
    # frame_start is the frame (0,1,2) of the first base of the exon (on the sense strand)
    cpdef list find_variants_aggbydist_groupbycodon( self, 
                                   int coz_target_start,
                                   int coz_target_end,
                                   int coz_exon_start,
                                   int coz_exon_end,
                                   int strand,
                                   int frame_start,
                                   int maxdist=0 ):

        cdef int nbreaks = 0
        cdef int coz_firstpos = 0
        cdef int frame_firstpos = 0
        # reference coordinates of breaks, i.e., when a contiguous variant intersects one of these at pos i, it should be broken
        # into two ...i-1, and i....
        # cdef int[:] coz_breaks = np.zeros( int((coz_target_end-coz_target_start+1)/3) + 1, dtype=np.int32 )
        #debug fudge
        cdef int[:] coz_breaks = np.zeros( int((coz_target_end-coz_target_start+1)/3) + 1000, dtype=np.int32 )

        cdef int coz_cur 

        cdef bint skipbreaks = 0

        if strand > 0:
            # for + exons, we want to break on frame 2 positions
            #  ------->
            #  12012012
            #  aaBBBccc

            if coz_exon_start < coz_target_start :
                #exon starts before target. 
                if coz_exon_end >= coz_target_start:
                    # if exon ends within or after target, then there is overlap.
                    # start counting frame 0 positions - these are where variants should break (i.e., to start at these positions)
                    coz_firstpos = coz_target_start
                    frame_firstpos = ( frame_start + coz_target_start - coz_exon_start ) % 3
                    coz_breaks[ nbreaks ] = coz_firstpos
                    nbreaks += 1
                else:
                    #don't do anything
                    skipbreaks = 1

            elif coz_exon_start <= coz_target_end:
                #exon starts within target. 
                coz_firstpos = coz_exon_start
                frame_firstpos = frame_start
                coz_breaks[ nbreaks ] = coz_firstpos
                nbreaks += 1

            else:
                skipbreaks = 1
                    

            if ~skipbreaks:
                # advance to next frame-0 position
                coz_cur = coz_firstpos + (3-frame_firstpos)
                while coz_cur <= min(coz_target_end, coz_exon_end):
                    coz_breaks[ nbreaks ] = coz_cur
                    nbreaks += 1
                    coz_cur += 3

                # add a break for the exon end if not already there
                if coz_exon_end <= coz_target_end:
                    if coz_breaks[ nbreaks - 1 ] != coz_exon_end:
                        coz_breaks[ nbreaks ] = coz_exon_end
                        nbreaks += 1

        elif strand < 0 :
            # for - exons, we want to break on frame 2 positions
            #  <-------
            #  21021021  
            #  aaaBBBcc

            if coz_exon_start > coz_target_end:
                # exon ends (starts, on sense strand) beyond target
                if coz_exon_end <= coz_target_end:
                    coz_firstpos = coz_target_end
                    frame_firstpos = ( frame_start + coz_exon_start - coz_target_end ) % 3
                else:
                    skipbreaks = 1

            elif coz_exon_end >= coz_target_start:
                # exon ends within target
                coz_first_pos = coz_exon_start
                frame_firstpos = frame_start

            else:
                skipbreaks = 1

            if ~skipbreaks:
                # advance (left) to next frame-2 position 
                coz_cur = coz_firstpos - ((3-frame_firstpos) - 1)
                while coz_cur >= max(coz_exon_end, coz_target_start):
                    coz_breaks[ nbreaks ] = coz_cur
                    nbreaks += 1
                    coz_cur -= 3

                # add a break for the exon end if not already there
                if coz_exon_end >= coz_target_start:
                    if coz_breaks[ nbreaks - 1 ] != coz_exon_end:
                        coz_breaks[ nbreaks ] = coz_exon_end
                        nbreaks += 1


        # print('ref coord breaks: '+
        #      ','.join( ['%d:%d'%(i,coz_breaks[i]) for i in range(nbreaks)] ))

        cdef list lvars=[]
        cdef tuple gintv = self.ref_to_gapped_interval( coz_target_start, coz_target_end )

        if not gintv[0]:
            return None

        cdef int gi = gintv[1]
        cdef int gj
        cdef int gimax = gintv[2]

        # where does current variant start?  (<0 --> not in a variant)
        cdef int gi_start_curvar = -1
        cdef int gi_end_curvar 
        
        cdef int i_nextbreak
        cdef int coz_ref_nextbreak
        # cdef int coz_cur


        if strand > 0:
            # breaks are sorted from small to large
            i_nextbreak = 0
        elif strand < 0:   
            # breaks are sorted from large to small 
            i_nextbreak = nbreaks - 1

        while gi <= gimax:

            # not in a variant. does one start here?
            if (self.gali_nongapRdRefMatch[gi,0]!=1) | (self.gali_nongapRdRefMatch[gi,1]!=1) | (self.gali_nongapRdRefMatch[gi,2]!=1):

                gi_start_curvar = gi
                gi_end_curvar = gi

                # look for the first break that is AFTER the current position
                if strand>0:
                    while i_nextbreak < nbreaks:
                        coz_ref_nextbreak = coz_breaks[ i_nextbreak ]
                        if coz_ref_nextbreak <= self.gali_ugofsRef[gi]+self.coz_ref_start:
                            i_nextbreak+=1
                        else:
                            break
                else:
                    while i_nextbreak >= 0:
                        coz_ref_nextbreak = coz_breaks[ i_nextbreak ]
                        if coz_ref_nextbreak <= self.gali_ugofsRef[gi]+self.coz_ref_start:
                            i_nextbreak-=1
                        else:
                            break

                gj = gi_end_curvar + 1
                while gj <= gi_end_curvar + 1 + maxdist:
                    if gj > gimax:
                        break

                    # are we at a break?
                    if ((strand > 0) and (i_nextbreak < nbreaks)) or \
                       ((strand < 0) and (i_nextbreak >= 0)):

                        if coz_ref_nextbreak == self.gali_ugofsRef[gj] + self.coz_ref_start:
                            # we've reached a break! stop here.
                            break

                    if (self.gali_nongapRdRefMatch[gj,0]!=1) | (self.gali_nongapRdRefMatch[gj,1]!=1) | (self.gali_nongapRdRefMatch[gj,2]!=1):
                        gi_end_curvar = gj
                
                    gj+=1

                # print('found variant gapped ofs range []':,gi_start_curvar,gi_end_curvar)

                v = Variant( self, gi_start_curvar, gi_end_curvar )

                lvars.append(v)

                gi = gi_end_curvar + 1  

            else:
                gi+=1

        return lvars
        

    # returns tuple
    # ( does alignment cover a given range, without having indels @ edges? 
    #   gapped start ofs
    #   gapped end ofs )
    #  

    # those gapped coords will be -9999999 if the alignment never covers these potisions

    cpdef tuple ref_to_gapped_interval( self,
                             int check_coz_ref_start,
                             int check_coz_ref_end ):

        cdef int gistart=-9999999, giend=-9999999

        cdef bint covers_start = 0, covers_end = 0

        if self.coz_ref_start + self.gali_ugofsRef[ 0 ] > check_coz_ref_start:
            covers_start = 0

        # is requested start covered?
        for gi in range(0,self.gali_len): # 0 <= gi < self.gali_len:
            # print('checkstart',gi,self.coz_ref_start + self.gali_ugofsRef[ gi ])

            # is there an indel at the queried start position?
            if self.coz_ref_start + self.gali_ugofsRef[ gi ] == check_coz_ref_start:
                if ( self.gali_nongapRdRefMatch[ gi, 0 ] != 1 or 
                     self.gali_nongapRdRefMatch[ gi, 1 ] != 1 ):
                    covers_start = 0
                    gistart = gi
                    break
                else:
                    covers_start = 1
                    gistart = gi
                    break


        if self.coz_ref_start + self.gali_ugofsRef[ self.gali_len - 1 ] < check_coz_ref_end:
            covers_end = False
        else:
            # is requested end covered?
            for gi in range(self.gali_len - 1, -1, -1):
                # print('checkend',gi,self.coz_ref_start + self.gali_ugofsRef[ gi ])    

                # is there an indel at the queried end position?
                if self.coz_ref_start + self.gali_ugofsRef[ gi ] == check_coz_ref_end:
                    if ( self.gali_nongapRdRefMatch[ gi, 0 ] != 1 or 
                         self.gali_nongapRdRefMatch[ gi, 1 ] != 1 ):
                        covers_end = 0
                        giend = gi 
                        break
                    else:
                        covers_end = 1
                        giend=gi
                        break

        if covers_start and covers_end:
            return (True, gistart, giend)
        else:
            return (False, gistart, giend)


