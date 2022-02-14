#!/usr/bin/python

"""
    RecoverMI
    Author: Simon Hackl
    Contact: simon.hackl@uni-tuebingen.de

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from itertools import groupby

import argparse
import os
import sys

def complementSeq( seq ) :
    """
    Returns the reverse complement of a nucleotide sequence.
    """
    compSeq = ""
    compDict = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }
    for c in list( seq )[ :: -1 ] :
        compSeq += compDict[ c ]
    return compSeq

def complementCIGAR( CIGAR ) :
    """
    Returns a reversed CIGAR string of the input CIGAR.
    """
    splitCIGAR = [ ''.join( group ) for _, group in groupby( CIGAR, str.isalpha ) ]
    revCIGAR = ""
    times = 0
    instruction = None
    for i in [ len( splitCIGAR ) - j - 1 for j in range( 0, len( splitCIGAR ) ) ] :
        if i % 2 == 0 :
            times = splitCIGAR[ i ]
            revCIGAR += times
            revCIGAR += instruction
        else :
            instruction = splitCIGAR[ i ]
    return revCIGAR
    

if __name__ == '__main__' :
    """
    Definition of parsable arguments.
    """
    parser = argparse.ArgumentParser(
        description = 'RecoverMI: Fix and recover equally well multi-mapped read information from .sam files.'
    )
    parser.add_argument( 'inputSam', metavar='I', type=str, help='The .sam file used as input. Has to be sorted with `samtools sort -n`')
    parser.add_argument( 'outputDir', metavar='O', type=str, help='The output directory at which the fixed file shall be stored.')
    args = parser.parse_args()
    """
    Validate arguments.
    """
    if not os.path.isfile( args.inputSam ) :
        raise FileNotFoundError( "The specified input file " + args.inputSam + " does not exist." )
    elif os.path.isfile( args.inputSam ) and not args.inputSam.endswith( ".sam" ) :
        raise IOError( "The specified input file " + args.inputSam + " seems not to be a .sam file." )
    elif not os.path.isdir( args.outputDir ) :
        raise IOError( "The specified output directory " + args.outputDir + " is no directory." )
    inputSam = args.inputSam
    outDir = args.outputDir
    if not outDir.endswith( "/" ) :
        outDir += "/"
    outputSam = outDir + os.path.basename( inputSam ).replace( ".sam", "" ) + ".fxd.sam"
    outputLog = outDir + os.path.basename( inputSam ).replace( ".sam", "" ) + ".fxd.log"
    """
    Initialize stats variables.
    """
    cnt_unmapped = 0
    cnt_multiMapped_conf = 0
    cnt_multiMapped_fixed = 0
    cnt_multiMapped_rejected = 0
    cnt_uniquelyMapped_conf = 0
    cnt_uniquelyMapped_rejected = 0
    """
    Initialize variables to store read records with common query name.
    """
    currentQName = ""
    currentQNameEntries = [ ]
    """
    Parse/write input/output file.
    """
    with open( inputSam, "r" ) as samInFile :
        with open( outputSam, "w+", newline = "\n" ) as samOutFile :
            with open( outputLog, "w+", newline = "\n" ) as log :

                log.write( "Command: python " + " ".join( sys.argv ) + "\n\n" )
                log.write( "Issues:\n" )

                """
                Definition of helper functions.
                """
                def getBitwiseFlag( qNameEntry ) :
                    """
                    Returns the FLAG field of a SAM format line as list of boolean values.
                    """
                    return [ True if l == "1" else False for l in list( "{0:b}".format( int( qNameEntry[ 1 ] ) ).zfill( 12 ) ) ]

                def scoreCIGARS( cigar1, cigar2 ) :
                    """
                    Scores two cigars based on clipped positions in cigar1 that were matched in cigar2.

                    Returns the number of such positions relative to the length represented by the cigar.

                    Internally cigar1 should represent the CIGAR string of a secondary alignment and cigar2
                    should represent the CIGAR string of the primary alignment.
                    """
                    cigar1Split = [ ''.join( group ) for _, group in groupby( cigar1, str.isalpha ) ]
                    cigar2Split = [ ''.join( group ) for _, group in groupby( cigar2, str.isalpha ) ]
                    cigar1Length = 0
                    cigar2Length = 0
                    cigar1Full = ""
                    cigar2Full = ""
                    cigarScore = 0
                    times = 0
                    instruction = None
                    for i in range( 0, len( cigar1Split ) ) :
                        if i % 2 == 0 :
                            times = int( cigar1Split[ i ] )
                        else :
                            instruction = cigar1Split[ i ]
                            if instruction == "D" :
                                continue
                            cigar1Length += times
                            cigar1Full += instruction * times
                    times = 0
                    instruction = None
                    for i in range( 0, len( cigar2Split ) ) :
                        if i % 2 == 0 :
                            times = int( cigar2Split[ i ] )
                        else :
                            instruction = cigar2Split[ i ]
                            if instruction == "D" :
                                continue
                            cigar2Length += times
                            cigar2Full += instruction * times
                    if cigar1Length != cigar2Length :
                        return -1
                    for c1, c2 in zip( cigar1Full, cigar2Full ) :
                        if c1 == "H" and c2 in [ "M", "I", "D", "=", "X" ] :
                            cigarScore += 1
                    return cigarScore / cigar1Length

                def processEntries( qNameEntries ) :
                    """
                    Processes a list of SAM format lines.
                    """
                    global cnt_multiMapped_conf
                    global cnt_multiMapped_fixed
                    global cnt_multiMapped_rejected
                    global cnt_uniquelyMapped_conf
                    global cnt_uniquelyMapped_rejected
                    # Process information, if at least one entry was collected.
                    if len( qNameEntries ) >= 1 :
                        # CASE 1: Only one entry was processed, write if mapq is above 0.
                        if len( qNameEntries ) == 1 and int( qNameEntries[ 0 ][ 4 ] ) > 0 :
                            samOutFile.write( "\t".join( qNameEntries[ 0 ] ).strip( ) + "\n" )
                            cnt_uniquelyMapped_conf += 1
                        elif len( qNameEntries ) == 1 and int( qNameEntries[ 0 ][ 4 ] ) == 0 :
                            cnt_uniquelyMapped_rejected += 1
                            return
                        # CASE 2: Read was mapped different ways/multiple times, additional processing is applied.
                        elif len( qNameEntries ) > 1 :
                            # Search for primary entry.
                            primaryEntry = None
                            primaryFlag = None
                            for qNameEntry in qNameEntries :
                                qNameEntryFlag = getBitwiseFlag( qNameEntry )
                                if not qNameEntryFlag[ -9 ] and not qNameEntryFlag[ -12 ] :
                                    if int( qNameEntry[ 4 ] ) == 0 :
                                        qNameEntry[ 4 ] = str( 60 )
                                    primaryEntry = qNameEntry
                                    primaryFlag = qNameEntryFlag
                                    break
                            if primaryEntry == None :
                                log.write( "+\n> No primary entry exists for read " + qNameEntries[ 0 ][ 0 ] + ".\n" )
                                log.write( "> Entries:\n" )
                                for qNameEntry in qNameEntries :
                                    log.write( ">\t" + "\t".join( qNameEntry ) )
                                return
                            # Process each qNameEntry dependent on its content relative to the primary entry.
                            for qNameEntry in qNameEntries :
                                # If the entry is the primary entry it is accepted.
                                if qNameEntry == primaryEntry :
                                    samOutFile.write( "\t".join( primaryEntry ).strip( ) + "\n" )
                                    cnt_multiMapped_conf += 1
                                # Else, if the entry has a SEQ record itself it is accepted.
                                elif qNameEntry[ 9 ] != "*" :
                                    # Remove possible secondary alignment flag.
                                    qNameEntryFlag = getBitwiseFlag( qNameEntry )
                                    qNameEntryFlag[ -9 ] = False
                                    qNameEntry[ 1 ] = str( sum( [ 2**i if qNameEntryFlag[ :: -1 ][ i ] else 0 for i in range( 0, len( qNameEntryFlag ) ) ] ) )
                                    # Adjust mapping quality.
                                    if int( qNameEntry[ 4 ] ) == 0 :
                                        qNameEntry[ 4 ] = str( 60 )
                                    samOutFile.write( "\t".join( qNameEntry ).strip( ) + "\n" )
                                    cnt_multiMapped_conf += 1
                                # Else, the entry has to be compared to the primary entry.
                                elif qNameEntry != primaryEntry :
                                    # Remove possible secondary alignment flag.
                                    qNameEntryFlag = getBitwiseFlag( qNameEntry )
                                    qNameEntryFlag[ -9 ] = False
                                    # Initzialize a fixed entry.
                                    fxdEntry = primaryEntry.copy( )[ : 11 ] + qNameEntry[ 11 : ]
                                    # CASE: Current and primary entry are mapped equally.
                                    # Reversed informaion -> entryFlag[ -5 ]
                                    # CIGAR information -> entry[ 5 ]
                                    if primaryFlag[ -5 ] == qNameEntryFlag[ -5 ] and primaryEntry[ 5 ] == qNameEntry[ 5 ] :
                                        fxdEntry[ 3 ] = qNameEntry[ 3 ] 
                                        samOutFile.write( "\t".join( fxdEntry ).strip( ) + "\n" )
                                        cnt_multiMapped_fixed += 1
                                        continue
                                    # If current and primary entries have different orientation, records have to be adjusted.
                                    if primaryFlag[ -5 ] != qNameEntryFlag[ -5 ] :
                                        fxdEntrySeq = complementSeq( primaryEntry[ 9 ] )
                                        fxdEntryBaseQual = "".join( list( primaryEntry[ 10 ] )[ :: -1 ] )
                                        primaryEntryCIGAR = complementCIGAR( primaryEntry[ 5 ] )
                                    else :
                                        fxdEntrySeq = primaryEntry[ 9 ]
                                        fxdEntryBaseQual = primaryEntry[ 10 ]
                                        primaryEntryCIGAR = primaryEntry[ 5 ]
                                    # If current and primary entry have different CIGAR strings, sequences may be adjusted or the entry is rejected.
                                    cigarsScore = scoreCIGARS( qNameEntry[ 5 ], primaryEntryCIGAR )
                                    if cigarsScore == -1 :
                                        log.write( "+\n> CIGARs indicate different sequence length, unable to reconstruct sequence for entry of read " + qNameEntries[ 0 ][ 0 ] + ".\n" )
                                        log.write( "> Primary entry:\t" + "\t".join( primaryEntry ) )
                                        log.write( "> Affected entry:\t" + "\t".join( qNameEntry ) )
                                        continue
                                    elif cigarsScore < 0.05 :
                                        qNameEntrySplitCIGAR = [ ''.join( group ) for _, group in groupby( qNameEntry[ 5 ], str.isalpha ) ]
                                        fxdEntrySeqList = list( fxdEntrySeq )
                                        fxdEntrySeq = ""
                                        fxdEntryBaseQualList = list( fxdEntryBaseQual )
                                        fxdEntryBaseQual = ""
                                        times = 0
                                        idx = 0
                                        try :
                                            for i in range( 0, len( qNameEntrySplitCIGAR ) ) :
                                                if i % 2 == 0 :
                                                    times = int( qNameEntrySplitCIGAR[ i ] )
                                                else :
                                                    if qNameEntrySplitCIGAR[ i ] == "H" :
                                                        idx += times
                                                    elif qNameEntrySplitCIGAR[ i ] == "D" :
                                                        continue
                                                    elif qNameEntrySplitCIGAR[ i ] in [ "M", "I", "S", "=", "X" ] :
                                                        for j in range( 0, times ) :
                                                            fxdEntrySeq += fxdEntrySeqList[ idx ]
                                                            fxdEntryBaseQual += fxdEntryBaseQualList[ idx ]
                                                            idx += 1
                                        except IndexError :
                                            log.write( "+\n> Processing failed for entry of read " + qNameEntries[ 0 ][ 0 ] + ".\n" )
                                            log.write( "> Primary entry:\t" + "\t".join( primaryEntry ) )
                                            log.write( "> Affected entry:\t" + "\t".join( qNameEntry ) )
                                            continue
                                        fxdEntry[ 1 ] = str( sum( [ 2**i if qNameEntryFlag[ :: -1 ][ i ] else 0 for i in range( 0, len( qNameEntryFlag ) ) ] ) )
                                        fxdEntry[ 3 ] = qNameEntry[ 3 ]  
                                        fxdEntry[ 5 ] = qNameEntry[ 5 ]
                                        fxdEntry[ 9 ] = fxdEntrySeq
                                        fxdEntry[ 10 ] = fxdEntryBaseQual
                                        samOutFile.write( "\t".join( fxdEntry ).strip( ) + "\n" )
                                        cnt_multiMapped_fixed += 1
                                    else :
                                        cnt_multiMapped_rejected += 1

                line = samInFile.readline( )
                while line :
                    if line.startswith( "@" ) :
                        samOutFile.write( line )
                        line = samInFile.readline( )
                        continue
                    # Field values: qName, flag, rName, pos, mapq, cigar, rnext, pnext, tlen, seq, qual
                    fields = line.split( "\t" )
                    bitwiseFlag = getBitwiseFlag( fields )
                    if bitwiseFlag[ -3 ] : # Skip unmapped reads.
                        line = samInFile.readline( )
                        cnt_unmapped += 1
                        continue
                    # Collect all reads with equal query name.
                    if fields[ 0 ] == currentQName :
                        currentQNameEntries.append( fields )
                    else :
                        processEntries( currentQNameEntries )
                        # Start collecting information about next read.
                        currentQName = fields[ 0 ]
                        currentQNameEntries = [ fields ]
                    line = samInFile.readline( )
                processEntries( currentQNameEntries )

                log.write( "\nStatistics:\n" )
                log.write( "Removed reads, unmapped:\t\t\t\t\t\t" + str( cnt_unmapped ) + "\n" )
                log.write( "Removed reads, uniquely mapped with low quality:\t\t\t" + str( cnt_uniquelyMapped_rejected ) + "\n" )
                log.write( "Adopted reads, uniquely mapped with high quality:\t\t\t" + str( cnt_uniquelyMapped_conf ) + "\n" )
                log.write( "Removed multi-map entries, failed similarity check to primary entry:\t" + str( cnt_multiMapped_rejected ) + "\n" )
                log.write( "Adopted multi-map entries, mapped with high quality:\t\t\t" + str( cnt_multiMapped_conf ) + "\n" )
                log.write( "Fixed multi-map entries, passed similarity check to primary entry:\t" + str( cnt_multiMapped_fixed ) + "\n" )
