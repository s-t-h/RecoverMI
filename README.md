# **RecoverMI**

`Author:    Simon Hackl`
`Contact:   simon.hackl@uni-tuebingen.de`

---
### **Usage**
    python RecoverMI.py [-h] I O
    
    RecoverMI: Fix and recover equally well multi-mapped read information from .sam files.
    
    positional arguments:
      I           The .sam file used as input. Has to be sorted with `samtools sort -n`
      O           The output directory at which the fixed file shall be stored.
    
    optional arguments:
      -h, --help  show this help message and exit

---

### **Description**
`RecoverMI` expects a SAM format file, sorted by the read name (cf. QNAME field in SAM format), as input and will generate a SAM file with multi-mapped alignment information being recovered in return.

To proceed successfully the input file should contain information about all found alignments, i.e. secondary alignments in SAM format, of the respective mapping tool.

While iterating over the file all alignments per read are collected and treated as follows:
- Entries that are marked as being unmapped in the SAM bitwise flag are ignored.
- Those alignments that are marked as mapped non-secondary alignments containing a sequence and per-base quality record as well as a mapping quality greater than 0 are written as-is into the output file.
- In case of the collected alignments containing one primary alignment (PA) and one or more secondary alignments (SAs), first, the quality field of the PA is set to 60, if it yields a value of 0 in the input, and adopted unchanged otherwise. - Next, each SA is validated for mapping equally well as the PA. This is done by comparing the CIGAR string of the SA and PA with respect to the indicated reverse complementation in the bitwise flag (i.e. the CIGAR string of the PA is reversed if the alignments do not share the same orientation).
- If the SA’s CIGAR string contains clipped positions, i.e. unaligned positions that were retained for completeness, that are not clipped in the PA, summing up to more than 5% of the length of the PA, it is rejected. Otherwise it is accepted. 
- For all accepted SAs the secondary alignment bit in the bitwise flag is set to 0 and the mapping quality of the PA is copied.
- The sequence and per-base quality records for the SAs are inferred from the respective PA record, based on the CIGAR string of the SA with respect to the orientation of the alignments. Finally, all collected accepted alignments are written to the output file.

---
### **Performance**
The script was implemented in the course of my master thesis `The OMPeome of Treponema pallidum` (at the University Tübingen) to allow variant calling in repeat regions. Thereby, samples with 10% to 80% of disregarded positions (in terms of variant calling) without `RecoverMI` being applied had below 1% disregarded positions after application.
