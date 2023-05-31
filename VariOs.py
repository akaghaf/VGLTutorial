 
from collections import OrderedDict, defaultdict
from datetime import datetime
from tkinter.filedialog import askopenfilename

import matplotlib
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib_venn import venn3, venn3_circles

matplotlib.use('Agg')
import concurrent.futures

import venn


def fbgbfileparsing(filename):
    print("Reading FreeBayes Genbank file...")
    with open("samples.txt") as fp:
        fp = fp.read()
    sampleslist = fp.split()
    samplessublists = []
    def divide_chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]
    samplessublists = list(divide_chunks(sampleslist, 3))
    cases = []
    for item in samplessublists:
        if item[-1].lower() == "c":
            cases.append(item[0])
    filefbgb = filename
    with open(filefbgb) as f1a:
        f1a = f1a.read()
        f1a = f1a.split("\n")
        f1a = (OrderedDict.fromkeys(f1a))
        f1a = list(f1a)
        while ("" in f1a):
            f1a.remove("")
    f1 = []
    fbgbheader = []
    for item in f1a:
        if item.startswith("#CHROM"):
            fbgbheader.append(item)
        if item[0] != "#":
            f1.append(item)
    columnheader = fbgbheader[0].split("\t")
    columnheader = columnheader[9:]
    caseindices = []
    for i in range(0,len(columnheader)):
        if columnheader[i] in cases:
            caseindices.append(i)

    fbgbvariants = []
    fbgballeles = []
    fbgbannotation = []
    fbgballelefrequnfiltered = []
    for item in f1:
        item = item.split("\t")
        fbgbvariants.append(item[2])
        allele = str(item[3] +  " -> " + item[4])
        fbgballeles.append(allele)

        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        fbgbannotation.append(annotation)
        fbgballelefrequnfiltered.append(item[9:])
    fbgballelefreqpreprocessedsublist = []
    fbgballelefreqpreprocessedwholelist = []
    for item in fbgballelefrequnfiltered:
        fbgballelefreqpreprocessedsublist = []
        for x in item:
            fbgballelefreqpreprocessedsublist.append(x[:x.index(":")])
        fbgballelefreqpreprocessedwholelist.append(fbgballelefreqpreprocessedsublist)
    fbgbfreqslist = []
    for item in fbgballelefreqpreprocessedwholelist:
        substr = ""
        for elem in item:
            if len(elem) == 1: substr += "./."
            else: substr += elem
# Here is where we subset cases versus controls!
        substrsub = ""
        for i in caseindices:
            substrsub += substr[i*3:i*3+3]
# Above^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        substr = substr.replace("/","").replace(".","")
        substrsub = substrsub.replace("/","").replace(".","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substr)))
        sublist = []
        for item in substrchars:
            sublist += [str(substrsub.count(item)/len(substrsub))]
        fbgbfreqslist.append(sublist)
    fbgbfreqslistmaxval = max(map(len, fbgbfreqslist))
    fbgb0 = []
    fbgb1 = []
    fbgb2 = []
    for item in fbgbfreqslist:
        if len(item) < fbgbfreqslistmaxval:
            item.extend('' for _ in range(fbgbfreqslistmaxval - len(item)))
        fbgb0.append(item[0])
        fbgb1.append(item[1])
        if len(item) == 3:
            fbgb2.append(item[2])
    return fbgb0, fbgb1, fbgb2, fbgballeles, fbgbvariants, fbgbannotation, fbgbheader, "fbgb"

def fbrsfileparsing(filename):
    print("Reading FreeBayes RefSeq file...")
    filefbrs = filename
    with open(filefbrs) as f2a:
        f2a = f2a.read()
        f2a = f2a.split("\n")
        f2a = (OrderedDict.fromkeys(f2a))
        f2a = list(f2a)
        while ("" in f2a):
            f2a.remove("")
    f2 = []
    fbrsheader = []
    for item in f2a:
        if item.startswith("#CHROM"):
            fbrsheader.append(item)
        if item[0] != "#":
            f2.append(item)
    fbrsvariants = []
    fbrsalleles = []
    fbrsannotation = []
    fbrsallelefrequnfiltered = []
    for item in f2:
        item = item.split("\t")
        fbrsvariants.append(item[2])
        allele = str(item[3] +  " -> " + item[4])
        fbrsalleles.append(allele)
        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        fbrsannotation.append(annotation)
        fbrsallelefrequnfiltered.append(item[9:])
    fbrsallelefreqpreprocessedsublist = []
    fbrsallelefreqpreprocessedwholelist = []
    for item in fbrsallelefrequnfiltered:
        fbrsallelefreqpreprocessedsublist = []
        for x in item:
            fbrsallelefreqpreprocessedsublist.append(x[:x.index(":")])
        fbrsallelefreqpreprocessedwholelist.append(fbrsallelefreqpreprocessedsublist)
    fbrsfreqslist = []
    for item in fbrsallelefreqpreprocessedwholelist:
        substr = ""
        for elem in item:
            substr += elem
        substr = substr.replace("/","").replace(".","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substr)))
        sublist = []
        for item in substrchars:

            sublist += [str(substr.count(item)/len(substr))]
        fbrsfreqslist.append(sublist)

    fbrsfreqslistmaxval = max(map(len, fbrsfreqslist))
    fbrs0 = []
    fbrs1 = []
    fbrs2 = []
    for item in fbrsfreqslist:
        if len(item) < fbrsfreqslistmaxval:
            item.extend('' for _ in range(fbrsfreqslistmaxval - len(item)))
        fbrs0.append(item[0])
        fbrs1.append(item[1])
        if len(item) == 3:
            fbrs2.append(item[2])
    return fbrs0, fbrs1, fbrs2, fbrsalleles, fbrsvariants, fbrsannotation, fbrsheader, "fbrs"

def gagbfileparsing(filename):
    
    print("Reading GATK GenBank file...")
    with open("samples.txt") as fp:
        fp = fp.read()
    sampleslist = fp.split()
    samplessublists = []
    def divide_chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]
    samplessublists = list(divide_chunks(sampleslist, 3))
    cases = []
    for item in samplessublists:
        if item[-1].lower() == "c":
            cases.append(item[0])
    filegagb = filename
    with open(filegagb) as f3a:
        f3a = f3a.read()
        f3a = f3a.split("\n")
        f3a = (OrderedDict.fromkeys(f3a))
        f3a = list(f3a)
        while ("" in f3a):
            f3a.remove("")
    f3, gagbheader = [], []
    for item in f3a:
        if item.startswith("#CHROM"):
            gagbheader.append(item)
        if item[0] != "#":
            f3.append(item)
    columnheader = gagbheader[0].split("\t")
    columnheader = columnheader[9:]
    caseindices = []
    for i in range(0,len(columnheader)):
        if columnheader[i] in cases:
            caseindices.append(i)

    gagbvariants, gagballeles, gagbannotation, gagballelefrequnfiltered = [],[],[],[]
    for item in f3:
        item = item.split("\t")
        variant = item[0] + "_" + item[1]
        gagbvariants.append(variant)
        allele = str(item[3] +  " -> " + item[4])
        gagballeles.append(allele)
        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        gagbannotation.append(annotation)
        gagballelefrequnfiltered.append(item[9:])
    gagballelefreqpreprocessedsublist = []
    gagballelefreqpreprocessedwholelist = []
    for item in gagballelefrequnfiltered:
        gagballelefreqpreprocessedsublist = []
        for x in item:
            gagballelefreqpreprocessedsublist.append(x[:x.index(":")])
        gagballelefreqpreprocessedwholelist.append(gagballelefreqpreprocessedsublist)
    gagbfreqslist = []
    for item in gagballelefreqpreprocessedwholelist:
        # if "./." in item:
        #     item.remove("./.")
        substr = ""
        for elem in item:
            if len(elem) == 1: substr += "./."
            else: substr += elem
        substrsub = ""
        for i in caseindices:
            substrsub += substr[i*3:i*3+3]
        substrsub = substrsub.replace("/","").replace(".","").replace("|","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substrsub)))
        sublist = []
        for item in substrchars:
            sublist += [str(substrsub.count(item)/len(substrsub))]
        gagbfreqslist.append(sublist)
    gagbfreqslistmaxval = max(map(len, gagbfreqslist))
    gagb0 = []
    gagb1 = []
    gagb2 = []
    for item in gagbfreqslist:
        if len(item) < gagbfreqslistmaxval:
            item.extend('' for _ in range(gagbfreqslistmaxval - len(item)))
        gagb0.append(item[0])
        gagb1.append(item[1])
        if len(item) == 3:
            gagb2.append(item[2])
    return gagb0, gagb1, gagb2, gagballeles, gagbvariants, gagbannotation, gagbheader, "gagb"

def garsfileparsing(filename): 
    print("Reading GATK RefSeq file...")
    filegars = filename
    with open(filegars) as f4a:
        f4a = f4a.read()
        f4a = f4a.replace("NC_009144.3", "chr1_").replace("NC_009145.3", "chr2_").replace("NC_009146.3", "chr3_")
        f4a = f4a.replace("NC_009147.3", "chr4_").replace("NC_009148.3", "chr5_").replace("NC_009149.3", "chr6_")
        f4a = f4a.replace("NC_009150.3", "chr7_").replace("NC_009151.3", "chr8_").replace("NC_009152.3", "chr9_")
        f4a = f4a.replace("NC_009153.3", "chr10_").replace("NC_009154.3", "chr11_").replace("NC_009155.3", "chr12_")
        f4a = f4a.replace("NC_009156.3", "chr13_").replace("NC_009157.3", "chr14_").replace("NC_009158.3", "chr15_")
        f4a = f4a.replace("NC_009159.3", "chr16_").replace("NC_009160.3", "chr17_").replace("NC_009161.3", "chr18_")
        f4a = f4a.replace("NC_009162.3", "chr19_").replace("NC_009163.3", "chr20_").replace("NC_009164.3", "chr21_")
        f4a = f4a.replace("NC_009165.3", "chr22_").replace("NC_009166.3", "chr23_").replace("NC_009167.3", "chr24_")
        f4a = f4a.replace("NC_009168.3", "chr25_").replace("NC_009169.3", "chr26_").replace("NC_009170.3", "chr27_")
        f4a = f4a.replace("NC_009171.3", "chr28_").replace("NC_009172.3", "chr29_").replace("NC_009173.3", "chr30_")
        f4a = f4a.replace("NC_009174.3", "chr31_").replace("NC_009175.3", "chrX_")
        f4a = f4a.split("\n")
        f4a = OrderedDict.fromkeys(f4a)
        f4a = list(f4a)
        while ("" in f4a):
            f4a.remove("")    
    f4 , garsheader = [], []
    for item in f4a:
        if item.startswith("#CHROM"):
            garsheader.append(item)
        if item[0] != "#":
            f4.append(item)
    garsvariants = []
    garsalleles = []
    garsannotation = []
    garsallelefrequnfiltered = []
    for item in f4:
        item = item.split("\t")
        variant = item[0] + item[1]
        garsvariants.append(variant)
        allele = str(item[3] +  " -> " + item[4])
        garsalleles.append(allele)
        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        garsannotation.append(annotation)
        garsallelefrequnfiltered.append(item[9:])
    garsallelefreqpreprocessedsublist = []
    garsallelefreqpreprocessedwholelist = []
    for item in garsallelefrequnfiltered:
        garsallelefreqpreprocessedsublist = []
        for x in item:
            garsallelefreqpreprocessedsublist.append(x[:x.index(":")])
        garsallelefreqpreprocessedwholelist.append(garsallelefreqpreprocessedsublist)
    garsfreqslist = []
    for item in garsallelefreqpreprocessedwholelist:
        if "./." in item:
            item.remove("./.")
        substr = ""
        for elem in item:
            substr += elem
        substr = substr.replace("/","").replace(".","").replace("|","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substr)))
        sublist = []
        for item in substrchars:

            sublist += [str(substr.count(item)/len(substr))]
        garsfreqslist.append(sublist)
    garsfreqslistmaxval = max(map(len, garsfreqslist))
    gars0 = []
    gars1 = []
    gars2 = []
    for item in garsfreqslist:
        if len(item) < garsfreqslistmaxval:
            item.extend('' for _ in range(garsfreqslistmaxval - len(item)))
        gars0.append(item[0])
        gars1.append(item[1])
        if len(item) == 3:
            gars2.append(item[2])
    return gars0, gars1, gars2, garsalleles, garsvariants, garsannotation, garsheader, "gars"

def stgbfileparsing(filename):
    print("Reading Samtools GenBank file...")
    with open("samples.txt") as fp:
        fp = fp.read()
    sampleslist = fp.split()
    samplessublists = []
    def divide_chunks(l, n):
        for i in range(0, len(l), n):
            yield l[i:i + n]
    samplessublists = list(divide_chunks(sampleslist, 3))
    cases = []
    for item in samplessublists:
        if item[-1].lower() == "c":
            cases.append(item[0])
    filestgb = "samtools_genbank/STGB.vcf"
    with open(filestgb) as f5a:
        f5a = f5a.read()
        f5a = f5a.split("\n")
        f5a = (OrderedDict.fromkeys(f5a))
        f5a = list(f5a)
        while ("" in f5a):
            f5a.remove("")
    f5, stgbheader = [], []
    for item in f5a:
        if item.startswith("#CHROM"):
            stgbheader.append(item)
        if item[0] != "#":
            f5.append(item)
    columnheader = stgbheader[0].split("\t")
    columnheader = columnheader[9:]
    caseindices = []
    for i in range(0,len(columnheader)):
        if columnheader[i] in cases:
            caseindices.append(i)

    stgbvariants = []
    stgballeles = []
    stgbannotation = []
    stgballelefrequnfiltered = []
    for item in f5:
        item = item.split("\t")
        stgbvariants.append(item[2])
        allele = str(item[3] +  " -> " + item[4])
        stgballeles.append(allele)
        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        stgbannotation.append(annotation)
        stgballelefrequnfiltered.append(item[9:])
    stgballelefreqpreprocessedsublist = []
    stgballelefreqpreprocessedwholelist = []
    for item in stgballelefrequnfiltered:
        stgballelefreqpreprocessedsublist = []
        for x in item:
            stgballelefreqpreprocessedsublist.append(x[:x.index(":")])
        stgballelefreqpreprocessedwholelist.append(stgballelefreqpreprocessedsublist)
    stgbfreqslist = []
    for item in stgballelefreqpreprocessedwholelist:
        substr = ""
        for elem in item:
            substr += elem
        substrsub = ""
        for i in caseindices:
            substrsub += substr[i*3:i*3+3]
        substr = substr.replace("/","").replace(".","").replace("|","")
        substrsub = substrsub.replace("/","").replace(".","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substr)))
        sublist = []
        for item in substrchars:
            sublist += [str(substrsub.count(item)/len(substrsub))]
        stgbfreqslist.append(sublist)
    stgbfreqslistmaxval = max(map(len, stgbfreqslist))
    stgb0 = []
    stgb1 = []
    stgb2 = []
    for item in stgbfreqslist:
        if len(item) < stgbfreqslistmaxval:
            item.extend('' for _ in range(stgbfreqslistmaxval - len(item)))
        stgb0.append(item[0])
        stgb1.append(item[1])
        if len(item) == 3:
            stgb2.append(item[2])
    return stgb0, stgb1, stgb2, stgballeles, stgbvariants, stgbannotation, stgbheader, "stgb"

def strsfileparsing(filename):
    print("Reading Samtools Refseq file...")
    filestrs = "samtools_refseq/STRS.vcf"
    with open(filestrs) as f6a:
        f6a = f6a.read()
        f6a = f6a.split("\n")
        f6a = (OrderedDict.fromkeys(f6a))
        f6a = list(f6a)
        while ("" in f6a):
            f6a.remove("")
    f6, strsheader = [], []
    for item in f6a:
        if item.startswith("#CHROM"):
            strsheader.append(item)
        if item[0] != "#":
            f6.append(item)
    strsvariants = []
    strsalleles = []
    strsannotation = []
    strsallelefrequnfiltered = []
    for item in f6:
        item = item.split("\t")
        strsvariants.append(item[2])
        allele = str(item[3] +  " -> " + item[4])
        strsalleles.append(allele)
        annotationstart = item[7].index("ANN=")
        annotation = item[7][annotationstart:]
        strsannotation.append(annotation)
        strsallelefrequnfiltered.append(item[9:])

    strsallelefreqpreprocessedsublist = []
    strsallelefreqpreprocessedwholelist = []
    for item in strsallelefrequnfiltered:
        strsallelefreqpreprocessedsublist = []
        for x in item:
            strsallelefreqpreprocessedsublist.append(x[:x.index(":")])
        strsallelefreqpreprocessedwholelist.append(strsallelefreqpreprocessedsublist)
    strsfreqslist = []
    for item in strsallelefreqpreprocessedwholelist:
        if "./." in item:
            item.remove("./.")
        substr = ""
        for elem in item:
            substr += elem
        substr = substr.replace("/","").replace(".","").replace("|","")
        substrchars = "".join(sorted(OrderedDict.fromkeys(substr)))
        sublist = []
        for item in substrchars:

            sublist += [str(substr.count(item)/len(substr))]
        strsfreqslist.append(sublist)

    strsfreqslistmaxval = max(map(len, strsfreqslist))
    strs0 = []
    strs1 = []
    strs2 = []
    for item in strsfreqslist:
        if len(item) < strsfreqslistmaxval:
            item.extend('' for _ in range(strsfreqslistmaxval - len(item)))
        strs0.append(item[0])
        strs1.append(item[1])
        if len(item) == 3:
            strs2.append(item[2])
    return strs0, strs1, strs2, strsalleles, strsvariants, strsannotation, strsheader, "strs"



if __name__  == '__main__':
        
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = [executor.submit(f, filen) for f, filen in zip(
             [fbgbfileparsing,
              fbrsfileparsing,
              gagbfileparsing, 
              garsfileparsing, 
              stgbfileparsing, 
              strsfileparsing],
            ["freebayes_genbank/FBGB.vcf",
             "freebayes_refseq/FBRS.vcf",
             "gatk_genbank/GAGB.vcf",
             "gatk_refseq/GARS.vcf",
             "samtools_genbank/STGB.vcf",
             "samtools_refseq/STRS.vcf"])]
        for f in concurrent.futures.as_completed(results):
            temp = f.result()
            if temp[-1] == "strs":
                strs0, strs1, strs2, strsalleles, strsvariants, strsannotation, strsheader = temp[:-1]
            if temp[-1] == "stgb":
                stgb0, stgb1, stgb2, stgballeles, stgbvariants, stgbannotation, stgbheader = temp[:-1]
            if temp[-1] == "gars":
                gars0, gars1, gars2, garsalleles, garsvariants, garsannotation, garsheader = temp[:-1]
            if temp[-1] == "gagb":
                gagb0, gagb1, gagb2, gagballeles, gagbvariants, gagbannotation, gagbheader = temp[:-1]
            if temp[-1] == "fbrs":
                fbrs0, fbrs1, fbrs2, fbrsalleles, fbrsvariants, fbrsannotation, fbrsheader = temp[:-1]
            if temp[-1] == "fbgb":
                fbgb0, fbgb1, fbgb2, fbgballeles, fbgbvariants, fbgbannotation, fbgbheader = temp[:-1]




    print("Combining variants by caller and removing duplicates...")
    FreeBayesVariants = OrderedDict.fromkeys(fbgbvariants + fbrsvariants)
    GATKVariants = OrderedDict.fromkeys(gagbvariants + garsvariants)
    SamtoolsVariants = OrderedDict.fromkeys(strsvariants + stgbvariants)
    FreeBayesVariants = list(FreeBayesVariants)
    GATKVariants = list(GATKVariants)
    SamtoolsVariants = list(SamtoolsVariants)
    print("Finding and sorting common variants between files...")
    data = [FreeBayesVariants,SamtoolsVariants,GATKVariants]

    common = set(data[0])
    for line in data:
        common = common.intersection(set(line))

    res = defaultdict(int)
    for line in data:
        for idx, item in enumerate(line):
            if item in common:
                res[item] += idx
    CommonToAllCallersVariants = [item[0] for item in sorted(res.items(), key=lambda x: x[1])]

    setFreeBayesVariants = set(FreeBayesVariants)
    setGATKVariants = set(GATKVariants)
    setSamtoolsVariants = set(SamtoolsVariants)

    FreeBayesGATK = [x for x in FreeBayesVariants if x in setGATKVariants]
    FreeBayesSamtools = [x for x in FreeBayesVariants if x in setSamtoolsVariants]
    GATKSamtools = [x for x in GATKVariants if x in setSamtoolsVariants]

    FreeBayesGATKVariants = []
    for item in FreeBayesGATK:
        if item not in CommonToAllCallersVariants:
            FreeBayesGATKVariants.append(item)
    FreeBayesSamtoolsVariants = []
    for item in FreeBayesSamtools:
        if item not in CommonToAllCallersVariants:
            FreeBayesSamtoolsVariants.append(item)
    GATKSamtoolsVariants = []
    for item in GATKSamtools:
        if item not in CommonToAllCallersVariants:
            GATKSamtoolsVariants.append(item)
    FreeBayesUniqueVariants = []
    for item in FreeBayesVariants:
        if item not in GATKVariants and item not in SamtoolsVariants:
            FreeBayesUniqueVariants.append(item)
    GATKUniqueVariants = []
    for item in GATKVariants:
        if item not in FreeBayesVariants and item not in SamtoolsVariants:
            GATKUniqueVariants.append(item)
    SamtoolsUniqueVariants = []
    for item in SamtoolsVariants:
        if item not in GATKVariants and item not in FreeBayesVariants:
            SamtoolsUniqueVariants.append(item)

    HighVariants, ModerateVariants, LowVariants, ModifierVariants = [] , [] , [] , []

    
    genename = []
    effect = []
    Allvariants = CommonToAllCallersVariants + FreeBayesGATKVariants + FreeBayesSamtoolsVariants + GATKSamtoolsVariants + FreeBayesUniqueVariants + SamtoolsUniqueVariants + GATKUniqueVariants
    FreeBayesGenBankAnnotation, FreeBayesGenBankAlleles, FreeBayesRefSeqAnnotation, FreeBayesRefSeqAlleles, GATKGenBankAnnotation, GATKGenBankAlleles = [],[],[],[],[],[]
    GATKRefSeqAnnotation, GATKRefSeqAlleles, SamtoolsGenBankAnnotation, SamtoolsGenBankAlleles, SamtoolsRefSeqAnnotation, SamtoolsRefseqAlleles, CallerCount = [],[],[],[],[],[],[]
    FreebayesRefAlleleFreq, FreebayesAlt1Freq, FreebayesAlt2Freq, GATKRefAlleleFreq, GATKAlt1Freq, GATKAlt2Freq, SamtoolsRefAlleleFreq, SamtoolsAlt1Freq, SamtoolsAlt2Freq = [],[],[],[],[],[],[],[],[]
    print("Extracting all annotations and alleles...")
    for item in Allvariants:
        callercountstring = ""
        if item in fbgbvariants:
            ItemIndex = fbgbvariants.index(item)
            FreeBayesGenBankAnnotation.append(fbgbannotation[ItemIndex])
            FreeBayesGenBankAlleles.append(fbgballeles[ItemIndex])
            FreebayesRefAlleleFreq.append(fbgb0[ItemIndex])
            FreebayesAlt1Freq.append(fbgb1[ItemIndex])
            if len(fbgb2) > 1:
                FreebayesAlt2Freq.append(fbgb2[ItemIndex])
            else:
                FreebayesAlt2Freq.append("N/A")
            callercountstring += "FreeBayes Genbank;"
        else:
            FreeBayesGenBankAnnotation.append("N/A")
            FreeBayesGenBankAlleles.append("N/A")
            FreebayesAlt2Freq.append("N/A")
            FreebayesAlt1Freq.append("N/A")
            FreebayesRefAlleleFreq.append("N/A")
        if item in fbrsvariants:
            ItemIndex = fbrsvariants.index(item)
            FreeBayesRefSeqAnnotation.append(fbrsannotation[ItemIndex])
            FreeBayesRefSeqAlleles.append(fbrsalleles[ItemIndex])
            callercountstring += "FreeBayes Ref Seq;"
        else:
            FreeBayesRefSeqAnnotation.append("N/A")
            FreeBayesRefSeqAlleles.append("N/A")
        if item in gagbvariants:
            ItemIndex = gagbvariants.index(item)
            GATKGenBankAnnotation.append(gagbannotation[ItemIndex])
            GATKGenBankAlleles.append(gagballeles[ItemIndex])
            GATKRefAlleleFreq.append(gagb0[ItemIndex])
            GATKAlt1Freq.append(gagb1[ItemIndex])
            if len(gagb2) > 1:
                GATKAlt2Freq.append(gagb2[ItemIndex])
            else:
                GATKAlt2Freq.append("N/A")
            callercountstring += "GATK Genbank;"
        else:
            GATKGenBankAnnotation.append("N/A")
            GATKGenBankAlleles.append("N/A")
            GATKAlt2Freq.append("N/A")
            GATKAlt1Freq.append("N/A")
            GATKRefAlleleFreq.append("N/A")
        if item in garsvariants:
            ItemIndex = garsvariants.index(item)
            GATKRefSeqAnnotation.append(garsannotation[ItemIndex])
            GATKRefSeqAlleles.append(garsalleles[ItemIndex])
            callercountstring += "GATK Ref Seq;"
        else:
            GATKRefSeqAnnotation.append("N/A")
            GATKRefSeqAlleles.append("N/A")
        if item in stgbvariants:
            ItemIndex = stgbvariants.index(item)
            SamtoolsGenBankAnnotation.append(stgbannotation[ItemIndex])
            SamtoolsGenBankAlleles.append(stgballeles[ItemIndex])
            SamtoolsRefAlleleFreq.append(stgb0[ItemIndex])
            SamtoolsAlt1Freq.append(stgb1[ItemIndex])
            if len(stgb2) > 1:
                SamtoolsAlt2Freq.append(stgb2[ItemIndex])
            else:
                SamtoolsAlt2Freq.append("N/A")
            callercountstring += "Samtools Genbank;"
        else:
            SamtoolsGenBankAnnotation.append("N/A")
            SamtoolsGenBankAlleles.append("N/A")
            SamtoolsAlt1Freq.append("N/A")
            SamtoolsAlt2Freq.append("N/A")
            SamtoolsRefAlleleFreq.append("N/A")
        if item in strsvariants:
            ItemIndex = strsvariants.index(item)
            SamtoolsRefSeqAnnotation.append(strsannotation[ItemIndex])
            SamtoolsRefseqAlleles.append(strsalleles[ItemIndex])
            callercountstring += "Samtools Ref Seq"
        else:
            SamtoolsRefSeqAnnotation.append("N/A")
            SamtoolsRefseqAlleles.append("N/A")
        CallerCount.append(callercountstring)
    print("Comparing annotations and alleles for each variant and extracting effect prediction level...")
    Predictedeffectlevel = []   
    AlleleDisplay = []
    for x in range(0, len(Allvariants)):
        if FreeBayesGenBankAnnotation[x] == GATKGenBankAnnotation[x] and FreeBayesGenBankAnnotation[x] == SamtoolsGenBankAnnotation[x]:
            GATKGenBankAnnotation[x] = ""
            SamtoolsGenBankAnnotation[x] = ""
        else:
            pass
        if FreeBayesRefSeqAnnotation[x] == GATKRefSeqAnnotation[x] and FreeBayesRefSeqAnnotation[x] == SamtoolsRefSeqAnnotation[x]:
            GATKRefSeqAnnotation[x] = ""
            SamtoolsRefSeqAnnotation[x] = ""
        else:
            pass
        if FreeBayesGenBankAlleles[x] == GATKGenBankAlleles[x] and FreeBayesGenBankAlleles[x] == SamtoolsGenBankAlleles[x]:
            AlleleDisplay.append(FreeBayesGenBankAlleles[x])
        elif FreeBayesGenBankAlleles[x] == "N/A" and GATKGenBankAlleles[x] == SamtoolsGenBankAlleles[x]:
            AlleleDisplay.append(GATKGenBankAlleles[x])
        elif GATKGenBankAlleles[x] == "N/A" and FreeBayesGenBankAlleles[x] == SamtoolsGenBankAlleles[x]:
            AlleleDisplay.append(SamtoolsGenBankAlleles[x])
        elif SamtoolsGenBankAlleles[x] == "N/A" and FreeBayesGenBankAlleles[x] == GATKGenBankAlleles[x]:
            AlleleDisplay.append(FreeBayesGenBankAlleles[x])
        elif FreeBayesGenBankAlleles[x] == "N/A" and GATKGenBankAlleles[x] == "N/A":
            AlleleDisplay.append(SamtoolsGenBankAlleles[x])
        elif FreeBayesGenBankAlleles[x] == "N/A" and SamtoolsGenBankAlleles[x] == "N/A":
            AlleleDisplay.append(GATKGenBankAlleles[x])
        elif GATKGenBankAlleles[x] == "N/A" and SamtoolsGenBankAlleles[x] == "N/A":
            AlleleDisplay.append(FreeBayesGenBankAlleles[x])
        else:
            AlleleDisplay.append(FreeBayesGenBankAlleles[x] + " | " + GATKGenBankAlleles[x] + " | " + SamtoolsGenBankAlleles[x])


        EffectPredictionString = FreeBayesGenBankAnnotation[x] + GATKGenBankAnnotation[x] + SamtoolsGenBankAnnotation[x] + FreeBayesRefSeqAnnotation[x] + GATKRefSeqAnnotation[x] + SamtoolsRefSeqAnnotation[x]
        
        spliteps = EffectPredictionString.split("|")
        genename.append(spliteps[3])
        effect.append(spliteps[1])
        
        strsmod = ""
        if "HIGH" in EffectPredictionString:
            strsmod += "HIGH, "
            HighVariants.append(Allvariants[x])
        if "MODERATE" in EffectPredictionString:
            strsmod += "MODERATE, "
            ModerateVariants.append(Allvariants[x])
        if "LOW" in EffectPredictionString:
            strsmod += "LOW, "
            LowVariants.append(Allvariants[x])
        if "MODIFIER" in EffectPredictionString:
            strsmod += "MODIFIER"
            ModifierVariants.append(Allvariants[x])
        Predictedeffectlevel.append(strsmod)


    SortingVariants = []
    for x in range(0,len(Allvariants)):
        TemporaryItem = [x, (Allvariants[x]).replace("chr","").replace("_","").replace("X","32")]
        SortingVariants.append(TemporaryItem)

    SortingVariants.sort(key= lambda SortingVariants:SortingVariants[1])
    for x in range(0,len(SortingVariants)):
        SortingVariants[x].append(0)
    NearByDistance = 10
    VariantNearBy = []

    for x in range(0,len(SortingVariants)-1):
        if abs(int(SortingVariants[x+1][1]) - int(SortingVariants[x][1])) < int(NearByDistance):
            SortingVariants[x].append(1)
            SortingVariants[x+1].append(1)
        else:
            SortingVariants[x+1].append(0)

    SortingVariants.sort(key= lambda SortingVariants:SortingVariants[0])
    VariantNearByType = []
    for x in range(0,len(SortingVariants)):
        VariantNearBy.append(SortingVariants[x][-1])
        if CallerCount[x].count(";") > 2 and VariantNearBy[x] == 1:
            VariantNearByType.append("B")
        elif CallerCount[x].count(";") <= 2 and VariantNearBy[x] == 1:
            VariantNearByType.append("A")
        else:
            VariantNearByType.append("")
        if Predictedeffectlevel[x][-1] == ", ":
            Predictedeffectlevel[x] = Predictedeffectlevel[x][:-2]
    HighVariantsSet = set(HighVariants)
    ModerateVariantsSet = set(ModerateVariants)
    LowVariantsSet = set(LowVariants)
    ModifierVariantsSet = set(ModifierVariants)
    now = datetime.now()
    SheetName = (str(now) + "compdata.xlsx").replace("/","_")


    print("Creating spreadsheet...")
    df = pd.DataFrame({'Chromosome Location': Allvariants,
                    "Found by Callers:": CallerCount,
                    "Ref to Alt":AlleleDisplay ,
                    "Freebayes GB Ann": FreeBayesGenBankAnnotation,
                    "GATK GB Ann": GATKGenBankAnnotation, 
                    "Samtools GB Ann": SamtoolsGenBankAnnotation,
                    "Freebayes RS Ann": FreeBayesRefSeqAnnotation, 
                    "GATK RS Ann": GATKRefSeqAnnotation, 
                    "Samtools RS Ann":SamtoolsRefSeqAnnotation, 
                    "Predicted Effect Level":Predictedeffectlevel,
                    "Variant Nearby?":VariantNearBy, 
                    "Variant Nearby Type":VariantNearByType, 
                    "Freebayes Case Reference Allele Frequency": FreebayesRefAlleleFreq, 
                    "Freebayes Case Alt 1 Frequency":FreebayesAlt1Freq, 
                    "Freebayes Case Alt 2 Frequency":FreebayesAlt2Freq, 
                    "GATK Case Reference Allele Frequency":GATKRefAlleleFreq, 
                    "GATK Case Alt 1 Freq":GATKAlt1Freq, 
                    "GATK Case Alt 2 Freq":GATKAlt2Freq,
                    "Samtools Case Reference Allele Frequency":SamtoolsRefAlleleFreq, 
                    "Samtools Case Alt 1 Freq":SamtoolsAlt1Freq, 
                    "Samtools Case Alt 2 Freq":SamtoolsAlt2Freq,
                    "Gene Name": genename,
                    "Effect": effect
                    })
    
    with open("SNPsiftPredictions.txt", "w") as sift:
        for i in range(len(Allvariants)):
            if "MODERATE" in Predictedeffectlevel[i] and len(AlleleDisplay[i])== 6:
                addy = Allvariants[i].replace("chr","")
                allele = AlleleDisplay[i].split(" -> ")
                addy = addy.split("_")
                siftline = str(addy[0]+" " +addy[1]+" "+ addy[1]+" "+ allele[0]+"/"+allele[1]+ " 1"+"\n")
                sift.write(siftline)
    dfvars = df['Chromosome Location'].values.tolist()
    fbdfvars = list(set(dfvars) & set(FreeBayesVariants))
    gadfvars = list(set(dfvars) & set(GATKVariants))
    smdfvars = list(set(dfvars) & set(SamtoolsVariants))
    df.to_excel(SheetName, sheet_name='Main Data', index=False)
    print("Creating Venn Diagram...")
    venn3([set(fbdfvars), set(gadfvars), set(smdfvars)], ['FreeBayes', 'GATK', 'Samtools'])
    plt.title("Variants-Per-Caller Venn Diagram")
    Venname = (str(now) + "venndiagram.jpg")
    plt.savefig(Venname, dpi = 300)
    labels = venn.get_labels([HighVariantsSet, ModerateVariantsSet,
                               LowVariantsSet, ModifierVariantsSet], 
                               fill=['number'])
    fig, ax = venn.venn4(labels, names=['High', 'Moderate', 'Low', 'Modifier'])
    LevelVennName = (str(now) + "levelvenn.jpg")
    plt.title("Occurence of SnpEff Predicted Levels")
    plt.savefig(LevelVennName, dpi = 300)
    print("Done!")
