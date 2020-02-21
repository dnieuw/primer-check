#!/usr/bin/env python3

import os
import time
import argparse
import re
import multiprocessing
try:
    import cPickle as pickle
except ImportError:
    import pickle
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Align import PairwiseAligner

def score_oligo(query, fasta):
    AMBIGUOUS_DNA_LETTERS = 'RYWSMKHBDV'
    
    testing = query.id
    
    query = query.seq
    fasta = fasta.seq
    
    #Score for a perfect match
    maxScore = len(query)
    
    amb_pos = []
    #Find ambiguous positions in query
    for i in range(len(query)):
        if query[i] in AMBIGUOUS_DNA_LETTERS:
            amb_pos.append(i)

    amb_pos_ref = []
    #Find ambiguous positions in fasta
    for i in range(len(fasta)):
        if fasta[i] in AMBIGUOUS_DNA_LETTERS:
            amb_pos_ref.append(i)
            
    #Setup aligner object
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -100
    aligner.extend_gap_score = -100
    
    #No ambiguous bases in ref or query
    if (len(amb_pos) == 0 and len(amb_pos_ref) == 0):
        #Perform pairwise alignment between probe and reference   
        alignment = aligner.align(fasta, query)[0]

    #No ambiguous bases in ref
    elif (len(amb_pos_ref) == 0):
        #Handle all the possible variants of the primer and find the best match
        current_alignment = aligner.align(fasta, query)[0]
        variant = query
        #Loop through all variants and find the best alignment
        for pos in amb_pos:
            for base in ambiguous_dna_values[query[pos]]:
                variant = str(variant)[:pos] + base + str(variant)[pos + 1:]
                alignment = aligner.align(fasta, variant)[0]
                if alignment.score > current_alignment.score:
                    current_alignment = alignment
                    break
        alignment = current_alignment

    #No ambiguous bases in query
    elif (len(amb_pos) == 0):
        #Handle all the possible variants of the primer and find the best match
        current_alignment = aligner.align(fasta, query)[0]
        variant = fasta
        #Loop through all variants and find the best alignment
        for pos in amb_pos_ref:
            for base in ambiguous_dna_values[fasta[pos]]:
                variant = str(variant)[:pos] + base + str(variant)[pos + 1:]
                alignment = aligner.align(variant, query)[0]
                if alignment.score > current_alignment.score:
                    current_alignment = alignment
                    break
        alignment = current_alignment
        
    #Ambiguous bases in ref and query
    else:
        #Handle all the possible variants of the primer and the reference and find the best match
        current_alignment = aligner.align(fasta, query)[0]
        variant_ref = fasta
        variant_query = query
        for pos in amb_pos_ref:
            for base in ambiguous_dna_values[fasta[pos]]:
                variant_ref = str(variant_ref)[:pos] + base + str(variant_ref)[pos + 1:]
                for pos in amb_pos:
                    for base in ambiguous_dna_values[query[pos]]:
                        variant_query = str(variant_query)[:pos] + base + str(variant_query)[pos + 1:]
                        alignment = aligner.align(variant_ref, variant_query)[0]
                        if alignment.score > current_alignment.score:
                            current_alignment = alignment
                            break
        alignment = current_alignment
        
    return(alignment)

def get_references(filename):
    try:
        pfile = open(filename, "r")
    except Exception as e:
        print(e)
        return(None)
    else:
        with open(filename, 'r') as ref_file:
            ref_dict = {}
            for record in SeqIO.parse(ref_file, "fasta"):
                #Remove gaps from records (when fasta is aligned)
                record.seq = str(record.seq).replace('-','N')
                record.seq = record.seq.upper()
                ref_dict[record.id] = record
            return(ref_dict)

def get_probes(filename):
    #Try to open the probe file
    try:
        pfile = open(filename, "r")
    except Exception as e:
        print(e)
        return(None)

    #Read probes to list
    probes = {}
    # Name;Sequence;Type;Origin;Gene;Set
    _ = pfile.readline()
    for line in pfile:
        pname, pseq, ptype, porigin, pgene, pset = line.strip().split(";")
        probeid = '_'.join([porigin, pgene, ptype, pset])
        record = SeqRecord(Seq(pseq.upper(), IUPAC.ambiguous_dna), id=probeid, description=line.strip())
        if ptype == "RE":
            record.seq = record.seq.reverse_complement()
        probes[probeid] = record
    pfile.close()
    return(probes)

def check_old_result_and_create_worklist(probelist, reflist, database):
    #Load file containing results for previous primer checks
    # oligo: {family: "", type: "", records: [], stats: []}
    checked_oligos = {}
    if os.path.exists(database):
        with open(database, "rb") as fp:
            checked_oligos = pickle.load(fp)
            
    worklist = []
    #Check if probe is in the results already
    for probe in probelist.keys():
        if probe in checked_oligos.keys():
            #If there is a new reference, match it to the current probe
            for ref in reflist.keys():
                if not ref in checked_oligos[probe]['alignments'].keys():
                    worklist.append((probelist[probe],reflist[ref]))
        else:
            #If the probe is not present, match it to all references
            for ref in reflist.keys():
                worklist.append((probelist[probe],reflist[ref]))
    return(checked_oligos, worklist)

def run_primer_check(probes, fasta, out, database, cores, loaddb):
    #Load an old database and write to tsv file
    if loaddb:
        if os.path.exists(database):
            with open(database, "rb") as fp:
                checked_oligos = pickle.load(fp)
            i=0
            #write out all old and new results        
            with open(out,"w") as outfile:
                for probe, result in checked_oligos.items():
                    for ref, aln in result['alignments'].items():
                        i+=1
                        print('\t'.join(result['description'].split(';')+[ref, str(aln[0]), str(aln[1]), str(aln[2]), aln[3], aln[4], aln[5]]), file=outfile)
            print("Loaded "+str(i)+" alignments")
        else:
            print("Database "+database+" does not exist")
        return(None)
    
    probes = get_probes(probes)
    if not probes:
        return(None)
    references = get_references(fasta)
    if not references:
        return(None)
    
    print(len(probes), " primers loaded")
    print(len(references), " references loaded")
        
    old_results, worklist = check_old_result_and_create_worklist(probes, references, database)
    
    print(len(worklist)," alignments to make")
    
    if not worklist:
        #write out all old and new results 
        with open(out,"w") as outfile:
            for probe, result in old_results.items():
                for ref, aln in result['alignments'].items():
                    print('\t'.join(result['description'].split(';')+[ref, str(aln[0]), str(aln[1]), str(aln[2]), aln[3], aln[4], aln[5]]), file=outfile)
        print("All probes and references are already checked, no work to be done!")
        return(None)
        
    with multiprocessing.Pool(processes=int(cores)) as p:
        resultlist = p.starmap(score_oligo, iter(worklist))
        
    fail, succes = 0,0
    for work, res in zip(worklist, resultlist):
        probe, ref = work
        #If the score is the same at the primer length the alignment was successfull
        if res.score == len(probe):
            succes+=1
        else:
            fail+=1
              
        #Split alignment string
        aln_ref, aln_str, aln_probe,_ = str(res).split("\n")
            
        #Find the aligned region
        m = re.search("[ACTGRYWSMKHBDVXN-]+", aln_probe).span()
            
        #Creat list of aligned-region substrings and score
        aln_substr = [res.score, m[0], m[1], aln_ref[m[0]:m[1]],
                      aln_str[m[0]:m[1]],
                      str(probe.seq)]
    
        #Add new results to the old results or a create new entry for a new primer
        if probe.id in old_results:
            old_results[probe.id]['alignments'][ref.id] = aln_substr
        else:
            old_results[probe.id] = {'description': probe.description,
                                      'alignments': {ref.id:aln_substr}}
    #write out all old and new results 
    with open(out,"w") as outfile:
        for probe, result in old_results.items():
            for ref, aln in result['alignments'].items():
                print('\t'.join(result['description'].split(';')+[ref, str(aln[0]), str(aln[1]), str(aln[2]), aln[3], aln[4], aln[5]]), file=outfile)
        
    print(succes," successfull alignments")
    print(fail," failed alignments")

    #Write new completed primer checks to file
    with open(database, "wb") as fp:
        pickle.dump(old_results, fp)