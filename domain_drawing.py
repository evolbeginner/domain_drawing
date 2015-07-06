#! /bin/env python

# used to draw domain structure for protein sequences given CDD formated file or tab-delimited file (see sample files)
# written by Sishuo Wang from Department of Botany, The University of British Columbia
# You are welcome to write e-mails to sishuowang@hotmail.ca if you have any questions and/or suggestions. Your help will be highly appreciated.


##############################################################
import getopt
import os
import sys
import re
from biograpy import Panel, tracks, features
from Bio import SeqFeature
from Bio import SeqIO


##############################################################
input = None
output = None
outdir = None
domain_objs = []
hit_types = []
seq_file=None
given_seq_length = None
file_type = None
is_domain_on_same_line = True
is_domain_on_diff_line = False
is_gene_name_replace_protein_name = False
force = False
test_track = None

colors = ['green', 'yellow', 'blue', 'red', 'cyan', 'orange', 'purple', 'grey', 'black', 'pink', 'cyan']
domain_colors = {}


##############################################################
class Domain:
    def __init__(self):
       pass 


def read_domain_file(input, hit_types=[]):
    domain_objs = {}
    domain_info = {}
    for line in open(input, 'r'):
        line = line.rstrip('\n\r')
        line_array = line.split('\t')
        line_array = filter(lambda x: x, line_array) # Note that there is one more tab after 1st tab. So it needs to be deleted.
        # AT1G64840.1		384	HMMPfam	PF03478	DUF295	286	336	4.5000000000000035E-7	28-Oct-2010	IPR005174	Protein of unknown function DUF29
        protein_name, length, hit_type, domain_name, start, stop = line_array[0], line_array[1], line_array[2], line_array[4], line_array[5], line_array[6]

        if hit_types and (not hit_type in hit_types):
            continue

        if not protein_name in domain_objs:
            domain_objs[protein_name] = []
        else:
            pass

        domain_obj = Domain()
        domain_obj.start = int(start)
        domain_obj.stop  = int(stop)
        domain_obj.short_name  = domain_name
        domain_objs[protein_name].append(domain_obj)
        domain_info[protein_name] = domain_objs[protein_name]
    return domain_info


def read_CDD(input, hit_types=[]):
    CDD_domain_objs = []
    with open(input,'r') as fh:
        for line in fh:
            if line.startswith('Query'):
                for line in fh:
                    if re.search('^$',line):
                        break
                    line = line.rstrip('\n\r')
                    # Query	Hit type	PSSM-ID	From	To	E-Value	Bitscore	Accession	Short name	Incomplete	Superfamily
                    line_array = line.split('\t')
                    hit_type, start, stop, short_name = line_array[1], line_array[3], line_array[4], line_array[8]
                    if hit_types and (not hit_type in hit_types):
                        continue
                    domain_obj = Domain()
                    domain_obj.start = int(start)
                    domain_obj.stop  = int(stop)
                    domain_obj.short_name  = short_name
                    CDD_domain_objs.append(domain_obj)
    return CDD_domain_objs


def get_length_from_seq_file(seq_file, seq_file_format='fasta'):
    class Seq:
        def __init__(self):
            self.length = len(seq_record.seq)

    seq_objs={}
    handle = open(seq_file, "r")
    for seq_record in SeqIO.parse(handle, seq_file_format):
        seq_objs[seq_record.id] = Seq()
    handle.close()
    return(seq_objs)


def get_color_for_each_domain(domain_info):
    domainsh = {}
    for gene_name,domain_objs in domain_info.iteritems():
        for domain_obj in domain_objs:
            short_name = domain_obj.short_name
            if not short_name in domainsh:
                domainsh[short_name] = 1
            else:
                domainsh[short_name] += 1
    return(domainsh)


def make_dir_force(dir,force=False):
    import os
    class MyError( Exception ): pass
    if os.path.exists(dir):
        if force:
            os.system("rm -rf " + dir + " 2>/dev/null")
            os.system("mkdir -p " + dir)
        else:
            raise MyError('dir ' + dir + ' has already existed!')
    else:
        os.system("mkdir -p " + dir)


###############################################################
try:
    opts,args = getopt.getopt(
        sys.argv[1:],
        "hi:o:",
        ["in=","out=","outdir=","force","type=",
        "hit_type=","seq_title=","seq_file=","length=","seq_length=","file_type=",
        "domain_on_diff_line","domain_on_same_line","gene_name_replace_protein_name"]
    )
except getopt.GetoptError:
    print "illegal error!"
    show_help()

for op, value in opts:
    if op == "-i" or op == "--in":
        input = value
    elif op == '-o' or op == '--out':
        output = value
    elif op == "--outdir":
        outdir = value
    elif op == "--force":
        force = True
    elif op == '--type':
        type = value
    elif op == "--hit_type":
        hit_types.append(value)
    elif op == '--seq_title':
        seq_title = value
    elif op == '--seq_file':
        seq_file = value
    elif op == '--length' or op == "--seq_length":
        given_seq_length = int(value)
    elif op == '--file_type':
        file_type = value
    elif op == '--domain_on_same_line':
        is_domain_on_same_line = True
        is_domain_on_diff_line = False
    elif op == '--domain_on_diff_line':
        is_domain_on_same_line = False
        is_domain_on_diff_line = True
    elif op == 'gene_name_replace_protein_name':
        is_gene_name_replace_protein_name = True
    elif op == "-h":
        show_help()

if not outdir:
    outdir = "."

if outdir == ".":
    pass
else:
    make_dir_force(outdir,force)

###############################################################
if file_type:
    if file_type == "CDD": 
        domain_objs = read_CDD(input, hit_types)
    elif file_type == "TAIR" or file_type == "normal":
        domain_info = read_domain_file(input, hit_types)
else:
    print "file_type has to be given! Exiting ......"
    sys.exit()


domainsh = get_color_for_each_domain(domain_info)
domains = sorted(domainsh.keys(), key = lambda x: domainsh[x])
for index,domain in enumerate(domains):
    domain_colors[domain] = colors[index]


##########################################################
counter=0

for seq_title,domain_objs in domain_info.iteritems():
    counter+=1
    feats = []
    domain_names = []
    output = None

    if given_seq_length:
        seq_length = given_seq_length
    elif seq_file:
        seq_objs = get_length_from_seq_file(seq_file,'fasta')
        if seq_title in seq_objs:
            seq_length = seq_objs[seq_title].length
        else:
            raise "seq_length cannot be obtained"
    else:
        raise "length has to be given by either --length, --seq_length or --seq"

    if is_gene_name_replace_protein_name:
        m=re.search('([^.]+)', seq_title)
        seq_title = m.group(1)
    else:
        pass
    print seq_title

    for domain_obj in domain_objs:
        print domain_obj.start

    # draw domain structure using BioGraPy
    panel = Panel(fig_dpi=100)

    #if not test_track:
    test_track = tracks.BaseTrack(name = seq_title, sort_by = None)

    for index,domain_obj in enumerate(domain_objs):
        dfeat = None
        feat = SeqFeature.SeqFeature()
        start, stop = domain_obj.start, domain_obj.stop
        #feat.location = SeqFeature.FeatureLocation(start, stop)
        feat.location = SeqFeature.FeatureLocation(1, 100)
        feats.append(feat)
        domain_names.append(domain_obj.short_name)
        if is_domain_on_same_line:
            continue
        if is_domain_on_diff_line:
            colors_4_this_domain = [domain_colors[domain_obj.short_name]]
            colors_4_this_domain = 'y'
            if index == 0:
                dfeat = features.DomainFeature([feat], name=domain_obj.short_name, seq_line=(1,seq_length), fc=colors_4_this_domain)
                dfeat = features.DomainFeature([feat], name=domain_obj.short_name, fc=12)
            else:
                dfeat = features.DomainFeature([feat], name=domain_obj.short_name, fc=colors_4_this_domain)
            dfeat = features.DomainFeature([feat], name=domain_obj.short_name)
            test_track.append(dfeat)

    if is_domain_on_same_line:
        #dfeat = features.DomainFeature(feats, name=domain_names, seq_line=(1,seq_length), fc=colors[0:len(feats)])
        colors_4_this_domain = map(lambda x: domain_colors[x], domain_names)
        print domain_names
        dfeat = features.DomainFeature(feats, name=domain_names, seq_line=(1,seq_length),
                    #ec = 'None', fc=colors_4_this_domain, color_by_cm = False)
                    ec = 'None', color_by_cm = False)
        test_track.append(dfeat)

    panel.add_track(test_track)
    #panel.tracks = [test_track]

    if not output:
        output = os.path.join(outdir, seq_title+'.png')
    #if counter==6:  print output
    print output
    panel.save(output)


