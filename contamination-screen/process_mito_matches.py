import sys
# python process_mito_matches.py blast_input_file output_file_name
# blast_input_file is -outfmt 7, has comment lines removed, and is sorted
# (sort -k1,1 -k7n)
def openfile(filein):
    contigdict = {}
    with open(filein, 'U') as infile:
        for line in infile:
            line = line.strip()
            splitline = line.split("\t")
            #print splitline
            contig = str(splitline[0])
            st_range = int(splitline[6])
            end_range = int(splitline[7])
            if contig in contigdict:
                code = 0
                for list_elem in range(len(contigdict[contig])):
                    #print list_elem
                    #print contigdict[contig][list_elem][0]
                    #print contigdict[contig][list_elem][1]
                    if contigdict[contig][list_elem][0] <= st_range <= contigdict[contig][list_elem][1] and contigdict[contig][list_elem][0] <= end_range <= contigdict[contig][list_elem][1]:
                        #print st_range, end_range
                        code = 1
                        break
                    elif contigdict[contig][list_elem][0] <= st_range <= contigdict[contig][list_elem][1] and end_range > contigdict[contig][list_elem][1]:
                        contigdict[contig][list_elem][1] = end_range
                        code = 1
                        #print contigdict
                if code == 1:
                    pass
                elif code == 0:
                    contigdict[contig].append([st_range,end_range])
                    #print contigdict
            else:
                contigdict[contig]=[[st_range,end_range]]
                #print contigdict
    return contigdict

def printdict(indict):
    with open(sys.argv[2], 'w') as outfile:
        for elem in indict:
            for subelem in indict[elem]:
                subelem_form = '\t'.join([str(x) for x in subelem])
                myline = elem + "\t" + subelem_form + "\n"
                outfile.write(myline)
myfile = sys.argv[1]

getdict = openfile(myfile)
printdict(getdict)
