import os
import sys
file_1  = sys.argv[1]
file_2 = sys.argv[2]
f3 = open(sys.argv[3],'r')
f1 = open(file_1,'r')
f2 = open(file_2,'w')
head = f1.readline()
m = f1.readlines()
gene_dict = {}
end_dict = {}
strand_dict = {}
for line in m:
    line_ = line.strip("\n").split("\t")
    if line_[2] == "gene" or line_[2] == "pseudogene":
        gene_dict[line_[3]] = []
        end_dict[line_[3]] = line_[4]
        strand_dict[line_[3]] = line_[6]
for line in m:
    line_ = line.strip("\n").split("\t")
    if line_[2] == "gene" or line_[2] == "pseudogene":
        key = line_[3]
        gene_dict[key].append(line)
    else:
        gene_dict[key].append(line)

sorted_data_keys_str = {k: v for k, v in sorted(gene_dict.items(), key=lambda item: int(item[0]))}
caga_gff = {}
def pre_find(key_list,key_post):
    for key in key_list.keys():
        if int(key_post) > int(key):
            key_pre = key
        elif key_post == key:
            return(key_pre)
        else:
            return('wukong')
def search_gene(TSS,gene_body):
    for start in gene_body.keys():
        end = end_dict[start]
        strand_value = strand_dict[start]
        if strand_value == "+":
            if int(TSS) >= int(start) and int(TSS) <= int(end):
                return(start)
            elif int(TSS) < int(start):
                return(start)
        if strand_value == "-":
            if int(TSS) >= int(start) and int(TSS) <= int(end):
                return(start)
            elif int(TSS) < int(start):
                start_pre = pre_find(gene_body,start)
                return(start_pre)
def cp_annot(gene_struc):
    for key in gene_struc:
        f2.write(key)
def correct_annot(start1,gene_struc):
    #for line in gene_struc:
     #   line_ = line.strip("\n")
     #   if line_[2] == "gene" or line_[2] == "pseudogene":
     #       gene_start = int(line_[3])
     #       gene_end   = int(line_[4])
    #print(gene_end)
    TSS = caga_gff[start1]

    gene_start = int(start1)
    gene_end  = int(end_dict[start1])
    for line in gene_struc:
        line_ = line.strip('\n').split("\t")
        strand = line_[6]
        if strand == "+":
            start = line_[3]
            end   = line_[4]
            if int(TSS) >= int(start) and int(TSS) <= int(end) and strand_caga[TSS] == "+":
                line_[3] = TSS
                anotation = line_[8] + ";add_dominant"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            elif int(TSS) < int(start) and int(start) == gene_start and strand_caga[TSS] == "+":
                line_[3] = TSS
                anotation = line_[8] + ";add_dominant_out"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            else:
                f2.write(line)
        if strand == "-":
            start = line_[3]
            end =   line_[4]
            if int(TSS) >= int(start) and int(TSS) <= int(end) and strand_caga[TSS] == "-":
                line_[4] = TSS
                anotation = line_[8] + ";add_dominant"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            elif int(TSS) > int(end) and int(end) == gene_end and strand_caga[TSS] == "-":
                line_[4] = TSS
                anotation = line_[8] + ";add_dominant_out"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            else:
                f2.write(line)





def correct_annot_2_member(start1,TSS,gene_struc):
    
    gene_start = int(start1)
    gene_end  = int(end_dict[start1])
    for line in gene_struc:
        line_ = line.strip('\n').split("\t")
        strand = line_[6]
        if strand == "+":
            start = line_[3]
            end   = line_[4]
            if int(TSS) >= int(start) and int(TSS) <= int(end) and strand_caga[TSS] == "+":
                line_[3] = TSS
                anotation = line_[8] + ";add_dominant_maybe_two_start"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            elif int(TSS) < int(start) and int(start) == gene_start and strand_caga[TSS] == "+":
                line_[3] = TSS
                anotation = line_[8] + ";add_dominant_out_maybe_two_start"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            else:
                f2.write(line)
        if strand == "-":
            start = line_[3]
            end =   line_[4]
            if int(TSS) >= int(start) and int(TSS) <= int(end) and strand_caga[TSS] == "-":
                line_[4] = TSS
                anotation = line_[8] + ";add_dominant_maybe_two_start"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            elif int(TSS) > int(end) and int(end) == gene_end and strand_caga[TSS] == "-":
                line_[4] = TSS
                anotation = line_[8] + ";add_dominant_out_maybe_two_start"
                line_[8] = anotation
                line_correct = "\t".join(line_) + "\n"
                f2.write(line_correct)
            else:
                f2.write(line)


n = f3.readlines()
strand_caga={}
for line in n:
    line_ = line.strip("\n").split("\t")
    dom = line_[7]
    strand = line_[4]
    strand_caga[dom] = strand






for line_caga in n:
    key = int(line_caga.strip("\n").split("\t")[7])
    wind = search_gene(key,sorted_data_keys_str) 
    if not(wind in caga_gff.keys()):
        caga_gff[wind] = str(key)
    else:
        first_value = caga_gff[wind] 
        caga_gff[wind] = first_value + "_"+ str(key)
    



f2.write(head)    
for key in sorted_data_keys_str.keys():

    if key in caga_gff.keys():
        TSS = caga_gff[key]
        if  not("_" in TSS):
            correct_annot(key,sorted_data_keys_str[key])
        else:
            key_list = TSS.split('_')
            for key_item in key_list:
                correct_annot_2_member(key,key_item,sorted_data_keys_str[key]) 
    else:
        cp_annot(sorted_data_keys_str[key])
f1.close()
f2.close()
f3.close()




        






