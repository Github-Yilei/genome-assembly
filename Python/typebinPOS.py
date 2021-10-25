
import sys
input_file = sys.argv[1]

chrom = ""
pos_start = 0
pos_end = 0

ref = 0
alt = 0
mis = 0
het = 0

pos_list = list()
with open(input_file, 'r') as lines:
    for line in lines:
        line_spl = line.split()
        temp = line_spl[2]
        if line_spl[0] != chrom:
            if ref != 0:
                ref_start = pos_start
                ref_end = pos_end 
                pos_list.append(chrom + "\t" + str(ref_start) + "\t" + str(ref_end) + "\t" + "ref"+ "\t" + str(ref))
                ref = 0
                
            elif alt != 0:
                alt_start = pos_start
                alt_end =  pos_end
                pos_list.append(chrom + "\t" + str(alt_start) + "\t" + str(alt_end) + "\t" + "qey" + "\t" + str(alt))
                alt = 0
                
            elif mis != 0:
                mis_start = pos_start
                mis_end =  pos_end
                pos_list.append(chrom + "\t" + str(mis_start) + "\t" + str(mis_end) + "\t" + "mis" + "\t" + str(mis))
                mis = 0
            
            elif het != 0:
                het_start = pos_start
                het_end =  pos_end
                pos_list.append(chrom + "\t" + str(het_start) + "\t" + str(het_end) + "\t" + "het" + "\t" + str(het))
                het = 0
            
            # re-setting records
            chrom =  line_spl[0]
            pos_start = line_spl[1]
            chrom =  line_spl[0]

            if temp == 'ref':
                ref += 1
            elif temp == 'alt':
                alt += 1
            elif temp == 'mis':
                mis += 1
            elif temp == 'het':
                het += 1
            
        elif line_spl[0] == chrom:    
            if temp == 'ref':
                ref += 1
                # if it is the first element of ref category 
                if ref == 1:
                    # if perior-elemnt is alt
                    if alt != 0:
                        alt_start = pos_start
                        alt_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" + str(alt_start) + "\t" + str(alt_end) + "\t" + "qey" + "\t" + str(alt))
                        alt = 0
                    elif mis != 0:
                        mis_start = pos_start
                        mis_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" +  str(mis_start) + "\t" + str(mis_end) + "\t" + "mis" + "\t" + str(mis))
                        mis = 0
                    elif het != 0:
                        het_start = pos_start
                        het_end =  pos_end
                        pos_list.append(chrom + "\t" + str(het_start) + "\t" + str(het_end) + "\t" + "het" + "\t" + str(het))
                        het = 0
                        
                    pos_start = line_spl[1]
                    chrom =  line_spl[0]

                
   
            elif temp == 'alt':
                alt += 1
                if alt == 1:
                    if ref != 0:
                        ref_start = pos_start
                        ref_end = int(line_spl[1]) - 1 
                        pos_list.append(chrom + "\t" + str(ref_start) + "\t" + str(ref_end) + "\t" + "ref"+ "\t" + str(ref))
                        ref = 0
                    elif mis != 0:
                        mis_start = pos_start
                        mis_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" + str(mis_start) + "\t" + str(mis_end) + "\t" + "mis" + "\t" + str(mis))
                        mis = 0
                    elif het != 0:
                        het_start = pos_start
                        het_end =  pos_end
                        pos_list.append(chrom + "\t" + str(het_start) + "\t" + str(het_end) + "\t" + "het" + "\t" + str(het))
                        het = 0
            
                    pos_start = line_spl[1]
                    chrom =  line_spl[0]
 
                
            elif temp == 'mis':
                mis += 1
                if mis == 1:
                    if ref != 0:
                        ref_start = pos_start
                        ref_end = int(line_spl[1]) - 1 
                        pos_list.append(chrom + "\t" + str(ref_start) + "\t" + str(ref_end) + "\t" + "ref"+ "\t" + str(ref))
                        ref = 0
                    elif alt != 0:
                        alt_start = pos_start
                        alt_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" + str(alt_start) + "\t" + str(alt_end) + "\t" + "qey" + "\t" + str(alt))
                        alt = 0
                    elif het != 0:
                        het_start = pos_start
                        het_end =  pos_end
                        pos_list.append(chrom + "\t" + str(het_start) + "\t" + str(het_end) + "\t" + "het" + "\t" + str(het))
                        het = 0
                    
                    pos_start = line_spl[1]
                    chrom =  line_spl[0]
            
            elif temp == 'het':
                het += 1
                if het == 1:
                    if ref != 0:
                        ref_start = pos_start
                        ref_end = int(line_spl[1]) - 1 
                        pos_list.append(chrom + "\t" + str(ref_start) + "\t" + str(ref_end) + "\t" + "ref"+ "\t" + str(ref))
                        ref = 0
                    elif alt != 0:
                        alt_start = pos_start
                        alt_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" + str(alt_start) + "\t" + str(alt_end) + "\t" + "qey" + "\t" + str(alt))
                        alt = 0
                    elif mis != 0:
                        mis_start = pos_start
                        mis_end =  int(line_spl[1]) - 1
                        pos_list.append(chrom + "\t" + str(mis_start) + "\t" + str(mis_end) + "\t" + "mis" + "\t" + str(mis))
                        mis = 0                 
                
                    pos_start = line_spl[1]
                    chrom =  line_spl[0]               
            
            # re-fresh pos_end
            pos_end = line_spl[1]
                
print("\n".join(pos_list))
