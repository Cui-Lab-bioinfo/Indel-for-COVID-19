import os, argparse,csv,torch
import numpy as np

# def compare(resv_list = [3,10]):
#     """
#
#     :param resv_list: file name
#     :param res:
#     :return:
#     """
#     file_patrn = "avglineage.growthrates.obal-median.3000.dkind1.rsnp0.aa2muc2.indel{}.resv{}.pt"
#
#     file_name1 = file_patrn.format(1,resv_list[0])
#     file_name2 = file_patrn.format(0,resv_list[1])
#     resfile = "compare{}_{}".format(resv_list[0],resv_list[1])
#
#     with open(file_name1,"r") as filein1, open(file_name2,"r") as filein2, open(res_file,"w") as fileout:
#         spamreader1 = csv.reader(filein1)
#         spamreader2 = csv.reader(filein2)
#
#         spamwriter = csv.writer(fileout)
#
#         for row_ind, row_list in enumerate(spamreader):
#             if row_ind == 0:
#                 continue
#             print(f"{row_ind} {row_list}")
#             growth_list.append(row_list)



def growthrate_clade2lineage(file_in, raw_out,file_out,file_dir,out_dir):
    """

    :param file_in: clade growth rate file
    :param raw_out: sort raw clade
    :param file_out: lineage growth rate file
    :return:
    """
    growth_list = []
    suffix = ".".join(file_in.split('.')[2:-1])
    print(f"suffix {suffix}")
    mutrans_file = f"mutrans.{suffix}_new.pt"
    mutrans = torch.load(os.path.join(file_dir,'mutrans',mutrans_file))
    clade_ids = mutrans['clade_id']
    clade_list = [0]*len(clade_ids)
    for k,idx in clade_ids.items():
        clade_list[idx] = k

    with open(os.path.join(file_dir,'summary',file_in),"r") as csv_f:
        spamreader = csv.reader(csv_f)
        for row_ind, row_list in enumerate(spamreader):
            if row_ind == 0:
                continue
            # print(f"{row_ind} {row_list}")
            row_list[-1] = clade_list[row_ind-1]
            assert row_list[-1] == clade_list[row_ind-1], f"clade in summary{row_list[-1]} clade in mutran {clade_list[row_ind-1]}"
            growth_list.append(row_list)


    sorted_list = sorted(growth_list,key=lambda x: float(x[1]), reverse=True)
    if raw_out == None:
        file_out = "singlelineage."+".".join(file_in.split('.')[2:-1])+".csv"
    with open(os.path.join(out_dir,file_out),"w") as csv_f:
        spamwriter = csv.writer(csv_f)
        for rank,clade_growth in enumerate(sorted_list):
            spamwriter.writerow(["{0:10}".format(rank),"{0:10}".format(clade_growth[3]),\
                       ".%9f"%(float(clade_growth[1])),".%9f"%(float(clade_growth[2])),clade_growth[4]])


    lineage_growth = {}
    with open(os.path.join(file_dir,'summary',file_in),"r") as csv_f:
        spamreader = csv.reader(csv_f)
        for row_ind, row_list in enumerate(spamreader):
            if row_ind == 0:
                continue
            # print(f"{row_ind} {row_list}")
            clade_id, grate, std, lineage, clade_name = row_list
            # clade_id = clade_id
            grate = float(grate)
            std = float(std)

            if not (lineage in lineage_growth):
                lineage_growth[lineage] = {}
                lineage_growth[lineage]["growth_rate"] = []
                lineage_growth[lineage]["std"] = []
                lineage_growth[lineage]["clade"] = []
                lineage_growth[lineage]["clade_name"] = []

            lineage_growth[lineage]["growth_rate"].append(grate)
            lineage_growth[lineage]["std"].append(std)
            lineage_growth[lineage]["clade"].append(clade_id)
            lineage_growth[lineage]["clade_name"].append(clade_name)

    print(f"\t\t******computer growth rate for {len(lineage_growth)} lineages*************")
    stat_dict = {}
    for lineage in lineage_growth:
        growth_rate = np.mean(lineage_growth[lineage]["growth_rate"])
        std = np.mean(lineage_growth[lineage]["std"])
        stat_dict.update({lineage:[growth_rate,std]})

    #sort the stats
    print(f"\t\t******sort lineage sgrowth rate {len(lineage_growth)} lineages*************")
    sorted_stats = dict(sorted(stat_dict.items(),key=lambda item: item[1][0],reverse=True))

    partn = "lineage: {0:10} \t growth rate: {1:20} \t std: {2:30}"
    for rank,(key,value) in enumerate(sorted_stats.items()):
        if rank >=10:
            break
        print(partn.format(key,value[0],value[1]))


    file_out = "avglineage."+".".join(file_in.split('.')[:-1])+".csv"

    with open(os.path.join(out_dir,file_out),"w") as csv_f:
        spamwriter = csv.writer(csv_f)
        for rank, (lineage,static_) in enumerate(sorted_stats.items()):
            spamwriter.writerow([str(rank),'\t' "{0:10}".format(lineage),'\t%.9f'%static_[0],"\t%.9f"%static_[1]])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="combine clade to lineage growth rate")
    # parser.add_argument("--file_in",
    #                     default="/home/yanhongliang/projects/bvas/old/growthrates.obal-median.3000.dkind1.rsnp0.aa2muc2.indel0.resv10.csv")
    parser.add_argument("--file_in",
                        default="growthrates.obal-median.3000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight50.maxloss18558864.csv")
    parser.add_argument("--file_out",
                        default=None)
    args = parser.parse_args()
    file_dir = "/home/yanhongliang/projects/results"
    out_dir = "/home/yanhongliang/projects/results/summary"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    file_list = os.listdir(os.path.join(file_dir,'summary'))
    for file_in in file_list:
        if file_in.startswith("growthrates") and ("3000" in file_in):
            growthrate_clade2lineage(file_in, None, None,file_dir,out_dir)
    # growthrate_clade2lineage(args.file_in ,None,None)
