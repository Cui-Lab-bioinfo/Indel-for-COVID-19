import os.path
import pickle,csv
from params import *
# max_num_clades = 3000
# dkind = 1
# rsnp = 0
# aa2muc = 2#0 for aa 1 for nuc 2 for combine
# add_indel = 1
# res_thresh = 4
# INDEL_weight = 50
# # max_loss = 0
# # max_loss = 15810344 #for indel1 resv 10
# # max_loss = 15810344 #for indel0 resv 10
# # max_loss = 17656428 #for indel0 resv 6
# max_loss = 18558864 #for indel1 resv 4
#
#
# suffix = "{}.dkind{}.rsnp{}.aa2muc{}.indel{}.resv{}.indelweight{}.maxloss{}"\
#     .format(max_num_clades,dkind,rsnp,aa2muc,add_indel,res_thresh,INDEL_weight,max_loss)
# clade2id : dict{clade:{"seq_id": list, "growth_rate":growth_rate,"pango": pango}}
def clade2accesionID(single_file,columns_file,outfile):
    clade2accesion = {}
    pango2accesion = {}
    seq2clade={}
    print(f"growth file {single_file} \ncloumn file {columns_file}")
    with open(single_file,"r") as f_csv, open(columns_file, "rb") as f_column, open(outfile, "wb") as f_out:
        columns = pickle.load(f_column)  # clade, index,day,location,lineage
        sample_count = 0
        for s_id in range(len(columns["clade"])):
            clade = columns["clade"][s_id].strip()#remove space
            indx = columns["index"][s_id]
            if not clade in clade2accesion:
                clade2accesion.update({clade: {}})
                clade2accesion[clade].update({"seq_id":[]})
            clade2accesion[clade]["seq_id"].append(indx)
            if not indx in seq2clade:
                seq2clade.update({indx:clade})
            else:
                pass
                # print(f"index:{indx} clade0 {seq2clade[indx]} clade1 {clade}")
            sample_count += 1

        spamreader = csv.reader(f_csv)
        f_csv2 = open(single_file, "r")
        len_csv = len(f_csv2.readlines())
        f_csv2.close()
        len_column = len(set(columns["clade"]))
        print(f"len column {len_column} len growth file {len_csv}")
        for row_ind, row_list in enumerate(spamreader):
            # if row_ind == 0:
            #     print(f"ignore row {row_list}")
            #     continue
            # print(f"{row_ind} {row_list}")
            pango = row_list[1].strip()
            growth_rate = row_list[2][1:]
            clade = row_list[4].strip()
            # assert clade in clade2accesion, f"Wrong clade, {clade} not in column file but in growth file"
            if clade in clade2accesion:
                clade2accesion[clade].update({"growth_rate":growth_rate})
                clade2accesion[clade].update({"pango": pango})
            # clade2accesion[clade].update({"growth_rate": growth_rate})
            # growth_list.append(row_list)

        #convert to pango->clade list
        clade_list = []
        for clade, clade_info in clade2accesion.items():
            # clade_dict = {}
            # print(clade_info.keys())
            if not "pango" in clade_info:
                print(f"clade: {clade} not exist pango")
                print(f"seq id {clade_info['seq_id']}")
                continue
            pango = clade_info["pango"]
            if not pango in pango2accesion:
                pango2accesion.update({pango:{}})
            clade_info.pop("pango")
            # clade_info.update({"clade":clade})
            clade_dict = {clade:clade_info}
            pango2accesion[pango].update(clade_dict)
            clade_list.append(clade)

        data = {"clade2accesion":clade2accesion,"pango2accesion":pango2accesion,'seq2clade':seq2clade}
        print(f"num clade {len(clade_list)}, num sample {sample_count}, num seq {len(seq2clade)}")
        pickle.dump(data, f_out)

if __name__ == "__main__":
    single_dir =    "/home/yanhongliang/projects/results/summary"
    column_dir =    "/home/yanhongliang/projects/results/rawdata"
    out_dir =       "/home/yanhongliang/projects/results/clade2accesionid"
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    file_list = os.listdir(single_dir)
    for single_file in file_list:
        if single_file.startswith("singlelineage") and not single_file.startswith("singlelineage.growthrates"):
            # suffix = single_file[14:-4]
            suffix = ".".join(single_file.split('.')[1:-1])
            print(f"clade2id in single file {single_file} suffix {suffix}")
            # growth_rate = f"singlelineage.growthrates.lobal-median.{suffix}.csv"
            columns_file = os.path.join(column_dir,f"columns.{suffix}.pkl")

            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            out_file = os.path.join(out_dir,f"clade2accesion.{suffix}.pkl")

            clade2accesionID(os.path.join(single_dir,single_file), columns_file, out_file)

