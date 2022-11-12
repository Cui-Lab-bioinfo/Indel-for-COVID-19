import os.path
import argparse,pickle,csv,torch

def select_clades(indel0map,indel1map,indel0summary,indel1summary,mut_file,out_file1,out_file2):
    """

    :param indel0map:
    :param indel1map:
    :param mut_file:
    :param out_file:
    :return:
    """

    with open(indel0map,'rb') as p_indel0,open(indel1map,'rb') as p_indel1, \
            open(indel0summary,'r') as csv_summary0, open(indel1summary,'r') as csv_summary1,\
            open(out_file1,"w") as csv_file1,open(out_file2,"w") as csv_file2:
        data0 = pickle.load(p_indel0)
        data1 = pickle.load(p_indel1)

        spawreader0 = csv.reader(csv_summary0)
        spawreader1 = csv.reader(csv_summary1)
        clade2gr0 = {}
        clade2gr1 = {}
        clade2line0 = {}
        clade2line1 = {}
        for raw_ind, raw_list in enumerate(spawreader0):

            growthrate = raw_list[2][1:]
            line = raw_list[1]
            clade = raw_list[4]
            clade2gr0.update({clade:growthrate})
            clade2line0.update({clade:line})

        for raw_ind, raw_list in enumerate(spawreader1):
            growthrate = raw_list[2][1:]
            line = raw_list[1]
            clade = raw_list[4]
            clade2gr1.update({clade:growthrate})
            clade2line1.update({clade:line})

        clade2seq0 = data0['clade2accesion']
        clade2seq1 = data1['clade2accesion']
        # seq2clade0 = data0['seq2clade']
        seq2clade1 = data1['seq2clade']

        len_can = lambda clade2seq, clade_can: [len(clade2seq[clade]['seq_id']) for clade in clade_can]
        growthrate_can = lambda clade2seq, clade_can: [clade2seq[clade]['growth_rate'] for clade in clade_can]

        # spamreader = csv.reader(c_diff)
        csvwriter1 = csv.writer(csv_file1)
        csvwriter2 = csv.writer(csv_file2)
        header = ['clade0', 'clade1', 'clade2', 'count0', 'count1', 'count2', 'lineage0', 'lineage1', 'lineage2','diff_mut1', 'diff_mut2','fitness0','fitness1','fitness2','all_mut1','all_mut2']
        csvwriter1.writerow(header)
        csvwriter2.writerow(header)
        # for raw_ind, raw_list in enumerate(spamreader):
        for clade_ind, fine0 in enumerate(clade2seq0):
            # fine0 = None
            # if raw_ind == 0:
            #     continue

            clade_list = []
            # count_list = []
            comp_Flag = False
            # print(raw_list)
            # fine_name, count0, count1= raw_list[0].split('\t')
            # fine_name = fine_name.strip()
            seq0 = clade2seq0[fine0]['seq_id']
            # growthrate0 = clade2seq0[fine0]['growth_rate']
            growthrate0 = clade2gr0[fine0]
            line0 = clade2line0[fine0]
            count0 = len(seq0)
            # count1 = int(count1.strip())
            # assert count0 != count1, f"clade {fine_name} count equals in indel0 and indel1"
            # if count1== 0:#not exist in indel1 fine name from clade0
            # fine0 = fine_name
            # seq0 = clade2seq0[fine0]['seq_id']
            # assert count0 == len(seq0), f"clade0 {fine0}: count0 {count0} not equal seq0 {len(seq0)}"
            clade_can1 = []
            for seq in seq0:
                clade_can1.append(seq2clade1[seq])
            clade_can1 = list(set(clade_can1))
            count_list = len_can(clade2seq1,clade_can1)
            # growthrate_list = growthrate_can(clade2gr1,clade_can1)
            growthrate_list = [clade2gr1[clade1] for clade1 in clade_can1]
            line_list = [clade2line1[clade1] for clade1 in clade_can1]
            len_can1 = sum(count_list)
            if len_can1 == count0:
                clade_list = clade_can1
                count_list.insert(0,count0)
                growthrate_list.insert(0,growthrate0)
                line_list.insert(0,line0)
                comp_Flag = True
            # else:
            #     #fine name from clade1
            #     fine1 = fine_name
            #     seq1 = clade2seq1[fine1]['seq_id']
            #     clade_can0 = []
            #     for seq in seq1:
            #         # print(seq1)
            #         clade_can0.append(seq2clade0[seq])
            #     clade_can0 = list(set(clade_can0))
            #     if not len(clade_can0) == 1:
            #         print(f"fine1 {fine1} not included in one indel0 clade {clade_can0} count: {len(seq1)}, {len_can(clade2seq0,clade_can0)}")
            #         continue
            #     fine0 = clade_can0[0]
            #     seq0 = clade2seq0[fine0]['seq_id']
            #     assert len(seq0) >= count1, f"seq in indel0 {len(seq0)} not greater than seq in indel1 {count1}"
            #     clade_can1 = []
            #     for seq in seq0:
            #         clade_can1.append(seq2clade1[seq])
            #     clade_can1 = list(set(clade_can1))
            #     count_list = len_can(clade2seq1, clade_can1)
            #     len_can1 = sum(count_list)
            #     if len_can1 == len(seq0):
            #         clade_list = clade_can1
            #         count_list.insert(0,len(seq0))
            #         comp_Flag = True

            if len(clade_list) > 1 and comp_Flag:# cross not included
                assert len(clade_list)+1 == len(count_list), f"clade len +1 {len(clade_list)+1} not equal {len(count_list)}"
                comp_clades(fine0,clade_list,count_list,growthrate_list,line_list,mut_file,csvwriter1,csvwriter2)

def comp_clades(fine0,clades,counts,growthrates,linelist,mut_file,csvwriter1,csvwriter2):
    """

    :param clades: clade list to compare
    :param mut_file: mutation file path
    :return:
    """

    mutrans = torch.load(mut_file)

    clades2lineage = mutrans["clade_to_lineage"]
    clades_all = list(clades2lineage.keys())
    for clade in clades:
        assert clade in clades_all, f"clade error, {clade} not in clades dict"
    clade_id = mutrans['clade_id']

    clades_idx = [clade_id[clade] for clade in clades]
    # clades_idx = [clades_all.index(clade) for clade in clades] #old
    feat = mutrans["features"]
    mutrans = mutrans["mutations"]

    # feat_clades = feat[clades_idx,:]
    # res_patrn = "clade pair {} \nfeat diff {}"

    clade2count = dict(zip(clades,counts[1:]))
    clade2growth = dict(zip(clades,growthrates[1:]))
    clade2line = dict(zip(clades,linelist[1:]))
    # with open(out_file,'w') as c_file:
    #     csv_writer = csv.writer(c_file)
    #     header = ['clade0','clade1','clade2','count0','count1','count2','mut1','mut2']
    #     csv_writer.writerow(header)
    if len(clades) == 1:#if clade is the same in indel0 and indel1
        row = [fine0, clades[0], " ", counts[0], counts[1], " ", linelist[0],linelist[1]," ",\
               " ", " ", growthrates[0], growthrates[1], " "]

        csvwriter1.writerow(row)
        csvwriter2.writerow(row)
    else:
        for i_clade, clade1 in enumerate(clades):
            clade1_indx = clade_id[clade1]
            feat1 = feat[clade1_indx,:]
            for clade2 in clades[i_clade+1:]:
            # for j_clade in range(i_clade+1,len(clades)):
                print_flag = False
                feat_diff = {'1': [], '2': []}
                feat_all = {'1': [], '2': []}

                # clade2 = clades[j_clade]
                clade2_indx = clade_id[clade2]
                feat2 = feat[clade2_indx,:]
                # feat1 = feat_clades[i_clade,:]
                # feat2 = feat_clades[j_clade,:]

                for i_feat, mut in enumerate(mutrans):
                    if not feat1[i_feat] == 0:
                        feat_all['1'].append(mut)
                    if not feat2[i_feat] == 0:
                        feat_all['2'].append(mut)

                    if feat1[i_feat] != feat2[i_feat]:
                        if 'Del' in mut or 'Ins' in mut:
                            print_flag = True
                        if feat1[i_feat] != 0:
                            feat_diff['1'].append(mut)
                        elif feat2[i_feat] != 0:
                            feat_diff['2'].append(mut)
                        # feat_diff.append(mut)
                clade_pair = [clade1, clade2]

                # if print_flag == True:
                row = [fine0, clade1, clade2, counts[0], clade2count[clade_pair[0]], clade2count[clade_pair[1]], \
                       linelist[0],clade2line[clade_pair[0]],clade2line[clade_pair[1]],\
                       feat_diff['1'], feat_diff['2'], growthrates[0], clade2growth[clade_pair[0]],clade2growth[clade_pair[1]],\
                       feat_all['1'],feat_all['2']]
                csvwriter2.writerow(row)
                if print_flag == True:
                    csvwriter1.writerow(row)


if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='compare clades')
    # parser.add_argument('--suffix', type=str, default="3000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss15810344")

    results_dir = "/home/yanhongliang/projects/results"
    #max_loss = 18558864  # for indel1 resv 4

    # parser.add_argument('--mut_dir', type=str, default="/home/yanhongliang/projects/results")
    # parser.add_argument('--clade2seq')
    # parser.add_argument('--clades', nargs='+',type=str, default=['fine.31.0...123.0..170.11.0..0.0.0.1.0.0..0..0.0..4.0.1.1.0.1...0.0.2..1...3...24.0..0...5..2336.0.268....20.','fine.31.0...123.0..170.11.0..0.0.0.1.0.0..0..0.0..4.0.1.1.0.1...0.0.2..1...3...24.0..0...5..2336.0.268....20'])
    # parser.add_argument('--clades', nargs='+', default=['fine.31.0...123.0..170.11.0..0.0.0.1.0.0..0..0.0..4.0.1.1.0.1...0.0.2..1...3...24.0..0...5..16','fine.31.0...123.0..170.11.0..0.0.0.1.0.0..0..0.0..4.0.1.1.0.1...0.0.2..1...3...24.0..0...5..16.125.0'])
    args = parser.parse_args()

    # if len(args.clades) <2:
    #     print("return None, clades number no greater than 2")
    #     exit()
    # if args.suffix == "3000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0":
    #     print(f"using default suffix {args.suffix}")
    # diff_file = os.path.join(results_dir,"10_diff_clade.tsv")
    suffixes = [["2000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0","2000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss2606117544"],\
              ["3000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0","3000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss15810344"],
              ["4000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0","4000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss4449996"]]
    for suffix_pair in suffixes:
        suffix0, suffix1 = suffix_pair
        # suffix0 = "2000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0"
        # suffix1 = "2000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss2606117544"
        # suffix0 = "3000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0"
        # suffix1 = "3000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss15810344"
        # suffix0 = "4000.dkind1.rsnp0.aa2muc2.indel0.resv10.indelweight0.maxloss0"
        # suffix1 = "4000.dkind1.rsnp0.aa2muc2.indel1.resv4.indelweight1.maxloss4449996"
        indel0map = os.path.join(results_dir,'clade2accesionid',f"clade2accesion.{suffix0}.pkl")
        indel1map = os.path.join(results_dir,'clade2accesionid',f"clade2accesion.{suffix1}.pkl")
        indel0summary = os.path.join(results_dir,'summary',f"singlelineage.{suffix0}.csv")
        # indel0summary = os.path.join(results_dir,'summary',f"singlelineage.growthrates.obal-median.{suffix0}.csv")
        # indel1summary = os.path.join(results_dir,'summary',f"singlelineage.growthrates.obal-median.{suffix1}.csv")
        indel1summary = os.path.join(results_dir,'summary',f"singlelineage.{suffix1}.csv")

        mut_file = os.path.join(results_dir,'mutrans',f"mutrans.{suffix1}_new.pt")
        if not os.path.exists(os.path.join(results_dir,'diff_clade')):
            os.makedirs(os.path.join(results_dir,'diff_clade'))
        out_file1 = os.path.join(results_dir, 'diff_clade',f'DEL_comp_clade{suffix0[:4]}.csv')
        out_file2 = os.path.join(results_dir, 'diff_clade',f'comp_clade{suffix0[:4]}.csv')
        if not (os.path.exists(indel0map)):
            print(f"file error: suffix {indel0map}")
            continue

        if not os.path.exists(indel1map):
            print(f"file error: suffix {indel1map}")
            continue

        if not os.path.exists(indel0summary):
            print(f"file error: suffix {indel0summary}")
            continue

        if not os.path.exists(indel1summary):
            print(f"file error: suffix {indel1summary}")
            continue

        if not os.path.exists(mut_file):
            print(f"file error: suffix {imut_file}")
            continue

        select_clades(indel0map, indel1map,indel0summary, indel1summary,mut_file,out_file1, out_file2)

#select clades
    #compare clades
    # comp_clades(args.clades, mut_file)



