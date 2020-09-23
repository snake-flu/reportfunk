
from collections import defaultdict
import pandas as pd
import os

def make_custom_table(query_dict, table_fields, snp_data_in_table):

    df_dict_indb = defaultdict(list)
    df_dict_seqprovided = defaultdict(list)
    indb = 0
    seqprovided = 0

    incogs = False
    seqprovideds = False

    for query in query_dict.values():
        if query.in_db:
            indb += 1
            df_dict = df_dict_indb
        else:
            df_dict = df_dict_seqprovided
            seqprovided += 1

        df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        if query.in_db:
            df_dict["Sequence name in tree"].append(query.name)
        else:
            df_dict["Closest sequence in tree"].append(query.closest)
            if snp_data_in_table or "SNPs" in table_fields:
                df_dict["Distance to closest sequence"].append(query.closest_distance)
                df_dict["SNPs"].append(query.snps)

        for field in table_fields:
            if field.lower() != "snps" or "distance":
                df_dict[field].append(query.table_dict[field])

        df_dict["Tree"].append(query.tree)
    
                
    if indb != 0:
        df_indb = pd.DataFrame(df_dict_indb)
        df_indb.set_index("Query ID", inplace=True)
        incogs = True
    
    if seqprovided != 0:
        df_seqprovided = pd.DataFrame(df_dict_seqprovided)
        df_seqprovided.set_index("Query ID", inplace=True)
        seqprovideds = True


    if seqprovideds and incogs:
        return df_indb, df_seqprovided
    elif seqprovideds and not incogs:
        return df_seqprovided
    elif incogs and not seqprovideds:
        return df_indb



def make_full_civet_table(query_dict, tree_fields, label_fields, input_column, outdir, table_fields, snp_data_in_table):

    df_dict_incog = defaultdict(list)
    df_dict_seqprovided = defaultdict(list)

    incog = 0
    seqprovided = 0
    incogs = False
    seqprovideds = False

    for query in query_dict.values():

        if query.in_db:
            df_dict = df_dict_incog
            incog += 1
        else:
            df_dict = df_dict_seqprovided
            seqprovided += 1
        
        df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        if query.in_db: 
            df_dict["Sequence name in Tree"].append(query.name)        

        df_dict["Sample date"].append(query.sample_date)

        if not query.in_db: 
            df_dict["Closest sequence in Tree"].append(query.closest)
            df_dict["Distance to closest sequence"].append(query.closest_distance)
            df_dict["SNPs"].append(query.snps)

        df_dict["UK lineage"].append(query.uk_lin)
        df_dict["Global lineage"].append(query.global_lin)
        df_dict["Phylotype"].append(query.phylotype)

        if query.tree != "NA":
            tree_number = query.tree.split("_")[-1]
            pretty_tree = "Tree " + str(tree_number)
            df_dict["Tree"].append(pretty_tree)
        else:
            df_dict["Tree"].append("NA") #this should never happen, it's more error catching

        if tree_fields != []:
            for i in tree_fields:
                df_dict[i].append(query.attribute_dict[i])
        
        if label_fields != []:
            for i in label_fields: 
                if i not in tree_fields and i != "sample_date" and i != input_column:
                    df_dict[i].append(query.attribute_dict[i])

    if incog != 0:
        df_incog = pd.DataFrame(df_dict_incog)
        file_name = os.path.join(outdir,"Sequences_already_in_cog")
        df_incog.to_csv(file_name, index=False)
        df_incog.set_index("Query ID", inplace=True)
        incogs = True
    
    if seqprovided != 0:
        df_seqprovided = pd.DataFrame(df_dict_seqprovided)
        file_name = os.path.join(outdir,"Sequences_provided")
        df_seqprovided.to_csv(file_name, index=False)
        df_seqprovided.set_index("Query ID", inplace=True)
        seqprovideds = True

    output = make_custom_table(query_dict, table_fields, snp_data_in_table)

    # print(output)
    # print(len(output))

    # for i in output:
    #     print(type(i))


    # for i in output:
    #     print("item")
    #     print(i)

    return output




