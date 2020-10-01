
from collections import defaultdict
import pandas as pd
import os

def make_custom_table(query_dict, table_fields, include_snp_table):

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

        if query.display_name != query.name:
            df_dict["Query ID"].append(query.display_name.replace("|","\|"))
        else:
            df_dict["Query ID"].append(query.query_id.replace("|","\|"))
        
        
        if query.in_db:
            df_dict["Sequence name in tree"].append(query.name)
        else:
            df_dict["Closest sequence in tree"].append(query.closest)
            if include_snp_table or "SNPs" in table_fields or "snps" in table_fields:
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








