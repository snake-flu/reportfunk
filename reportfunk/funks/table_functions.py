
from collections import defaultdict
from collections import Counter
import pandas as pd
import os
import reportfunk.funks.parsing_functions as parse

def make_custom_table(query_dict, taxa_dict, table_fields, remove_snp_table):

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


        df_dict["Query ID"].append(query.display_name.replace("|","\|"))
        
        if query.in_db:
            df_dict["Sequence name in tree"].append(query.name)
        else:
            df_dict["Closest sequence in tree"].append(query.closest)
            if not remove_snp_table or "SNPs" in table_fields or "snps" in table_fields:
                df_dict["Distance to closest sequence"].append(query.closest_distance)
                df_dict["SNPs"].append(query.snps)

        for field in table_fields:
            if field.lower() != "snps" or "distance":
                if (field == "phylotype" or field == "uk_lineage" or field == "lineage") and not query.in_db:
                    df_dict[field].append(taxa_dict[query.closest].table_dict[field])
                else:
                    df_dict[field].append(query.table_dict[field])

        df_dict["Tree"].append(query.tree)
    
                
    if indb != 0:
        df_indb = pd.DataFrame(df_dict_indb)
        df_indb.set_index("Query ID", inplace=True)
        incogs = True
    
    if seqprovided != 0:
        df_seqprovided = pd.DataFrame(df_dict_seqprovided)
        if "phylotype" in df_seqprovided.columns:
            df_seqprovided.rename(columns={"phylotype": "phylotype of closest sequence"}, inplace=True)
        if "uk_lineage" in df_seqprovided.columns:
            df_seqprovided.rename(columns={"uk_lineage": "UK lineage of closest sequence"}, inplace=True) 
        if "lineage" in df_seqprovided.columns:
            df_seqprovided.rename(columns={"lineage": "lineage of closest sequence"}, inplace=True)
        df_seqprovided.set_index("Query ID", inplace=True)
        seqprovideds = True


    if seqprovideds and incogs:
        return df_indb, df_seqprovided
    elif seqprovideds and not incogs:
        return df_seqprovided
    elif incogs and not seqprovideds:
        return df_indb

def context_table(query_dict, taxa_dict, summarise_by):

    df_dict = defaultdict(list)
    values = []
    closest_values = {}
    no_values = []

    for query in query_dict.values():
        if not query.in_db and query.attribute_dict["context_table_summary_field"] == "NA" and taxa_dict[query.closest].attribute_dict["context_table_summary_field"]:
            to_add = taxa_dict[query.closest].attribute_dict["context_table_summary_field"]
            if to_add != "NA":
                closest_values[query.name] = to_add
            else:
                no_values.append(query.name)
        else:
            to_add = query.attribute_dict["context_table_summary_field"]
        
        if to_add != "NA":
            values.append(to_add)            
       
    value_count = Counter(values)

    summary = defaultdict(list)
    summary_dates = defaultdict(list)

    for value, count in value_count.items():
        for taxon in taxa_dict.values():
            if taxon.attribute_dict["context_table_summary_field"] == value and taxon.country == "UK":
                summary[value].append(taxon)
                summary_dates[value].append(parse.convert_date(taxon.sample_date))

    

    for summary, total_tax in summary.items():
        min_date = min(summary_dates[summary])
        pretty_min = min_date.strftime("%Y-%m-%d")
        max_date = max(summary_dates[summary])
        pretty_max = max_date.strftime("%Y-%m-%d")

        df_dict[summarise_by].append(summary)
        df_dict["Date range"].append(f"{pretty_min} to {pretty_max}")
        df_dict["Number in dataset"].append(value_count[summary])
        df_dict["Total in COG"].append(len(total_tax))

    df = pd.DataFrame(df_dict)
    df.set_index(summarise_by, inplace=True)

    return df, closest_values, no_values









