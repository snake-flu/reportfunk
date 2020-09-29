#!/usr/bin/env python3
from collections import defaultdict
from collections import Counter
import pandas as pd
import csv
from tabulate import tabulate
import reportfunk.funks.baltic as bt
import os
import datetime as dt
import math
import matplotlib.pyplot as plt

from reportfunk.funks.class_definitions import taxon,lineage

def convert_date(date_string):
    
    try:
        bits = date_string.split("-")
        date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
    except ValueError:
        print("The wrong date format was supplied. Please re-run with YYYY-MM-DD")
    
    return date_dt

def parse_tree_tips(tree_dir):

    tips = []
    tip_to_tree = {}

    for fn in os.listdir(tree_dir):
        if fn.endswith("tree"):
            tree_name = fn.split(".")[0]
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
            for k in tree.Objects:
                if k.branchType == 'leaf' and "inserted" not in k.name:
                    tips.append(k.name)
                    tip_to_tree[k.name] = tree_name

        elif fn.endswith(".txt"):
            with open(tree_dir + "/" + fn) as f:
                for l in f:
                    tip_string = l.strip("\n").split("\t")[1]
                    tip_list = tip_string.split(",")
                    tips.extend(tip_list)

    return tips, tip_to_tree

###MAYBE MAKE UK SPECIFIC METADATA FUNCTIONS, AND THEN GLOBAL ONES - or have UK=False or corona=False as defaults that we can set to true for civet?

def parse_filtered_metadata(metadata_file, tip_to_tree, label_fields, tree_fields, table_fields):
    
    query_dict = {}
    query_id_dict = {}

    tree_to_tip = defaultdict(list)

    uk_contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}


    with open(metadata_file, "r") as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            #can't assume adm1 is the country column? or maybe can in this function
            if "UK" in sequence["adm1"]:
                adm1_prep = sequence["adm1"].split("-")[1]
                adm1 = uk_contract_dict[adm1_prep]
            else:
                adm1 = sequence["adm1"]
            
            #all this specific stuff I think is only for the civet table, so maybe we make a new function for civet and then some general ones for date/tree/label fields?
            #maybe it's only the uk lineage and phylotype that's actually specific here to civet, and sample_date would need to genearlised somehow
            #glob_lin is specific to corona
            glob_lin = sequence['lineage']
            uk_lineage = sequence['uk_lineage']
            query_id = sequence['query_id']
            query_name = sequence['query']
            closest_name = sequence["closest"]
            closest_distance = sequence["closest_distance"]
            snps = sequence['snps']

                
            phylotype = sequence["phylotype"]
            sample_date = sequence["sample_date"]
            
            new_taxon = taxon(query_name, glob_lin, uk_lineage, phylotype, label_fields, tree_fields, table_fields)

            new_taxon.query_id = query_id

            #then this bit is only for civet - or just change it to "in_db"
            if query_name == closest_name: #if it's in COG, get it's sample date
                new_taxon.in_db = True
                new_taxon.sample_date = sample_date
                #new_taxon.attribute_dict["adm1"] = adm1
                new_taxon.closest = "NA"

            else:
                new_taxon.closest = closest_name
                new_taxon.closest_distance = closest_distance
                new_taxon.snps = snps
                
                if "adm1" in tree_fields:
                    for k,v in uk_contract_dict.items():
                        if k in query_name or v in query_name: #if any part of any country name is in the query name it will pick it up assign it
                            new_taxon.attribute_dict["adm1"] = v
            
            #this is civet specific - can't assume other programs are querying UK sequences
            # maybe add if country in col names here   
            new_taxon.country = "UK" 

            relevant_tree = tip_to_tree[query_name]
            new_taxon.tree = relevant_tree

            tree_to_tip[relevant_tree].append(new_taxon)
           
            query_dict[query_name] = new_taxon
            query_id_dict[query_id] = new_taxon
            
    return query_dict, query_id_dict, tree_to_tip

def Uk_adm1(input_value):
    
    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}
    cleaning = {"SCOTLAND":"Scotland", "WALES":"Wales", "ENGLAND":"England", "NORTHERN_IRELAND": "Northern_Ireland", "NORTHERN IRELAND": "Northern_Ireland"}

    if "UK" in input_value:
        adm1_prep = input_value.split("-")[1]
        adm1 = contract_dict[adm1_prep]
    else:
        if input_value.upper() in cleaning.keys():
            adm1 = cleaning[input_value.upper()]
        else:
            adm1 = input_value

    return adm1

def parse_input_csv(input_csv, query_id_dict, input_column, display_name, tree_fields, label_fields, adm2_adm1_dict, table_fields, date_fields=None): 
    
    new_query_dict = {}
    
    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            name = sequence[input_column]

            if name in query_id_dict.keys():
                taxon = query_id_dict[name]

                taxon.display_name = sequence[display_name]
                
                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "":
                            date_dt = convert_date(sequence[field])
                            taxon.date_dict[field] = date_dt 

                #keep this separate to above, because sample date is specifically needed
                if "sample_date" in col_names: #if it's not in database but date is provided (if it's in database, it will already have been assigned a sample date hopefully.)
                    if sequence["sample_date"] != "":
                        taxon.sample_date = sequence["sample_date"]
                elif "collection_date" in col_names:
                    if sequence["collection_date"] != "": #for the COG report 
                        taxon.sample_date = sequence["collection_date"]

                for col in col_names: #Add other metadata fields provided
                    if col in table_fields:
                        if sequence[col] != "":
                            taxon.table_dict[col] = sequence[col]
                    if col in label_fields:
                        if sequence[col] != "":
                            taxon.attribute_dict[col] = sequence[col]
                    else:
                        if col in tree_fields and col != "name" and col != "adm1":
                            if sequence[col] != "":
                                taxon.attribute_dict[col] = sequence[col]
                        
                        #uk specific stuff##########
                        if col == "adm1":
                            adm1 = Uk_adm1(sequence[col])
                            taxon.attribute_dict["adm1"] = adm1

                        if col == "adm2":
                            taxon.attribute_dict["adm2"] = sequence["adm2"]

                            if "adm1" not in col_names and "adm1" in tree_fields:
                                if sequence[col] in adm2_adm1_dict.keys():
                                    adm1 = adm2_adm1_dict[sequence[col]]
                                    taxon.attribute_dict["adm1"] = adm1

                        ########################

                new_query_dict[taxon.name] = taxon

      
    return new_query_dict 

def parse_full_metadata(query_dict, label_fields, tree_fields, table_fields, full_metadata, present_in_tree, node_summary_option, tip_to_tree, database_name_column, date_fields=None):

    full_tax_dict = query_dict.copy()

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    with open(full_metadata, 'r') as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            
            #these bits would need generalising to whatever database is the background to make it non-corona
            uk_lin = sequence["uk_lineage"]
            seq_name = sequence[database_name_column]

            date = sequence["sample_date"]
            adm2 = sequence["adm2"]
            country = sequence["country"]

            glob_lin = sequence["lineage"]
            phylotype = sequence["phylotype"]

            node_summary_trait = sequence[node_summary_option]

            if seq_name in present_in_tree and seq_name not in query_dict.keys():
                new_taxon = taxon(seq_name, glob_lin, uk_lin, phylotype, label_fields, tree_fields, table_fields)
                if date == "":
                    date = "NA"
                
                new_taxon.sample_date = date

                new_taxon.node_summary = node_summary_trait

                if seq_name in tip_to_tree.keys():
                    new_taxon.tree = tip_to_tree[seq_name]

                new_taxon.attribute_dict["adm2"] = adm2
                new_taxon.country = country 

                full_tax_dict[seq_name] = new_taxon

            #There may be sequences not in COG tree but that are in the full metadata, so we want to pull out the additional information if it's not in the input csv
            #Remember, then name has to match fully, so if it's not the x/y/z name this section won't work
            if seq_name in query_dict.keys(): 
                tax_object = query_dict[seq_name]
                if tax_object.sample_date == "NA" and date != "":
                    tax_object.sample_date = date
                    tax_object.all_dates.append(convert_date(date))
                
                if "adm2" not in tax_object.attribute_dict.keys() and adm2 != "":
                    tax_object.attribute_dict["adm2"] = adm2

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "" and field not in tax_object.date_dict.keys():
                            date_dt = convert_date(sequence[field])
                            tax_object.date_dict[field] = date_dt 
                    
                for field in label_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                                tax_object.attribute_dict[field] = sequence[field]

                for field in tree_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                            if field != "adm1":
                                tax_object.attribute_dict[field] = sequence[field]
                            else:
                                adm1 = Uk_adm1(sequence[field])
                                tax_object.attribute_dict[field] = adm1


                for field in table_fields:
                    if field in col_names:
                        if tax_object.table_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                                tax_object.table_dict[field] = sequence[field]

                    


                full_tax_dict[seq_name] = tax_object
                    
    return full_tax_dict
    

def parse_all_metadata(treedir, filtered_cog_metadata, full_metadata_file, input_csv, input_column, database_column, display_name, label_fields, tree_fields, table_fields, node_summary_option, adm2_to_adm1, date_fields=None):

    present_in_tree, tip_to_tree = parse_tree_tips(treedir)
    
    #parse the metadata with just those queries found in cog
    query_dict, query_id_dict, tree_to_tip = parse_filtered_metadata(filtered_cog_metadata, tip_to_tree, label_fields, tree_fields, table_fields) 

    #Any query information they have provided
    query_dict = parse_input_csv(input_csv, query_id_dict, input_column, display_name, tree_fields, label_fields, adm2_to_adm1, table_fields, date_fields)
    
    #parse the full background metadata
    full_tax_dict = parse_full_metadata(query_dict, label_fields, tree_fields, table_fields, full_metadata_file, present_in_tree, node_summary_option, tip_to_tree, database_column, date_fields)

    return full_tax_dict, query_dict, tree_to_tip    

def investigate_QC_fails(QC_file, tax_dict):

    fail_dict = {}

    with open(QC_file) as f:
        reader = csv.DictReader(f)
        in_data = [r for r in reader]
        for sequence in in_data:
            name = tax_dict[sequence["name"]].display_name
            reason = sequence["reason_for_failure"]

            if "seq_len" in reason:
                length = reason.split(":")[1]
                final_reason = "Sequence too short: only " + length + " bases."
            elif "N_content" in reason:
                n_content = reason.split(":")[1]
                final_reason = "Sequence has too many Ns: " + str(float(round(float(n_content)*100))) + "\% of bases"

            fail_dict[name] = final_reason


    return fail_dict


def find_new_introductions(query_dict, min_date): #will only be called for the COG sitrep, and the query dict will already be filtered to the most recent sequences

    lin_to_tax = defaultdict(list)
    intro_to_regions = defaultdict(dict)
    df_dict = defaultdict(list)
    new_intros = []
    lin_dict = {}

    df = "" #so that if there are no new ones then it returns an empty string

    for taxon in query_dict.values():
        lin = taxon.uk_lin
        lin_to_tax[lin].append(taxon)

    for k,v in lin_to_tax.items():
        new_lin = lineage(k,v)
        lin_dict[k] = new_lin    

    for key, value in lin_dict.items():
        if value.first_date > min_date:
            new_intros.append(value)


    for intro in new_intros:
        adm2s = []
        trees = set()
        for i in intro.taxa:
            if i.attribute_dict["adm2"] != "":
                adm2s.append(i.attribute_dict["adm2"])
                
            trees.add(i.tree)

        adm2_counts = Counter(adm2s)


        place_string = ""
        for place, count in adm2_counts.items():
            place_string += place + " (" + str(count) + ") "  
        

        df_dict["Name"].append(intro.name)
        df_dict["Size"].append(len(intro.taxa))
        df_dict["Locations"].append(place_string)
        df_dict["Global lineage"].append(intro.global_lins)
        df_dict["Trees"].append(trees)
        

    df = pd.DataFrame(df_dict)

    new_lins = []
    for i in df["Global lineage"]:
        new_lin = str(i).strip("{").strip("}").replace("'","")
        new_lins.append(new_lin)
    df["Global lineage"] = new_lins

    new_trees = []
    for i in df["Trees"]:
        new_tree = str(i).strip("{").strip("}").replace("'","")
        new_trees.append(new_tree)
    df["Trees"] = new_trees

    df.set_index("Name", inplace=True)

    return new_intros, df
                







