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
from epiweeks import Week,Year

from reportfunk.funks.class_definitions import taxon,lineage

def convert_date(date_string):
    try:
        bits = date_string.split("-")
        date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
        return date_dt
    except ValueError:
        print("The wrong date format was supplied. Please re-run with YYYY-MM-DD")
    

def parse_collapsed_nodes(collapsed_node_file):

    collapsed_node_dict = defaultdict(list)

    with open(collapsed_node_file) as f:
        next(f)
        for l in f:
            toks = l.strip("\n").split(",")
            name = toks[0]
            members = toks[2].replace("[","").replace("]","")
            
            mem_list = members.split(" ")
            
            collapsed_node_dict[name] = mem_list
    
    return collapsed_node_dict


def parse_tree_tips(tree_dir, collapsed_node_file):

    collapsed_node_dict = parse_collapsed_nodes(collapsed_node_file)

    present_in_tree = [] #for pulling out the correct sequences from the background metadata to make objects
    tip_to_tree = {} #for finding which subtree the queries are in
    tree_to_all_tip = defaultdict(list) #for summarising trees when they are too big
    inserted_node_dict = defaultdict(dict) #for node summaries
    protected_sequences = []

    for fn in os.listdir(tree_dir):
        tree_name = fn.split(".")[0]
        if tree_name not in tree_to_all_tip: #for summarising trees when they are too big - contains subtrees
            all_tips = []
        else:
            all_tips = tree_to_all_tip[tree_name]
        
        if fn.endswith("tree"):
            tree = bt.loadNewick(tree_dir + "/" + fn, absoluteTime=False)
            for k in tree.Objects: 
                if k.branchType == 'leaf' and "inserted" not in k.name and "subtree" not in k.name:
                    if "collapsed" not in k.name:
                        present_in_tree.append(k.name)
                        all_tips.append(k.name)
                        tip_to_tree[k.name] = tree_name
                        protected_sequences.append(k.name)
                    else:
                        in_collapsed = collapsed_node_dict[k.name]
                        present_in_tree.extend(in_collapsed)
                        all_tips.extend(in_collapsed)

                if k.branchType == 'leaf' and "subtree" in k.name:
                    all_tips.append(k.name)

        elif fn.endswith(".txt") and fn != "collapse_report.txt":
            node_dict = defaultdict(list)
            with open(tree_dir + "/" + fn) as f:
                next(f)
                for l in f:
                    list_of_tips = []
                    tip_string = l.strip("\n").split("\t")[1]
                    tip_list = tip_string.split(",")
                    node_name = l.strip("\n").split("\t")[0]
                    for tip in tip_list:
                        if "collapsed" not in tip:
                            present_in_tree.append(tip)
                            all_tips.append(tip)
                            list_of_tips.append(tip)
                        else:
                            in_collapsed = collapsed_node_dict[tip]
                            present_in_tree.extend(in_collapsed)
                            all_tips.extend(in_collapsed)
                            list_of_tips.extend(in_collapsed)

                    node_dict[node_name] = list_of_tips

            inserted_node_dict[tree_name] = node_dict

        tree_to_all_tip[tree_name] = all_tips

    return present_in_tree, tip_to_tree, tree_to_all_tip, inserted_node_dict, protected_sequences

def parse_filtered_metadata(metadata_file, tip_to_tree, label_fields, tree_fields, table_fields, database_date_column):
    
    query_dict = {}
    query_id_dict = {}

    closest_seqs = set()

    tree_to_tip = defaultdict(list)

    with open(metadata_file, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames   

    with open(metadata_file, "r", encoding="utf-8") as f:
        in_data = csv.DictReader(f)
        for sequence in in_data:
            
            country = sequence["country"]
            query_id = sequence['query_id']
            query_name = sequence['query']
            closest_name = sequence["closest"]
            
            sample_date = sequence[database_date_column] #this may need to be flexible if using a different background database

            closest_distance = sequence["SNPdistance"]
            snps = sequence['SNPs']

            if query_id not in query_id_dict: #it's in the fasta file and in the db, this should take the db
               
                new_taxon = taxon(query_name, country, label_fields, tree_fields, table_fields)

                new_taxon.query_id = query_id

                if query_name == closest_name: #if it's in database, get its sample date
                    new_taxon.in_db = True
                    new_taxon.sample_date = sample_date
                    new_taxon.epiweek = Week.fromdate(convert_date(sample_date))
                    new_taxon.closest = "NA"
                else:
                    new_taxon.closest = closest_name
                    new_taxon.closest_distance = closest_distance
                    new_taxon.snps = snps
                    closest_seqs.add(closest_name)
                    
                if query_name in tip_to_tree:
                    relevant_tree = tip_to_tree[query_name]
                else:
                    relevant_tree = "NA"
                new_taxon.tree = relevant_tree

                tree_to_tip[relevant_tree].append(new_taxon)
            
                query_dict[query_name] = new_taxon
                query_id_dict[query_id] = new_taxon
            
    return query_dict, query_id_dict, tree_to_tip, closest_seqs

def UK_adm1(query_name, input_value):
    
    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}
    cleaning = {"SCOTLAND":"Scotland", "WALES":"Wales", "ENGLAND":"England", "NORTHERN_IRELAND": "Northern_Ireland", "NORTHERN IRELAND": "Northern_Ireland"}

    adm1 = input_value

    if "UK" in input_value:
        adm1_prep = input_value.split("-")[1]
        adm1 = contract_dict[adm1_prep]
    else:
        if input_value.upper() in cleaning.keys():
            adm1 = cleaning[input_value.upper()]
        else:
            for k,v in contract_dict.items():
                if k in query_name or v in query_name: #if any part of any country name is in the query name it will pick it up assign it
                    adm1 = v

    return adm1

def parse_input_csv(input_csv, query_id_dict, input_column, display_name, sample_date_column, tree_fields, label_fields, table_fields, context_table_summary_field, date_fields=None, UK_adm2_dict=None, patient_id_col=None, reinfection=False): 
    
    full_query_count = 0
    new_query_dict = {}
    
    with open(input_csv, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())
    
    with open(input_csv, 'r') as f:
        in_data = csv.DictReader(f)
        #in_data = [r for r in reader]
        for sequence in in_data:
            full_query_count += 1
            name = sequence[input_column]

            if name in query_id_dict.keys():
                taxon = query_id_dict[name]

                if reinfection:
                    taxon.attribute_dict["patient"] = sequence[patient_id_col]
                    
                taxon.input_display_name = sequence[display_name]

                if reinfection:
                    taxon.attribute_dict["patient"] = sequence[patient_id_col]

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "" and sequence[field] != "NA":
                            date_dt = convert_date(sequence[field])
                            taxon.date_dict[field] = date_dt 

                if sample_date_column in col_names: #if it's not in the background database or there is no date in the background database but date is provided in the input query
                    if sequence[sample_date_column] != "":
                        taxon.sample_date = sequence[sample_date_column]
                        taxon.epiweek = Week.fromdate(convert_date(sequence[sample_date_column]))

                if context_table_summary_field and context_table_summary_field in col_names:
                    if sequence["context_table_summary_field"] != "":
                        taxon.attribute_dict["context_table_summary_field"] = sequence[context_table_summary_field]
                     
                for col in col_names: #Add other metadata fields provided
                    if col in table_fields:
                        if sequence[col] != "":
                            taxon.table_dict[col] = sequence[col]
                    
                    if col in label_fields:
                        if sequence[col] != "":
                            taxon.attribute_dict[col] = sequence[col]
                    
                    if col in tree_fields and col != input_column and col != "adm1":
                        if sequence[col] != "":
                            taxon.attribute_dict[col] = sequence[col]
                    
                    if taxon.country == "UK": 
                        if col == "adm1":
                            adm1 = UK_adm1(name, sequence[col])
                            taxon.attribute_dict["adm1"] = adm1

                        if col == "adm2":

                            adm2 = sequence["adm2"] 
                            if "|" in adm2:
                                adm2 = "|".join(sorted(adm2.split("|")))
                            
                            taxon.attribute_dict["adm2"] = adm2 

                            if "location" in col_names:
                               location_label = sequence["location"]
                            else:
                                location_label = adm2

                            taxon.attribute_dict["location_label"] = location_label
                            
                            if "adm1" not in col_names and "adm1" in tree_fields:
                                if sequence[col] in UK_adm2_dict.keys():
                                    adm1 = UK_adm2_dict[sequence[col]]
                                    taxon.attribute_dict["adm1"] = adm1               

                new_query_dict[taxon.name] = taxon

      
    return new_query_dict, full_query_count 

def parse_background_metadata(query_dict, label_fields, tree_fields, table_fields, background_metadata, present_in_tree, closest_sequences, node_summary_option, tip_to_tree, database_name_column, database_sample_date_column, protected_sequences,context_table_summary_field, date_fields, virus):

    full_tax_dict = query_dict.copy()

    with open(background_metadata, 'r') as f:
        reader = csv.DictReader(f)
        col_name_prep = next(reader)
        col_names = list(col_name_prep.keys())

    old_data = False
    with open(background_metadata, 'r') as f:
        in_data = csv.DictReader(f)
        for sequence in in_data:
            
            seq_name = sequence[database_name_column]
            date = sequence[database_sample_date_column] 
            country = sequence["country"]

            if "adm2_raw" not in col_names: ##for civet
                old_data = True

            if "adm2" in col_names:
                adm2 = sequence['adm2']
                if "|" in adm2:
                    adm2 = "|".join(sorted(adm2.split("|")))

                if "location" in col_names:
                    location_label = sequence["location"]
                else:
                    location_label = adm2

                adm2_present_in_background = True
            else:
                adm2 = ""
                location_label = ""
                adm2_present_in_background = False

            # if virus == "sars-cov-2":	
            #     uk_lineage = sequence["uk_lineage"]	
            #     global_lineage = sequence["lineage"]	
            #     phylotype = sequence["phylotype"]	

            if node_summary_option == "adm2":
                if country != "UK":
                    node_summary_trait = country 
                else:
                    node_summary_trait = sequence["adm2"] 
            else:
                node_summary_trait = sequence[node_summary_option]

            if (seq_name in present_in_tree or seq_name in closest_sequences) and seq_name not in query_dict.keys():
                
                # if virus == "sars-cov-2":	
                #     new_taxon = taxon(seq_name, country, label_fields, tree_fields, table_fields, global_lineage=global_lineage, uk_lineage=uk_lineage, phylotype=phylotype)	
                # else:	
                new_taxon = taxon(seq_name, country, label_fields, tree_fields, table_fields)

                if date == "":
                    date = "NA"
                
                new_taxon.sample_date = date
                new_taxon.node_summary = node_summary_trait
                new_taxon.epiweek = Week.fromdate(convert_date(date))

                if new_taxon.name in protected_sequences:
                    new_taxon.protected = True

                if seq_name in tip_to_tree.keys():
                    new_taxon.tree = tip_to_tree[seq_name]

                new_taxon.attribute_dict["adm2"] = adm2
                new_taxon.attribute_dict["location_label"] = location_label

                new_taxon.input_display_name = seq_name

                for field in label_fields:
                    if field in col_names:
                        if sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                            new_taxon.attribute_dict[field] = sequence[field]

                if context_table_summary_field and context_table_summary_field in col_names:
                    if sequence[context_table_summary_field] != "":
                        new_taxon.attribute_dict["context_table_summary_field"] = sequence[context_table_summary_field]

                for field in table_fields:
                    if field in col_names:
                        if sequence[field] != "NA" and sequence[field] != "":
                            new_taxon.table_dict[field] = sequence[field]

                full_tax_dict[seq_name] = new_taxon

            

            #There may be sequences not in COG tree but that are in the full metadata, so we want to pull out the additional information if it's not in the input csv
            if seq_name in query_dict.keys(): 
                tax_object = query_dict[seq_name]
                if tax_object.sample_date == "NA" and date != "" and date != "NA":
                    tax_object.sample_date = date
                    converted = convert_date(date)
                    tax_object.all_dates.append(converted)
                    tax_object.epiweek = Week.fromdate(converted)

                
                if "adm2" not in tax_object.attribute_dict.keys() and adm2 != "":
                    tax_object.attribute_dict["adm2"] = adm2
                if "location_label" not in tax_object.attribute_dict.keys() and location_label != "":
                    tax_object.attribute_dict["location_label"] = location_label

                if context_table_summary_field and context_table_summary_field in col_names:
                    if sequence[context_table_summary_field] != "" and tax_object.attribute_dict["context_table_summary_field"] == "NA":
                        tax_object.attribute_dict["context_table_summary_field"] = sequence[context_table_summary_field]

                for field in date_fields:
                    if field in reader.fieldnames:
                        if sequence[field] != "" and sequence[field] != "NA" and field not in tax_object.date_dict.keys():
                            date_dt = convert_date(sequence[field])
                            tax_object.date_dict[field] = date_dt 
                    

                for field in tree_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                            if field != "adm1":
                                tax_object.attribute_dict[field] = sequence[field]
                            else:
                                if country == "UK":
                                    adm1 = UK_adm1(tax_object.name,sequence[field])
                                else:
                                    adm1 = "Other"
                                tax_object.attribute_dict[field] = adm1

                for field in label_fields:
                    if field in col_names:
                        if tax_object.attribute_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                                tax_object.attribute_dict[field] = sequence[field]

                for field in table_fields:
                    if field in col_names:
                        if tax_object.table_dict[field] == "NA" and sequence[field] != "NA" and sequence[field] != "": #this means it's not in the input file
                                tax_object.table_dict[field] = sequence[field]


                # if virus == "sars-cov-2":
                #     tax_object.global_lineage = global_lineage
                #     tax_object.uk_lineage = uk_lineage
                #     tax_object.phylotype = phylotype


                full_tax_dict[seq_name] = tax_object
                    
    return full_tax_dict, adm2_present_in_background, old_data

def parse_all_metadata(treedir, collapsed_node_file, filtered_background_metadata, background_metadata_file, input_csv, input_column, database_column, database_sample_date_column, display_name, sample_date_column, label_fields, tree_fields, table_fields, node_summary_option, context_table_summary_field, date_fields=None, UK_adm2_adm1_dict=None, reinfection=False, patient_id_col=None, virus="sars-cov-2"):

    present_in_tree, tip_to_tree, tree_to_all_tip, inserted_node_dict, protected_sequences = parse_tree_tips(treedir, collapsed_node_file)
    
    #parse the metadata with just those queries found in cog
    query_dict, query_id_dict, tree_to_tip, closest_sequences = parse_filtered_metadata(filtered_background_metadata, tip_to_tree, label_fields, tree_fields, table_fields, database_sample_date_column) 

    #Any query information they have provided
    query_dict, full_query_count = parse_input_csv(input_csv, query_id_dict, input_column, display_name, sample_date_column, tree_fields, label_fields, table_fields, context_table_summary_field, date_fields=date_fields, UK_adm2_dict=UK_adm2_adm1_dict, patient_id_col=patient_id_col, reinfection=reinfection)
    
    #parse the full background metadata
    full_tax_dict, adm2_present_in_background, old_data = parse_background_metadata(query_dict, label_fields, tree_fields, table_fields, background_metadata_file, present_in_tree, closest_sequences, node_summary_option, tip_to_tree, database_column, database_sample_date_column, protected_sequences, context_table_summary_field, date_fields=date_fields, virus=virus)

    return full_tax_dict, query_dict, tree_to_tip, tree_to_all_tip, inserted_node_dict, adm2_present_in_background, full_query_count, old_data 

def investigate_QC_fails(QC_file, input_column):

    fail_dict = {}

    with open(QC_file) as f:
        in_data = csv.DictReader(f)
        for sequence in in_data:
            name = sequence[input_column]
            reason = sequence["reason_for_failure"]
            final_reason = reason
            write_reason = True
            if "seq_len" in reason:
                length = reason.split(":")[1]
                final_reason = "Sequence too short: only " + length + " bases."
            elif "N_content" in reason:
                n_content = reason.split(":")[1]
                final_reason = "Sequence has too many Ns: " + str(float(round(float(n_content)*100))) + "\% of bases"
            elif "not_in_query_csv" in reason:
                final_reason = "Sequence not given in -i/--input"
            elif "already_in_tree" in reason:
                write_reason=False
            
            if write_reason:
                fail_dict[name] = final_reason


    return fail_dict

def investigate_missing_sequences(missing_seq_file):

    missing_list = []

    with open(missing_seq_file) as f:
        for l in f:
            missing_list.append(l.strip("\n"))
    
    return missing_list


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
                







