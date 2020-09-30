#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import datetime 
import tempfile
import pkg_resources
import yaml

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'

def make_config_file(config_name, config):
    config_to_write = {}
    for k in config:
        if k != "generate_config":
            config_to_write[k] = config[k]
    config_out = os.path.join(config["outdir"],config_name)
    with open(config_out,"w") as fw:
        yaml.dump(config_to_write, fw)
        print(green(f"Config file written to {config_out}."))
        sys.exit()



def type_input_file(input_arg,cwd,config):

    query,configfile="",""
    if input_arg:
        if "," in input_arg or "." not in input_arg:
            id_list = input_arg.split(",")
            print(green(f"ID string detected"))
            config["ids"] = id_list
        else:
            input_file = os.path.join(cwd,input_arg)
            path_to_file = os.path.abspath(os.path.dirname(input_file))
            config["path_to_query"] = path_to_file

            ending = input_file.split(".")[-1]

            if ending in ["yaml","yml"]:
                print(green(f"Input config file:") + f" {input_file}")
                configfile  = input_file

            elif ending == "csv":
                print(green(f"Input file:") + f" {input_file}")
                query = input_file

            else:
                sys.stderr.write(cyan(f"Error: -i,--input accepts either a csv or yaml file, or a comma-separated string of IDs"))
                sys.exit(-1)

    return query,configfile
    
def make_csv_from_ids(id_list, config):
    query = os.path.join(config["outdir"], "query.csv")
    with open(query,"w") as fw:
        in_col = "name"
        config["input_column"] = in_col
        fw.write(f"{in_col}\n")
        c = 0
        for i in id_list:
            c +=1
            fw.write(i+'\n')
        print(green(f"Number of IDs:") + f" {c}")
    return query

def parse_yaml_file(configfile,config):
    with open(configfile,"r") as f:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
        for key in input_config:
            snakecase_key = key.replace("-","_")

            if not snakecase_key in config:
                config[snakecase_key] = input_config[key]

    return config

def check_query_file(query, cwd, config):
    queryfile = ""
    
    if query:
        queryfile = query

    elif "query" in config:
        queryfile = os.path.join(config["path_to_query"],config["query"])
    elif "ids" in config:
        if type(config["ids"]) == list:
            queryfile = make_csv_from_ids(config["ids"], config)
        else:
            id_list = config["ids"].split(",")
            queryfile = make_csv_from_ids(id_list, config)

        config["path_to_query"] = config["outdir"]
    else:
        sys.stderr.write(cyan(f"Error: no query input provided\nPlease specify query or from_metadata\n"))
        sys.exit(-1)

    if os.path.exists(queryfile):
        config["query"] = queryfile
    else:
        sys.stderr.write(cyan(f"Error: cannot find query file at {queryfile}\nCheck if the file exists, or if you're inputting a set of ids in config (e.g. EPI12345,EPI23456) please provide them under keyword `ids`\n."))
        sys.exit(-1)

def check_query_for_input_column(config,default_dict):

    input_column = config["input_column"]
    
    with open(config["query"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if input_column not in header:
            sys.stderr.write(cyan(f"{input_column} column not in query metadata\n"))
            sys.exit(-1)

    config["query_metadata_header"] = header


def get_snakefile(thisdir):
    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n'))
        sys.exit(-1)
    return snakefile

def get_query_fasta(fasta_arg,cwd,config):

    if fasta_arg:
        fasta = os.path.join(cwd, fasta_arg)

    elif "fasta" in config:
        fasta = os.path.join(config["path_to_query"], config["fasta"]) 

    else:
        fasta = ""

    if fasta:
        if not os.path.exists(fasta):
            sys.stderr.write(cyan(f'Error: cannot find fasta query at {fasta}\n'))
            sys.exit(-1)
        else:
            print(green(f"Input fasta file:") + f" {fasta}")
    
    config["fasta"] = fasta 

def get_outdir(outdir_arg,cwd,config):
    outdir = ''
    
    if outdir_arg:
        expanded_path = os.path.expanduser(outdir_arg)
        outdir = os.path.join(cwd,expanded_path)
        rel_outdir = os.path.relpath(outdir, cwd) 

    elif "outdir" in config:
        expanded_path = os.path.expanduser(config["outdir"])
        outdir = os.path.join(config["path_to_query"],expanded_path)
        rel_outdir = os.path.relpath(outdir, cwd) 

    else:
        timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
        outdir = os.path.join(cwd, timestamp)
        
        rel_outdir = os.path.join(".",timestamp)
        
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    print(green(f"Output dir:") + f" {outdir}")
    config["outdir"] = outdir 
    config["rel_outdir"] = os.path.join(".",rel_outdir) 
        
def get_temp_dir(tempdir_arg,no_temp_arg, cwd,config):
    tempdir = ''
    outdir = config["outdir"]
    if no_temp_arg:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
    elif "no_temp" in config:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
    elif tempdir_arg:
        to_be_dir = os.path.join(cwd, tempdir_arg)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    elif "tempdir" in config:
        to_be_dir = os.path.join(cwd, config["tempdir"])
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name
    
    config["tempdir"] = tempdir 
    return tempdir

def local_lineages_to_config(central, neighbouring, region, config):

    if config["local_lineages"] == True:
        lineage_tables = []
        for r,d,f in os.walk(os.path.join(config["outdir"], 'figures')):
            for fn in f:
                if fn.endswith("_lineageTable.md"):
                    lineage_tables.append(os.path.join(config["outdir"], 'figures', fn))

        config["lineage_tables"] = lineage_tables
        config["lineage_maps"] = [central, neighbouring, region]
    else:
        config["lineage_tables"] = []
        config["lineage_maps"] = []


def get_tree_name_stem(tree_dir,config):
    tree_name_stems = []
    for r,d,f in os.walk(tree_dir):
        for fn in f:
            if fn.endswith(".tree"):
                basename = ".".join(fn.split(".")[:-1])
                stem = "_".join(basename.split("_")[:-1])
                tree_name_stems.append(stem)
    tree_name_stems = list(set(tree_name_stems))
    
    if len(tree_name_stems) > 1:
        sys.stderr.write("Error: Multiple tree names found")
        sys.exit(-1)
    elif len(tree_name_stems) == 0:
        sys.stderr.write("Error: No trees found in tree directory.")
        sys.exit(-1)
    else:
        tree_name_stem = tree_name_stems[0]

    config["tree_name_stem"] = tree_name_stem


def check_args_and_config_list(config_key, argument, column_names,config,default_dict):

    input_to_check = check_arg_config_default(config_key,argument,config,default_dict)

    background_headers = config["background_metadata_header"]

    list_of_fields = []
    arg_list = []

    if type(input_to_check) != bool:

        if not type(input_to_check) is list:
            input_to_check = input_to_check.split(",")
        
        for field in input_to_check:
            field = field.replace(" ","")
            if field in column_names or field in background_headers:
                list_of_fields.append(field)
            else:
                sys.stderr.write(cyan(f"Error: '{field}' column not found in query metadata file or background metadata file for {config_key}\n"))
                sys.exit(-1)

        field_str = ",".join(input_to_check)


        print(green(f"{config_key} shown:") + f" {field_str}")

    else:
        field_str = input_to_check

    config[config_key] = field_str
    
    return field_str

def check_date_format(date_string,row_number,column_name):
    date_format = '%Y-%m-%d'
    check_date= ""
    try:
        check_date = datetime.strptime(date_string, date_format).date()
    except:
        sys.stderr.write(cyan(f"Error: Metadata field `{date_string}` [at column: {column_name}, row: {row_number}] contains unaccepted date format\nPlease use format {date_format}, i.e. `YYYY-MM-DD`\n"))
        sys.exit(-1)
    
    return check_date

def check_date_columns(query, metadata, date_column_list):

    with open(query) as f:
        reader = csv.DictReader(f)
        query_header = reader.fieldnames

    with open(metadata) as f:
        reader = csv.DictReader(f)
        metadata_header = reader.fieldnames
    
    with open(query) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for col in date_column_list:
            if col in query_header:
                row_num = 0
                for line in data:
                    row_num +=1
                    check_date_format(line[col],row_num, col)

    with open(metadata) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for col in date_column_list:
            if col in metadata_header:
                row_num = 0
                for line in data:
                    row_num +=1
                    check_date_format(line[col],row_num, col)

def check_metadata_for_search_columns(config,default_dict):

    data_column = config["data_column"]
    
    with open(config["background_metadata"],"r") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if data_column not in header:
            sys.stderr.write(cyan(f"{data_column} column not in metadata\n"))
            sys.exit(-1)

    config["background_metadata_header"] = header
    config["data_column"] = data_column

def data_columns_to_config(args,config,default_dict):

    ## input_column
    input_column = check_arg_config_default("input_column",args.input_column, config, default_dict)
    config["input_column"] = input_column

    ## data_column
    data_column = check_arg_config_default("data_column",args.data_column, config, default_dict)
    config["data_column"] = data_column

def check_args_and_config_dict(argument, config_key, default_key, default_value,column_names, value_check, config):

    background_metadata_headers = config["background_metadata_header"]

    output = []

    if argument: 
        sections = argument.split(",")
        for item in sections:
            splits = item.split("=")
            key = splits[0].replace(" ","")
            if key in column_names or key in background_metadata_headers:
                if len(splits) == 1:
                    output.append(key + ":" + default_value)
                else:
                    value = splits[1]
                    if value in value_check or value == default_value:
                        output.append(key + ":" + value)
                    else:
                        sys.stderr.write(cyan(f"Error: {value} not compatible\n"))
                        sys.stderr.write(cyan(f"Please use one of {value_check}"))
                        sys.exit(-1)
            else:
                sys.stderr.write(cyan(f"Error: {key} field not found in metadata file or background metadata file for {config_key}\n"))
                sys.exit(-1)

    elif config_key in config:
        for key, value in config[config_key].items():
            key = key.replace(" ","")
            if value not in value_check and value != default_value:
                sys.stderr.write(cyan(f"Error: {value} not compatible\n"))
                sys.stderr.write(cyan(f"Please use one of {value_check}"))
                sys.exit(-1)
            elif key not in column_names and key not in background_metadata_headers:
                sys.stderr.write(cyan(f"Error: {key} field not found in metadata file or background metadata file for {config_key}\n"))
                sys.exit(-1)
            else:
                output.append(key + ":" + value)        

    else:
        output.append(default_key + ":" + default_value)

    output = ",".join(output)

    return output


def node_summary(node_summary,config):
    column_names = config["background_metadata_header"]

    if not node_summary and "node_summary" not in config:
        summary = "country"
    else:
        if "node_summary" in config:
            option = config["node_summary"]
        else:
            option = node_summary
        
        if option in column_names:
            summary = option
        else:
            sys.stderr.write(cyan(f"Error: {option} field not found in metadata file\n"))
            sys.exit(-1)
    
    print(green(f"Summarise collapsed nodes by:") + f" {summary}")
    config["node_summary"] = summary

def check_label_and_tree_and_date_fields(tree_fields, label_fields, colour_by_arg, date_fields, input_column, display_name_arg, config, default_dict,metadata):

    acceptable_colours = get_colours()
    queries = []
    
    labels = []

    graphics_output = []
    column_names = []

    input_column = config["input_column"]
    data_column= config["data_column"]

    display_name = check_arg_config_default("display_name",display_name_arg,config,default_dict) 
    if not display_name:
        display_name = input_column

    with open(config["query"], newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames

        if input_column not in column_names: # Checking input column present in query file
            sys.stderr.write(cyan(f"Error: Query file missing header field {input_column}\n"))
            sys.exit(-1)
        if display_name not in column_names:
            sys.stderr.write(cyan(f"Error: Query file missing header field {display_name}\n"))
            sys.exit(-1)
        else:
            config["input_column"] = input_column
            config["display_name"] = display_name

        print(green("Input querys to process:"))
        queries = []
        for row in reader:
            queries.append(row[input_column])
            
        print(green(f"Number of queries:") + f" {len(queries)}")

    tree_field_str = check_args_and_config_list("tree_fields", tree_fields, column_names, config, default_dict)
    
    labels_str = check_args_and_config_list("label_fields", label_fields, column_names, config, default_dict)

    date_field_str = check_args_and_config_list("date_fields",date_fields, column_names,config, default_dict)
    if date_field_str:
        check_date_columns(config["query"], metadata, date_field_str.split(",")) 
    
    graphic_dict_output = check_args_and_config_dict(colour_by_arg, "graphic_dict", default_dict["graphic_dict"], "default",column_names, acceptable_colours, config)

    for i in graphic_dict_output.split(","):
        element = i.split(":")[0]
        if element not in tree_field_str.split(",") and element != "adm1":
            sys.stderr.write(cyan(f"Error: Field {element} in graphic dictionary but not in tree fields. Please add to tree fields if you want to colour by it.\n"))
            sys.exit(-1)

    print(green(f"Colouring by: ") + f"{graphic_dict_output}")
    config["graphic_dict"] = graphic_dict_output

def check_table_fields(table_fields, snp_data, config, default_dict):
    
    with open(config["query"], newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames

    table_field_str = check_args_and_config_list("table_fields",table_fields, column_names, config, default_dict)

    if snp_data:
        config["snps_in_seq_table"] = True
    elif not snp_data and "snps_in_seq_table" not in config:
        config["snps_in_seq_table"] = False
    #otherwise it's just specified in the config

def check_summary_fields(summary_field, config):

    column_names = config["background_metadata_header"]

    if not summary_field:
        summary = "lineage"
    else:
        if summary_field in column_names:
            summary = summary_field
        else:
            sys.stderr.write(cyan(f"Error: {summary_field} field not found in metadata file\n"))
            sys.exit(-1)
        
    print(green(f"Going to summarise collapsed nodes by: ") + f"{summary}")
    config["node_summary"] = summary

def check_arg_config_default(key,arg,config,default):
    new_str = ""
    if arg:
        new_str = arg
    elif key in config:
        new_str = config[key]
    else:
        new_str = default[key]
    return new_str

def input_file_qc(minlen_arg,maxambig_arg,config,default_dict):
    post_qc_query = ""
    qc_fail = ""
    fasta = config["fasta"]

    minlen = check_arg_config_default("min_length",minlen_arg,config,default_dict)
    maxambig = check_arg_config_default("max_ambiguity",maxambig_arg,config,default_dict)

    config["min_length"] = minlen
    config["max_ambiguity"] = maxambig

    if fasta != "":
        do_not_run = []
        run = []
        for record in SeqIO.parse(fasta, "fasta"):
            if len(record) <minlen:
                record.description = record.description + f" fail=seq_len:{len(record)}"
                do_not_run.append(record)
                print(cyan(f"    - {record.id}\tsequence too short: Sequence length {len(record)}"))
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > maxambig: 
                    record.description = record.description + f" fail=N_content:{prop_N}"
                    do_not_run.append(record)
                    print(cyan(f"    - {record.id}\thas an N content of {prop_N}"))
                else:
                    run.append(record)

        post_qc_query = os.path.join(config["outdir"], 'query.post_qc.fasta')
        with open(post_qc_query,"w") as fw:
            SeqIO.write(run, fw, "fasta")
        qc_fail = os.path.join(config["outdir"],'query.failed_qc.csv')

        input_column = config["input_column"]
        with open(qc_fail,"w") as fw:
            fw.write(f"{input_column},reason_for_failure\n")
            for record in do_not_run:
                desc = record.description.split(" ")
                for i in desc:
                    if i.startswith("fail="):
                        fw.write(f"{record.id},{i}\n")

    config["post_qc_query"] = post_qc_query
    config["qc_fail"] = qc_fail

def get_dict_of_metadata_filters(to_parse, metadata):
    column_names =""
    query_dict = {}

    with open(metadata, newline="") as f:
        reader = csv.DictReader(f)
        column_names = reader.fieldnames
        
        # get each of the factors for the query
        for factor in to_parse:
            
            # eg country=Ireland 
            column_name,to_search = factor.split("=")
            print(f"\t- {column_name}\t {to_search}")
            # if the factor is in the metadata file add to the query dict
            if column_name in column_names:
                query_dict[column_name] = to_search

            else:
                # exit and print what the valid column names are
                cols = "\n- ".join(column_names)
                cols = cols + "\n"
                sys.stderr.write(cyan(f"""Error: `from-metadata` argument contains a column {column_name} that is not found in the metadata file supplied.
Columns that were found:\n{cols}"""))
                sys.exit(-1)
    return query_dict,column_names

def parse_date_range(metadata,column_name,to_search,rows_to_search):
    date_range = to_search.split(":")
    start_date = datetime.strptime(date_range[0], "%Y-%m-%d").date()
    end_date = datetime.strptime(date_range[1], "%Y-%m-%d").date()

    if rows_to_search == []:
        with open(metadata, newline="") as f:
            reader = csv.DictReader(f)
            c =0
            for row in reader:
                c +=1
                row_date = row[column_name]
                
                check_date = check_date_format(row_date,c,column_name)

                if start_date <= check_date <= end_date:
                    rows_to_search.append((row,c))
    else:
        last_rows_to_search = rows_to_search
        new_rows_to_search = []
        for row,c in last_rows_to_search:
            row_date = row[column_name]

            check_date = check_date_format(row_date,c,column_name)

            if start_date <= check_date <= end_date:
                new_rows_to_search.append((row,c))

        rows_to_search = new_rows_to_search
    return rows_to_search

def parse_general_field(metadata,column_name,to_search,rows_to_search):
    if rows_to_search == []:
        with open(metadata, newline="") as f:
            reader = csv.DictReader(f)
            c =0
            for row in reader:
                c +=1
                row_info = row[column_name]
                
                if row_info.upper() == to_search:
                    rows_to_search.append((row,c))
    else:
        last_rows_to_search = rows_to_search

        new_rows_to_search = []
        for row,c in last_rows_to_search:
            row_info = row[column_name]
            
            if row_info.upper() == to_search:

                new_rows_to_search.append((row,c))

        rows_to_search = new_rows_to_search

    return rows_to_search

def generate_query_from_metadata(from_metadata, metadata, config):

    print(green("From metadata:"))
    to_parse = ""
    if from_metadata:
        to_parse = from_metadata

    elif "from_metadata" in config:
        to_parse = config["from_metadata"]

    data_column = config["data_column"]
    config["input_column"] = data_column
    
    # checks if field in metadata file and adds to dict: query_dict[country]=Ireland for eg
    query_dict,column_names = get_dict_of_metadata_filters(to_parse, metadata)
    
    # if this is empty, for each column to search it'll open the whole file and search them
    # if it's not empty, it'll only search this list of rows
    rows_to_search = []
    
    for column_name in query_dict:
        
        to_search = query_dict[column_name].upper()

        # assumes its a date range if it has a ':' and startswith 2020-, 2019- or 2021-
        if ':' in to_search:
            if to_search.startswith("2020-") or to_search.startswith("2019-") or to_search.startswith("2021-"):
                print(f"Date range detected: {to_search}")
                rows_to_search = parse_date_range(metadata,column_name,to_search,rows_to_search)
        else:
            # parse by exact match 
            rows_to_search = parse_general_field(metadata,column_name,to_search,rows_to_search)

    query = os.path.join(config["outdir"], "from_metadata_query.csv")

    with open(query,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=column_names,lineterminator='\n')
        writer.writeheader()
        count = 0

        query_ids = []
        for row,c in rows_to_search:
            writer.writerow(row)
            count +=1
            query_ids.append(row[data_column])

        if count == 0:
            sys.stderr.write(cyan(f"Error: No sequences meet the criteria defined with `--from-metadata`.\nExiting\n"))
            sys.exit(-1)
        print(green(f"Number of sequences matching defined query:") + f" {count}")
        if len(query_ids) < 100:
            for i in query_ids:
                print(f" - {i}")
    return query

def distance_config(distance,up_distance,down_distance,config,default_dict):

    distance = check_arg_config_default("distance",distance, config, default_dict)
    down_distance = check_arg_config_default("down_distance",down_distance, config, default_dict)
    up_distance = check_arg_config_default("up_distance",up_distance, config, default_dict)

    try:
        distance = int(distance)
    except:
        sys.stderr.write(cyan(f"Error: distance must be an integer\n"))
        sys.exit(-1)

    if down_distance:
        try:
            config["down_distance"] = int(down_distance)
        except:
            sys.stderr.write(cyan(f"Error: down_distance must be an integer\n"))
            sys.exit(-1)
    else:
        config["down_distance"] = distance

    if up_distance:
        try:
            config["up_distance"] = int(up_distance)
        except:
            sys.stderr.write(cyan(f"Error: up_distance must be an integer\n"))
            sys.exit(-1)

    config["distance"] = distance

    print(green(f"Extraction radius:\n")+f"\tUp distance: {up_distance}\n\tDown distance: {down_distance}\n")

def get_colours():
    colours = ['viridis', 'plasma', 'inferno', 'magma', 'cividis','Greys', 
            'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper',
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
            'twilight', 'twilight_shifted', 'hsv',
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c',
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar']
    return colours


def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    elif 'cyan' in text_colour:
        coloured_text = 'cyan'
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING

def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING
