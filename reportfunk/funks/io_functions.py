#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import datetime 
from datetime import date
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
        if k not in ["generate_config","tempdir","summary_dir","treedir",
                    "collapse_summary","filtered_background_metadata","outfile"]:
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
            
            check_file = os.path.join(cwd,input_arg)
            path_to_check_file = os.path.abspath(os.path.dirname(check_file))

            if os.path.isfile(path_to_check_file):
                sys.stderr.write(cyan(f"Error: -i,--input accepts either a `.csv` or `.yaml` file, or a comma-separated string of IDs\n"))
                sys.exit(-1)

            else:
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

            elif ending == "xls":
                sys.stderr.write(cyan(f"Error: it looks like you've provided an excel file as input.\nPlease don't do this.\n-i,--input accepts either a csv or yaml file, or a comma-separated string of IDs\n"))
                sys.exit(-1)
            else:
                sys.stderr.write(cyan(f"Error: -i,--input accepts either a csv or yaml file, or a comma-separated string of IDs\n"))
                sys.exit(-1)

    return query,configfile
    
def parse_yaml_file(configfile,config):
    with open(configfile,"r") as f:
        input_config = yaml.load(f, Loader=yaml.FullLoader)
        for key in input_config:
            snakecase_key = key.replace("-","_")
            config[snakecase_key] = input_config[key]
    
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

def check_background_for_queries(config):

    data_column = config["data_column"]
    input_column = config["input_column"]
    queries = []
    with open(config["query"],"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            queries.append(row[input_column])
    c = 0
    with open(config["background_metadata"], "r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row[data_column] in queries:
                c +=1
    if c == 0:
        sys.stderr.write(cyan(f'Error: no valid queries to process.\n') + f'\
0 queries from `{input_column}` column matched in background metadata to `{data_column}`.\n')
        sys.exit(-1) 


def check_query_for_input_column(config):

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

def get_cluster_snakefile(thisdir):
    snakefile = os.path.join(thisdir, 'scripts','cluster_civet.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}\n Check installation\n'))
        sys.exit(-1)
    return snakefile

def get_query_fasta(fasta_arg,cwd,config):

    if fasta_arg:
        fasta = os.path.join(cwd, fasta_arg)

    elif "fasta" in config:
        if config["fasta"]:
            expanded_path = os.path.expanduser(config["fasta"])
            fasta = os.path.join(config["path_to_query"], expanded_path) 
        else:
            fasta = ""

    else:
        fasta = ""

    if fasta:
        if not os.path.isfile(fasta):
            sys.stderr.write(cyan(f'Error: cannot find fasta query at {fasta}\n'))
            sys.exit(-1)
        else:
            print(green(f"Input fasta file:") + f" {fasta}")
    
    config["fasta"] = fasta 

def make_timestamped_outdir(cwd,outdir,config):

    output_prefix = config["output_prefix"]
    split_prefix = output_prefix.split("_")
    if split_prefix[-1].startswith("20"):
        output_prefix = '_'.join(split_prefix[:-1])
    config["output_prefix"] = output_prefix
    timestamp = str(datetime.now().isoformat(timespec='milliseconds')).replace(":","").replace(".","").replace("T","-")
    outdir = os.path.join(cwd, f"{output_prefix}_{timestamp}")
    rel_outdir = os.path.join(".",timestamp)

    return outdir, rel_outdir

def get_outdir(outdir_arg,output_prefix_arg,cwd,config):
    outdir = ''
    
    add_arg_to_config("output_prefix",output_prefix_arg, config)
    
    if outdir_arg:
        expanded_path = os.path.expanduser(outdir_arg)
        outdir = os.path.join(cwd,expanded_path)
        rel_outdir = os.path.relpath(outdir, cwd) 

    elif config["update"] or config["cluster"]:
        outdir, rel_outdir = make_timestamped_outdir(cwd,outdir,config)

    elif "outdir" in config:
        expanded_path = os.path.expanduser(config["outdir"])
        outdir = os.path.join(config["path_to_query"],expanded_path)
        rel_outdir = os.path.relpath(outdir, cwd) 

    else:
        outdir, rel_outdir = make_timestamped_outdir(cwd,outdir,config)
    
    today = date.today()
    d = today.strftime("%Y-%m-%d")
    output_prefix = config["output_prefix"]
    split_prefix = output_prefix.split("_")
    if split_prefix[-1].startswith("20"):
        output_prefix = '_'.join(split_prefix[:-1])
    config["output_prefix"] = f"{output_prefix}_{d}"

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    report_output = os.path.join(outdir, "report")
    if not os.path.exists(report_output):
        os.mkdir(report_output)
    figure_output = os.path.join(report_output, "figures")
    if not os.path.exists(figure_output):
        os.mkdir(figure_output)

    print(green(f"Output dir:") + f" {outdir}")
    config["outdir"] = outdir 
    config["rel_outdir"] = os.path.join(".",rel_outdir) 
        
def get_temp_dir(tempdir_arg,no_temp_arg, cwd,config):
    tempdir = ''
    outdir = config["outdir"]
    if no_temp_arg:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
        config["no_temp"] = no_temp_arg
    elif config["no_temp"]:
        print(green(f"--no-temp:") + f" All intermediate files will be written to {outdir}")
        tempdir = outdir
    elif tempdir_arg:
        expanded_path = os.path.expanduser(tempdir_arg)
        to_be_dir = os.path.join(cwd,expanded_path)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    elif "tempdir" in config:
        expanded_path = os.path.expanduser(config["tempdir"])
        to_be_dir = os.path.join(cwd,expanded_path)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name

    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name
    
    config["tempdir"] = tempdir 
    return tempdir


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


def qc_list_inputs(config_key, column_names,config):

    background_headers = config["background_metadata_header"]

    list_of_fields = []
    arg_list = []

    input_to_check = config[config_key]

    if type(input_to_check) != bool:

        if not type(input_to_check) is list:
            input_to_check = input_to_check.split(",")
        
        for field in input_to_check:
            field = field.replace(" ","")
            if field in column_names or field in background_headers:
                list_of_fields.append(field)
            elif config_key == "table_fields" and field.lower() == "tree":
                print(green("Tree has been provided as a field, but will always be included in the table anyway so there is no need to specify it\n"))
            else:
                sys.stderr.write(cyan(f"Error: '{field}' column not found in query metadata file or background metadata file for {config_key}\n"))
                sys.exit(-1)

        field_str = ",".join(input_to_check)

        print(green(f"{config_key} shown:") + f" {field_str}")

    else:
        field_str = input_to_check

    config[config_key] = field_str
    
    return field_str

def check_date_format(date_string,config_key=None,row_number=None,column_name=None):
    date_format = '%Y-%m-%d'
    check_date= ""
    if date_string != "" and date_string != "NA":
        try:
            check_date = datetime.strptime(date_string, date_format).date()
        except:
            if row_number and column_name:
                sys.stderr.write(cyan(f"Error: Metadata field `{date_string}` [at column: {column_name}, row: {row_number}] contains unaccepted date format\nPlease use format {date_format}, i.e. `YYYY-MM-DD`\n"))
            else:
                sys.stderr.write(cyan(f"Error: Input '{date_string}' is the wrong date format.\nPlease use format {date_format}, i.e. `YYYY-MM-DD`\n"))

            sys.exit(-1)
            
    
    return check_date

def check_date_columns(config, date_column_list):

    metadata_header = config["background_metadata_header"]
    query = config["query"]
    query_header = config["query_metadata_header"]
    metadata = config["background_metadata"]

    # with open(query) as f:
    #     reader = csv.DictReader(f)
    #     query_header = reader.fieldnames

    # with open(metadata, "r", encoding = "utf-8") as f:
    #     reader = csv.DictReader(f)
    #     metadata_header = reader.fieldnames
    
    with open(query) as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for col in date_column_list:
            if col in query_header:
                row_num = 0
                for line in data:
                    row_num +=1
                    check_date_format(line[col],row_num, col)

    with open(metadata, "r", encoding = "utf-8") as f:
        reader = csv.DictReader(f)
        data = [r for r in reader]
        for col in date_column_list:
            if col in metadata_header:
                row_num = 0
                for line in data:
                    row_num +=1
                    check_date_format(line[col],row_num, col)

def check_metadata_for_search_columns(config):

    data_column = config["data_column"]
    
    with open(config["background_metadata"],"r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        if data_column not in header:
            sys.stderr.write(cyan(f"{data_column} column not in metadata\n"))
            sys.exit(-1)

    config["background_metadata_header"] = header
    

def data_columns_to_config(args,config):
    ## input_column
    add_arg_to_config("input_column",args.input_column, config)

    ## data_column
    add_arg_to_config("data_column",args.data_column, config)

def qc_dict_inputs(config_key, value_check, config):

    background_metadata_headers = config["background_metadata_header"]
    column_names = config["query_metadata_header"]
    output = []

    input_to_check = config[config_key]

    default_value = "Paired"

    if type(input_to_check) == str: 
        sections = input_to_check.split(",") 
    else:
        sections = input_to_check
    
    for item in sections:
        if "=" in item:
            splits = item.split("=")
        elif ":" in item:
            splits = item.split(":")
        else:
            splits = [item]        
        key = splits[0].replace(" ","")
        if key in column_names or key in background_metadata_headers:
            if len(splits) == 1:
                output.append(key + ":" + default_value)
            else:
                value = splits[1]
                if value in value_check or value == default_value or value == "default":
                    output.append(key + ":" + value)
                else:
                    sys.stderr.write(cyan(f"Error: {value} not compatible\n"))
                    sys.stderr.write(cyan(f"Please use one of {value_check}\n"))
                    sys.exit(-1)
           
        else:
            sys.stderr.write(cyan(f"Error: {key} field not found in metadata file or background metadata file for {config_key}\n"))
            sys.exit(-1) 

    output = ",".join(output)

    return output


def check_label_and_tree_and_date_fields(config):

    metadata = config["background_metadata"]
    metadata_headers = config["background_metadata_header"]
    
    acceptable_colours = get_colours()
    queries = []
    
    labels = []

    graphics_output = []
    column_names = []

    input_column = config["input_column"]
    data_column= config["data_column"]

    if not config["display_name"]:
        display_name = input_column
    else:
        display_name = config["display_name"]

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

    tree_field_str = qc_list_inputs("tree_fields", column_names, config)
    labels_str = qc_list_inputs("label_fields", column_names, config)
    date_field_str = qc_list_inputs("date_fields", column_names,config)
    
    if date_field_str:
        check_date_columns(config, date_field_str.split(",")) 

    graphic_dict_output = qc_dict_inputs("colour_by", acceptable_colours, config)

    for i in graphic_dict_output.split(","):
        element = i.split(":")[0]
        tree_field_list = tree_field_str.split(",")
        if element not in tree_field_str.split(",") and element != "adm1":
            tree_field_list.append(element)
            tree_field_str = ",".join(tree_field_list)
            config["tree_fields"] = tree_field_str
            # sys.stderr.write(cyan(f"Error: Field {element} in graphic dictionary but not in tree fields. Please add to tree fields if you want to colour by it.\n"))
            # sys.exit(-1)

    print(green(f"Colouring by: ") + f"{graphic_dict_output}")
    config["colour_by"] = graphic_dict_output

    database_sample_date_column = config["database_sample_date_column"]
    if database_sample_date_column not in metadata_headers:
        sys.stderr.write(cyan(f"Error: Field {database_sample_date_column} not in background metadata for sample date indication.\n"))
        sys.exit(-1)
    else:
        print(green(f'Using {database_sample_date_column} as sample date in background metadata'))

    #Removing this for now because it checks in the report maker if it's present. Easier I think because they don't have to provide any date here.
    # sample_date_column = config["sample_date_column"]
    # if len(column_names) > 1 and sample_date_column not in column_names and sample_date_column != "sample_date": #if the input is a query string, or they've provided an actual date column to check. 
    # #If it's default I don't think we check it because they don't to provide that data (as long as it's somewhere)
    #     sys.stderr.write(cyan(f"Error: Field {sample_date_column} not in query for sample date indication.\n"))
    #     sys.exit(-1)
    # else:
    #     print(green(f'Using {sample_date_column} as sample date in query metadata'))

def check_table_fields(table_fields, snp_data, config):
    
    column_names = config["query_metadata_header"]

    table_field_str = qc_list_inputs("table_fields", column_names, config)

    print(green(f"Displaying {table_field_str} as well as Query ID and Tree in table\n"))
    
    if config["include_snp_table"]:
        print(green(f"Showing SNP distance data in table\n"))
    else:
        print(green(f"Not showing SNP information in table\n"))

def check_summary_field(config_key, config):

    column_names = config["background_metadata_header"]

    summary_field = config[config_key]

    if summary_field not in column_names: 
        sys.stderr.write(cyan(f"Error: {summary_field} field not found in metadata file\n"))
        sys.exit(-1)
        
    print(green(f"Going to summarise collapsed nodes by: ") + f"{summary_field}")

def collapse_summary_path_to_config(config):
    path = os.path.join(config["outdir"],"catchment_trees", "tree_collapsed_nodes.csv")
    config["collapse_summary"] = path

def add_arg_to_config(key,arg,config):
    if arg:
        config[key] = arg

def input_file_qc(minlen_arg,maxambig_arg,config):
    post_qc_query = ""
    qc_fail = ""
    fasta = config["fasta"]

    add_arg_to_config("min_length",minlen_arg,config)
    add_arg_to_config("max_ambiguity",maxambig_arg,config)

    num_seqs =0

    if fasta != "":

        queries = []
        with open(config["query"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                queries.append(row[config["input_column"]])

        do_not_run = []
        passed = []
        
        

        for record in SeqIO.parse(fasta, "fasta"):
            if len(record) < config["min_length"]:
                record.description = record.description + f" fail=seq_len:{len(record)}"
                do_not_run.append(record)
                print(cyan(f"    - {record.id}\tsequence too short: Sequence length {len(record)}"))
            else:
                num_N = str(record.seq).upper().count("N")
                prop_N = round((num_N)/len(record.seq), 2)
                if prop_N > config["max_ambiguity"]: 
                    record.description = record.description + f" fail=N_content:{prop_N}"
                    do_not_run.append(record)
                    print(cyan(f"    - {record.id}\thas an N content of {prop_N}"))
                else:
                    if record.id not in queries:
                        record.description = record.description + f" fail=not_in_query_csv"
                        do_not_run.append(record)
                        print(cyan(f"    - {record.id}\tis not in query"))
                    else:
                        passed.append(record)
        
        passed_ids = [i.id for i in passed]
        already_in_tree = []
        with open(config["background_metadata"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row[config["data_column"]] in passed_ids:
                    already_in_tree.append(row[config["data_column"]])
        run = []

        for record in passed:
            if record.id in already_in_tree:
                record.description = record.description + f" fail=already_in_tree"
                do_not_run.append(record)
                print(cyan(f"    - {record.id}\tis already in tree"))
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

        num_seqs = len(run)

    config["post_qc_query"] = post_qc_query
    config["qc_fail"] = qc_fail
    config["num_seqs"] = num_seqs

    return num_seqs

def get_dict_of_metadata_filters(arg_type,to_parse, metadata):
    column_names =""
    query_dict = {}

    if not type(to_parse)==list:
        to_parse = to_parse.split(" ")
    with open(metadata, newline="",encoding = "utf-8") as f:
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
                sys.stderr.write(cyan(f"""Error: `{arg_type}` argument contains a column {column_name} that is not found in the metadata file supplied.
Columns that were found:\n{cols}"""))
                sys.exit(-1)
    return query_dict,column_names

def parse_date_range(metadata,column_name,to_search,rows_to_search):
    date_range = to_search.split(":")
    start_date = datetime.strptime(date_range[0], "%Y-%m-%d").date()
    end_date = datetime.strptime(date_range[1], "%Y-%m-%d").date()

    if rows_to_search == []:
        with open(metadata, newline="", encoding = "utf-8") as f:
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
        with open(metadata, newline="", encoding="utf-8") as f:
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

def filter_down_metadata(query_dict,metadata):
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

    return rows_to_search


def from_metadata_checks(config):
    if "query" in config:
        if config["query"]:
            if not config["update"]:
                sys.stderr.write(cyan('Error: please specifiy either -fm/--from-metadata or an input csv/ID string.\n'))
                sys.exit(-1)
    elif "fasta" in config:
        if config["fasta"]:
            sys.stderr.write(cyan('Error: fasta file option cannot be used in conjunction with -fm/--from-metadata.\nPlease specifiy an input csv with your fasta file.\n'))
            sys.exit(-1)

def generate_query_from_metadata(query,from_metadata, metadata, config):

    print(green("From metadata:"))
    to_parse = ""
    if from_metadata:
        to_parse = from_metadata

    elif "from_metadata" in config:
        to_parse = config["from_metadata"]

    data_column = config["data_column"]
    config["input_column"] = data_column
    
    # checks if field in metadata file and adds to dict: query_dict[country]=Ireland for eg
    query_dict,column_names = get_dict_of_metadata_filters("from_metadata",to_parse, metadata)
    
    rows_to_search = filter_down_metadata(query_dict,metadata)
    # if this is empty, for each column to search it'll open the whole file and search them
    # if it's not empty, it'll only search this list of rows
    
    

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
            sys.stderr.write(cyan(f"Error: No sequences meet the criteria defined with `--from-metadata`.\nPlease check your query is in the correct format (e.g. sample_date=YYYY-MM-DD).\nExiting\n"))
            sys.exit(-1)
        print(green(f"Number of sequences matching defined query:") + f" {count}")
        if len(query_ids) < 50:
            for i in query_ids:
                print(f" - {i}")
    return query


def parse_protect(protect_arg,metadata,config):

    print(green("Protect sequences:"))
    to_parse = ""
    if protect_arg:
        to_parse = protect_arg

    elif "protect" in config:
        to_parse = config["protect"]
    else:
        config["protect"] = False
        to_parse = False

    if to_parse:
        data_column = config["data_column"]
        
        query_dict,column_names = get_dict_of_metadata_filters("protect",to_parse, metadata)

        rows_to_search = filter_down_metadata(query_dict,metadata)

        protect = os.path.join(config["outdir"], "protected_background.csv")

        with open(protect,"w") as fw:
            writer = csv.DictWriter(fw, fieldnames=column_names,lineterminator='\n')
            writer.writeheader()
            count = 0

            protect_ids = []
            for row,c in rows_to_search:
                writer.writerow(row)
                count +=1
                protect_ids.append(row[data_column])

            if count == 0:
                print(cyan(f"Note: No sequences meet the criteria defined with `protect`.\n"))
                config["protect"] = False
            else:
                config["protect"] = protect
                print(green(f"Number of background sequences to be protected:") + f" {count}")


def collapse_config(collapse_threshold,config):

    add_arg_to_config("collapse_threshold",collapse_threshold, config)

    try:
        collapse_threshold = int(config["collapse_threshold"])
    except:
        sys.stderr.write(cyan(f"Error: collapse_threshold must be an integer\n"))
        sys.exit(-1)

    config["collapse_threshold"] = collapse_threshold

    print(green(f"Collapse threshold: ")+f"{collapse_threshold}")


def distance_config(distance,up_distance,down_distance,config):

    add_arg_to_config("distance",distance, config)
    add_arg_to_config("down_distance",down_distance, config)
    add_arg_to_config("up_distance",up_distance, config)

    try:
        distance = int(config["distance"])
        config["distance"] = distance
    except:
        sys.stderr.write(cyan(f"Error: distance must be an integer\n"))
        sys.exit(-1)

    if config["down_distance"]:
        try:
            config["down_distance"] = int(config["down_distance"])
        except:
            sys.stderr.write(cyan(f"Error: down_distance must be an integer\n"))
            sys.exit(-1)
    else:
        down_distance = config["distance"]
        config["down_distance"] = down_distance

    if config["up_distance"]:
        try:
            config["up_distance"] = int(config["up_distance"])
        except:
            sys.stderr.write(cyan(f"Error: up_distance must be an integer\n"))
            sys.exit(-1)
    else:
        up_distance = config["distance"]
        config["up_distance"] = up_distance


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
