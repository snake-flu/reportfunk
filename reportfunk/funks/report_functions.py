#!/usr/bin/env python3

import os
import argparse
import yaml
import csv

def make_title(config):

    if "title" in config:
        title = config["title"]
        if "#" not in title:
            title = "# " + title

    else:
        if config["sequencing_centre"] == "DEFAULT":
            title = " # Cluster investigation "
        else:
            seq_centre = config["sequencing_centre_file"].split(".")[0]
            title = " # Cluster investigation for sequences generated by " + seq_centre

    config["title"] = title

def report_date(config):

    if report_date in config:
        date = config["report_date"]
        config["report_date"] = f'This investigation was started on {date}'

    else:
        config["report_date"] = ""

def multi_line_free_text(free_text_arg,config):

    if free_text_arg in config:
        
        start = config[free_text_arg].replace("'","")
        
        end = "print('''" + start + "''')"
             
    else:
        end = ""

    config[free_text_arg] = end


def authors(config):

    if "authors" in config:

        authors_str = config["authors"]

        config["authors"] = "**Authors**: " + authors_str

    else:
        config["authors"] = ""

def outbreak_id(config):

    if "outbreak_id" in config:

        outbreak_id = config["outbreak_id"]

        config["outbreak_id"] = "**Outbreak ID**: " + outbreak_id

    else:
        config["outbreak_id"] = ""

def appendix(omit_appendix_arg, config):

    if "omit_appendix" in config:
        value = config["omit_appendix"]
    else:
        value = omit_appendix_arg

    config["omit_appendix"] = value    

def free_text_args(config):
    
    free_text_args = ["description", "conclusions"]
    
    report_date(config)
    authors(config)
    outbreak_id(config)

    for i in free_text_args:
        multi_line_free_text(i, config)


def bars(bar_arg, config):

    if "add_bars" in config:
        final = config["add_bars"]
    else:
        final = bar_arg

    config["add_bars"] = final



    

    