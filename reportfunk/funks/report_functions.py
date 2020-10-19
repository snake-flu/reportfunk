#!/usr/bin/env python3

import os
import argparse
import yaml
import csv

def make_title(config):

    title = config["title"]
    if "#" not in title:
        title = "## " + title

    config["title"] = title

def report_date(config):
    date = config["report_date"]
    if str(date).startswith("Investigation started: "):
        report_date = config["report_date"]
    else:
        report_date = f'Investigation started: {date}'
    
    config["report_date"] = report_date


def multi_line_free_text(free_text_arg,config):

    if free_text_arg in config:
        if config[free_text_arg]:
            start = config[free_text_arg].replace("'","")
            
            end = "print('''" + start + "''')"
                
            config[free_text_arg] = end

def authors(config):

    if "authors" in config:
        if config["authors"]:

            authors_str = config["authors"]

            config["authors"] = "**Authors**: " + authors_str


def outbreak_id(config):

    if "outbreak_id" in config:
        if config["outbreak_id"]:
            outbreak_id = config["outbreak_id"]

            config["outbreak_id"] = "**Outbreak ID**: " + outbreak_id

def free_text_args(config):
    
    free_text_args = ["description", "conclusions"]
    
    report_date(config)
    authors(config)
    outbreak_id(config)

    for i in free_text_args:
        multi_line_free_text(i, config)




    

    