
import datetime as dt


def convert_date(date_string):
    bits = date_string.split("-")
    date_dt = dt.date(int(bits[0]),int(bits[1]), int(bits[2]))
    
    return date_dt

class taxon():
    #will need to change the args here to make data parsing general
    def __init__(self, name, country, label_fields, tree_fields, table_fields, global_lineage="NA", uk_lineage="NA",phylotype="NA"):
        self.name = name
        self.input_display_name = name
        # self.display_name = name

        self.sample_date = "NA"
        self.epiweek = "NA"

        self.date_dict = {}
        self.table_dict = {}       
        self.attribute_dict = {}

        self.protected = False

        self.country = country
       
        self.in_db = False
 
        for i in label_fields:
            self.attribute_dict[i.replace(" ","")] = "NA"
        for i in tree_fields:
            self.attribute_dict[i.replace(" ","")] = "NA"
        for i in table_fields:
            self.table_dict[i.replace(" ","")] = "NA"
        
        self.tree = "NA"

        self.closest_distance = "NA"
        self.snps = "NA"

        self.global_lineage = global_lineage
        self.uk_lineage = uk_lineage
        self.phylotype = phylotype

class lineage():
    
    def __init__(self, name, taxa):
        
        self.name = name
        self.taxa = taxa
        self.dates = []
        self.global_lins = set()
        
        for tax in taxa:
            if tax.sample_date != "NA":
                tax.date_dt = convert_date(tax.sample_date)
                self.dates.append(tax.date_dt)
            self.global_lins.add(tax.global_lin)
                
        if self.dates == []:
            self.first_date = "NA"
        else:
            self.first_date = min(self.dates)