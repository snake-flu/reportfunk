
import csv

def prepping_civet_arguments(name_stem_input, tree_fields_input, graphic_dict_input, label_fields_input, date_fields_input, table_fields_input):

    tree_fields = prep_argument_list(tree_fields_input, "NONE")
    label_fields = prep_argument_list(label_fields_input, "NONE")
    date_fields = prep_argument_list(date_fields_input, "NONE")
    table_fields = prep_argument_list(table_fields_input, "NONE")

    if "/" in name_stem_input:
        name_stem = name_stem_input.split("/")[-1]
    else:
        name_stem = name_stem_input

    graphic_dict = {}
    splits = graphic_dict_input.split(",")
    for element in splits:
        key = element.split(":")[0]
        value = element.split(":")[1]
        graphic_dict[key] = value
            
    for key in graphic_dict.keys():
        if key not in tree_fields:
            tree_fields.append(key)
  
    return name_stem, tree_fields, graphic_dict, label_fields, date_fields, table_fields

def prepping_llama_arguments(colour_fields_input, label_fields_input, name_stem_input):
    
    label_fields = prep_argument_list(label_fields_input, "NONE")
    colour_fields = prep_argument_list(colour_fields_input, "NONE")

    if "/" in name_stem_input:
        name_stem = name_stem_input.split("/")[-1]
    else:
        name_stem = name_stem_input

    return label_fields, colour_fields, name_stem


def prep_argument_list(list_input, absence):
    lst = []
    if list_input != absence:
        options = list_input.split(",")
        for i in options: 
            lst.append(i)
    return lst



def analyse_inputs(inputs):

    tree_fields, label_fields, graphic_dict, node_summary_option, map_sequences, mapping_trait, map_inputs = inputs 

    print("Showing " + ",".join(tree_fields) + " on the tree.")
    print(",".join(list(graphic_dict.keys())) + " fields are displayed graphically using " + ",".join(list(graphic_dict.values())) + " colour schemes respectively.")
    
    if label_fields != "NONE":
        print("Labelling by " + ",".join(label_fields) + " on the tree.")
    
    print("Summarising nodes by " + node_summary_option)

    if map_sequences:
        map_args = map_inputs.split(",")
        if len(map_args) == 2:
            if mapping_trait:
                print("Mapping sequences using columns " + map_args[0] + " " + map_args[1] + " for x values and y values respectively, and colouring by " + mapping_trait)
            else:
                print("Mapping sequences using columns " + map_args[0] + " " + map_args[1] + " for x values and y values respectively.")
        else:
            if mapping_trait:
                print("Mapping sequences using columns " + map_args[0] + " for outer postcodes, and colouring by " + mapping_trait)
            else:
                print("Mapping sequences using columns " + map_args[0] + " for outer postocdes.")


def prepping_adm2_adm1_data(full_metadata):

    official_adm2_adm1 = {'BARNSLEY': 'England', 'BATH AND NORTH EAST SOMERSET': 'England', 'BEDFORDSHIRE': 'England', 'BIRMINGHAM': 'England', 'BLACKBURN WITH DARWEN': 'England', 'BLACKPOOL': 'England', 'BOLTON': 'England', 'BOURNEMOUTH': 'England', 'BRACKNELL FOREST': 'England', 'BRADFORD': 'England', 'BRIGHTON AND HOVE': 'England', 'BRISTOL': 'England', 'BUCKINGHAMSHIRE': 'England', 'BURY': 'England', 'CALDERDALE': 'England', 'CAMBRIDGESHIRE': 'England', 'CENTRAL BEDFORDSHIRE': 'England', 'CHESHIRE EAST': 'England', 'CHESHIRE WEST AND CHESTER': 'England', 'CORNWALL': 'England', 'COVENTRY': 'England', 'CUMBRIA': 'England', 'DARLINGTON': 'England', 'DERBY': 'England', 'DERBYSHIRE': 'England', 'DEVON': 'England', 'DONCASTER': 'England', 'DORSET': 'England', 'DUDLEY': 'England', 'DURHAM': 'England', 'EAST RIDING OF YORKSHIRE': 'England', 'EAST SUSSEX': 'England', 'ESSEX': 'England', 'GATESHEAD': 'England', 'GLOUCESTERSHIRE': 'England', 'GREATER LONDON': 'England', 'HALTON': 'England', 'HAMPSHIRE': 'England', 'HARTLEPOOL': 'England', 'HEREFORDSHIRE': 'England', 'HERTFORDSHIRE': 'England', 'ISLE OF WIGHT': 'England', 'ISLES OF SCILLY': 'England', 'KENT': 'England', 'KINGSTON UPON HULL': 'England', 'KIRKLEES': 'England', 'KNOWSLEY': 'England', 'LANCASHIRE': 'England', 'LEEDS': 'England', 'LEICESTER': 'England', 'LEICESTERSHIRE': 'England', 'LINCOLNSHIRE': 'England', 'LUTON': 'England', 'MANCHESTER': 'England', 'MEDWAY': 'England', 'MIDDLESBROUGH': 'England', 'MILTON KEYNES': 'England', 'NEWCASTLE UPON TYNE': 'England', 'NORFOLK': 'England', 'NORTH LINCOLNSHIRE': 'England', 'NORTH SOMERSET': 'England', 'NORTH TYNESIDE': 'England', 'NORTH YORKSHIRE': 'England', 'NORTHAMPTONSHIRE': 'England', 'NORTHUMBERLAND': 'England', 'NOTTINGHAM': 'England', 'NOTTINGHAMSHIRE': 'England', 'OLDHAM': 'England', 'OXFORDSHIRE': 'England', 'PETERBOROUGH': 'England', 'PLYMOUTH': 'England', 'POOLE': 'England', 'PORTSMOUTH': 'England', 'READING': 'England', 'REDCAR AND CLEVELAND': 'England', 'ROCHDALE': 'England', 'ROTHERHAM': 'England', 'RUTLAND': 'England', 'SAINT HELENS': 'England', 'SALFORD': 'England', 'SANDWELL': 'England', 'SEFTON': 'England', 'SHEFFIELD': 'England', 'SHROPSHIRE': 'England', 'SLOUGH': 'England', 'SOLIHULL': 'England', 'SOMERSET': 'England', 'SOUTH GLOUCESTERSHIRE': 'England', 'SOUTH TYNESIDE': 'England', 'SOUTHAMPTON': 'England', 'SOUTHEND-ON-SEA': 'England', 'STAFFORDSHIRE': 'England', 'STOCKPORT': 'England', 'STOCKTON-ON-TEES': 'England', 'STOKE-ON-TRENT': 'England', 'SUFFOLK': 'England', 'SUNDERLAND': 'England', 'SURREY': 'England', 'SWINDON': 'England', 'TAMESIDE': 'England', 'TELFORD AND WREKIN': 'England', 'THURROCK': 'England', 'TORBAY': 'England', 'TRAFFORD': 'England', 'WAKEFIELD': 'England', 'WALSALL': 'England', 'WARRINGTON': 'England', 'WARWICKSHIRE': 'England', 'WEST BERKSHIRE': 'England', 'WEST SUSSEX': 'England', 'WIGAN': 'England', 'WILTSHIRE': 'England', 'WINDSOR AND MAIDENHEAD': 'England', 'WIRRAL': 'England', 'WOKINGHAM': 'England', 'WOLVERHAMPTON': 'England', 'WORCESTERSHIRE': 'England', 'YORK': 'England', 'ANTRIM AND NEWTOWNABBEY': 'Northern Ireland', 'ARMAGH, BANBRIDGE AND CRAIGAVON': 'Northern Ireland', 'BELFAST': 'Northern Ireland', 'CAUSEWAY COAST AND GLENS': 'Northern Ireland', 'DERRY AND STRABANE': 'Northern Ireland', 'FERMANAGH AND OMAGH': 'Northern Ireland', 'LISBURN AND CASTLEREAGH': 'Northern Ireland', 'MID AND EAST ANTRIM': 'Northern Ireland', 'MID ULSTER': 'Northern Ireland', 'NEWRY, MOURNE AND DOWN': 'Northern Ireland', 'NORTH DOWN AND ARDS': 'Northern Ireland', 'ABERDEEN': 'Scotland', 'ABERDEENSHIRE': 'Scotland', 'ANGUS': 'Scotland', 'ARGYLL AND BUTE': 'Scotland', 'CLACKMANNANSHIRE': 'Scotland', 'DUMFRIES AND GALLOWAY': 'Scotland', 'DUNDEE': 'Scotland', 'EAST AYRSHIRE': 'Scotland', 'EAST DUNBARTONSHIRE': 'Scotland', 'EAST LOTHIAN': 'Scotland', 'EAST RENFREWSHIRE': 'Scotland', 'EDINBURGH': 'Scotland', 'EILEAN SIAR': 'Scotland', 'FALKIRK': 'Scotland', 'FIFE': 'Scotland', 'GLASGOW': 'Scotland', 'HIGHLAND': 'Scotland', 'INVERCLYDE': 'Scotland', 'MIDLOTHIAN': 'Scotland', 'MORAY': 'Scotland', 'NORTH AYRSHIRE': 'Scotland', 'NORTH LANARKSHIRE': 'Scotland', 'ORKNEY ISLANDS': 'Scotland', 'PERTHSHIRE AND KINROSS': 'Scotland', 'RENFREWSHIRE': 'Scotland', 'SCOTTISH BORDERS': 'Scotland', 'SHETLAND ISLANDS': 'Scotland', 'SOUTH AYRSHIRE': 'Scotland', 'SOUTH LANARKSHIRE': 'Scotland', 'STIRLING': 'Scotland', 'WEST DUNBARTONSHIRE': 'Scotland', 'WEST LOTHIAN': 'Scotland', 'ANGLESEY': 'Wales', 'BLAENAU GWENT': 'Wales', 'BRIDGEND': 'Wales', 'CAERPHILLY': 'Wales', 'CARDIFF': 'Wales', 'CARMARTHENSHIRE': 'Wales', 'CEREDIGION': 'Wales', 'CONWY': 'Wales', 'DENBIGHSHIRE': 'Wales', 'FLINTSHIRE': 'Wales', 'GWYNEDD': 'Wales', 'MERTHYR TYDFIL': 'Wales', 'MONMOUTHSHIRE': 'Wales', 'NEATH PORT TALBOT': 'Wales', 'NEWPORT': 'Wales', 'PEMBROKESHIRE': 'Wales', 'POWYS': 'Wales', 'RHONDDA, CYNON, TAFF': 'Wales', 'SWANSEA': 'Wales', 'TORFAEN': 'Wales', 'VALE OF GLAMORGAN': 'Wales', 'WREXHAM': 'Wales'}

    contract_dict = {"SCT":"Scotland", "WLS": "Wales", "ENG":"England", "NIR": "Northern_Ireland"}

    illegal_values = ["", "NOT FOUND", "NONE","OTHER", "WALES", "UNKNOWN", "UNKNOWN SOURCE"]

    test_set = set()

    adm2_adm1 = official_adm2_adm1.copy()

    with open(full_metadata) as f:
        r = csv.DictReader(f)
        in_data = [x for x in r]
        for seq in in_data:
            if seq["country"] == "UK":
                adm1 = contract_dict[seq["adm1"].split("-")[1]]
                adm2 = seq["adm2"]
            
                if adm2.upper() not in illegal_values and adm2 not in official_adm2_adm1:
                    adm2_adm1[adm2] = adm1

    return adm2_adm1

        