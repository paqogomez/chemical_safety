import requests
import urllib.parse
import json
import subprocess
import os
import csv
from datetime import datetime, timedelta
import pandas as pd
import re
import Levenshtein as lev
from io import BytesIO
import base64
from rdkit import Chem as rdkChem
from rdkit.Chem import Draw as rdkDraw
from scipy import stats
import numpy as np

from .molecule import molecule

"""
The Important Thing About My_Class

The important thing about My_Class is that it isn't garbage.

It is designed with care,
each method serving a purpose,
each attribute neatly in place.
It handles data with grace,
and errors never leave a trace,
because they are caught and managed.
It is efficient,
reliable,
and its code is clean and readable.

But the important thing about My_Class is that it isn't garbage.
"""
 
WSU_PCODES_LIST = {'P222','P223','P250','P262','P263','P411','P412'}

class chemical:
    def __init__(self,chemical_name,cid=None):
        if chemical_name is None and cid is not None:
            self.cid = cid
            self.name = None
        else:
            self.name = chemical_name
            self.cid = self._get_cid(chemical_name)

        self.LD50_oral = None
        self.LD50_dermal = None
        self.LC50 = None

        self._full_json = self._get_pubchem_data(self.cid)
        self.full_name = self._full_json["Record"]["RecordTitle"]
        self.name =  self.full_name if self.name is None else self.name

        self.name_difference = self._parse_name_difference()

        self.dp_molecule, self.SMILES, self.IUPAC_name = self._parse_molecular_info()

        self.signal_word , self.pictograms, self.hazard_statements, self.hazard_codes, self.p_statements, self.p_codes = self._parse_GHS()

        self.WSU_particularly_hazardous, self.WSU_carcinogen, self.WSU_reproductive_toxin, self.WSU_highly_acute_toxin, self.WSU_No_GHS, self.WSU_PHC_info = self._check_if_particularly_hazardous()

        self.flash_point, self.boiling_point, self.melting_point, self.density = self._parse_physical_properties() 
        self.flammability_class, self.flammability_class_info = self._parse_flammability_info()

        self.peroxide_class, self.peroxide_class_info = self._parse_peroxide_info()

        self.hazardous_waste, self.hazardous_waste_info = self._parse_hazardous_waste_info()

        self.disposal_info = self._parse_disposal_info()

    def _get_cid(self,compound_name, recursive=True):
        try:
            compound_name = str(compound_name).lower()
        except:
            raise ValueError(f'Compound name could not be made a string. dtype = {type(compound_name)}')
        
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cache_folder = os.path.join(script_dir, 'data')
        csv_path = os.path.join(cache_folder, 'cid.csv')

        if os.path.exists(csv_path) and os.path.getsize(csv_path) > 0:
            with open(csv_path, mode='r', newline='', encoding='utf-8') as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    if row['name'].lower() == compound_name.lower():
                        date_retrieved = datetime.strptime(row['date_retrieved'], '%Y-%m-%d')
                        if (datetime.now() - date_retrieved).days < 30:
                            # The CID is recent enough to use
                            return int(row['cid'])
        else:
            # We are initializing the cid.csv list. create it and write the header
            os.makedirs(os.path.dirname(csv_path), exist_ok=True)
            with open(csv_path, mode='w', newline='', encoding='utf-8') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=['name', 'CAS', 'cid', 'date_retrieved'])
                writer.writeheader()

        # If the compound is not in the CSV or the date is too old, make the API call
        encoded_compound_name = urllib.parse.quote(compound_name)
        
        search_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_compound_name}/cids/JSON'
        response = requests.get(search_url)

        if response.status_code == 200:
            data = response.json()
            cid = int(data['IdentifierList']['CID'][0])
            with open(csv_path, mode='a', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([compound_name, '', cid, datetime.now().strftime('%Y-%m-%d')])
            return cid
        elif recursive:
            # initial lookup failed. need to get a little fuzzy with the autocomplete API
            autocomplete_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/autocomplete/Compound/{encoded_compound_name}/json'
            response = requests.get(autocomplete_url)

            if response.status_code == 200:
                data = response.json()
                suggestions = data['dictionary_terms']['compound']
                suggestion = sorted(suggestions, key=lambda x: lev.distance(compound_name, x))[0]
                print(f'No results for "{compound_name}". Trying "{suggestion}"')

                try:
                    cid = self._get_cid(suggestion)
                except:
                    cid = None
                if cid:
                    with open(csv_path, mode='a', newline='', encoding='utf-8') as csvfile:
                        writer = csv.writer(csvfile)
                        writer.writerow([suggestion, '', cid, datetime.now().strftime('%Y-%m-%d')])
                    return cid
            else:
                response.raise_for_status()

        else:
            response.raise_for_status()
        
    def _get_pubchem_data(self,cid):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        cache_folder = os.path.join(script_dir, 'data')
        cache_file_path = os.path.join(cache_folder, f'{self.cid}.json')

        data = None

        if os.path.exists(cache_file_path) and self._file_is_recent(cache_file_path):
            #print(f"Using cached data for CID: {cid}")
            with open(cache_file_path, 'r') as file:
                data = json.load(file)
        else:
            search_url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON?SourceID=GHS+Classification'

            try:
                # Make the request to the PubChem API
                response = requests.get(search_url)
                response.raise_for_status()  # Raise an error for bad status codes
                data = response.json()

                with open(cache_file_path, 'w') as file:
                    json.dump(data, file, indent=4)
        
            except requests.HTTPError as http_err:
                return f'HTTP error occurred: {http_err}'
            except Exception as err:
                return f'Other error occurred: {err}'
        return data

    def _file_is_recent(self, file_path, max_age_days=30):
        file_mod_time = datetime.fromtimestamp(os.path.getmtime(file_path))
        if (datetime.now() - file_mod_time) < timedelta(days=max_age_days):
            return True
        else:
            return False
        

    def _parse_molecular_info(self):
        molec = None
        IUPAC_name_string = ''
        Identifiers = next((item for item in self._full_json['Record']['Section'] if item['TOCHeading'] == 'Names and Identifiers'), None)
        if Identifiers:
            MolecFormula = next((item for item in Identifiers['Section'] if item['TOCHeading'] == 'Molecular Formula'), None)
            if MolecFormula: 
                for i in MolecFormula['Information']:
                    mf_string = i['Value']['StringWithMarkup'][0]['String']
                    mf_string = mf_string.replace('.','')
                    mf_string = mf_string.replace('[','(')
                    mf_string = mf_string.replace(']',')')
                    if len(mf_string) > 0:
                        molec = molecule(mf_string)
                        break

            CDs = next((item for item in Identifiers['Section'] if item['TOCHeading'] == 'Computed Descriptors'), None)
            if CDs:
                SMILES_entry = next((item for item in CDs['Section'] if item['TOCHeading'] == 'Canonical SMILES'), None)
                if SMILES_entry: 
                    for i in SMILES_entry['Information']:
                        SMILES_string = i['Value']['StringWithMarkup'][0]['String']
                        if len(SMILES_string) > 0:
                            break

                IUPAC_name_entry = next((item for item in CDs['Section'] if item['TOCHeading'] == 'IUPAC Name'), None)
                if IUPAC_name_entry:
                    for i in IUPAC_name_entry['Information']:
                        IUPAC_name_string = i['Value']['StringWithMarkup'][0]['String']
                        if len(IUPAC_name_string) > 0:
                            break 
        
        return molec, SMILES_string, IUPAC_name_string

    def _parse_GHS(self):
        """
        go into json
        get the GHS entry
        return the parameters we want
        """
        
        #get the section
        GHS = None
        SafetyInfo = next((item for item in self._full_json['Record']['Section'] if item['TOCHeading'] == 'Safety and Hazards'), None)
        if SafetyInfo:
            HazardID = next((item for item in SafetyInfo['Section'] if item['TOCHeading'] == 'Hazards Identification'), None)
            if HazardID:
                GHS = next((item for item in HazardID['Section'] if item['TOCHeading'] == 'GHS Classification'), None)       
        
        if GHS is None:
            return "Danger" , [], ['Warning: No GHS Record Available'], ['GHS_404'], [], []


        reference_options = []

        for item in GHS['Information']:
            if item['Name'] == 'Pictogram(s)':
                reference_number = item['ReferenceNumber']
                reference_options.append(reference_number)

        reference_num = self._decide_reference_number_GHS(reference_options)

        #parse the section
        pictogram_list = []
        signal_word = None
        hazard_statements = []
        hazard_codes = []
        p_codes = []
        p_codes_WSU = []
        p_statements = []

        Pictograms = next((item for item in GHS['Information'] if item['Name'] == 'Pictogram(s)' and item['ReferenceNumber'] == reference_num), None)
        if Pictograms:
            Markup = Pictograms['Value']['StringWithMarkup'][0]['Markup']
            if Markup:
                for m in Markup:
                    pictogram_list.append(m['Extra'])

        Signal = next((item for item in GHS['Information'] if item['Name'] == 'Signal' and item['ReferenceNumber'] == reference_num), None)
        if Signal:
            signal_word = Signal['Value']['StringWithMarkup'][0]['String']

        # Get the hazard statements
        HazStatements = next((item for item in GHS['Information'] if item['Name'] == 'GHS Hazard Statements' and item['ReferenceNumber'] == reference_num), None)
        if HazStatements:
            Statements = HazStatements['Value']['StringWithMarkup']
            if Statements:
                # Regex to extract H-code, percentage, and statement
                regex = r"H(\d{3,4}[a-zA-Z]*) \((\d+\.?\d*)%\): ([^[]+)\[.*\]"
                for s in Statements:
                    statement = s['String']
                    # Use regex to find matches
                    match = re.match(regex, statement)
                    if match:
                        h_code = "H"+match.group(1)       # H-code like '302'
                        percentage = float(match.group(2))  # percentage like '29.97'
                        text_statement = match.group(3).strip()  # statement text
                        if "Danger" in statement:
                            statement_signal_word = "Danger"
                        else:
                            statement_signal_word = "Warning"

                        if percentage > 20:
                            hazard_codes.append(h_code)
                            hazard_statements.append(f'{statement_signal_word}: {text_statement} ({h_code})')
        danger_statements = [s for s in hazard_statements if "Danger:" in s]
        warning_statements = [s for s in hazard_statements if "Warning:" in s]
        hazard_statements = danger_statements + warning_statements

        # Get the p statements
        pStatements = next((item for item in GHS['Information'] if item['Name'] == 'Precautionary Statement Codes' and item['ReferenceNumber'] == reference_num), None)
        if pStatements:
            StatementString = pStatements['Value']['StringWithMarkup'][0]['String']
            if StatementString:
                clean_string = re.sub(r'\band\b', '', StatementString).replace(' ', '')
                p_codes = clean_string.split(',')

        for pc in p_codes:
            if pc in WSU_PCODES_LIST:
                p_codes_WSU.append(pc)

        if p_codes_WSU:
            p_statements_dict = {}
            
            dir_path = os.path.dirname(os.path.abspath(__file__))  # Get the directory of the current script
            csv_path = os.path.join(dir_path, 'lists', 'p_statements.csv')  # Path to the CSV file within the same directory

            with open(csv_path, mode='r', encoding='utf-8') as file:
                reader = csv.reader(file)
                next(reader, None)  # Skip the header if there is one
                for row in reader:
                    if len(row) >= 2:
                        p_statements_dict[row[0]] = row[1]  # Assume P code is in the first column, statement in the second

            for pc in p_codes_WSU:
                p_statements.append(f'{pc}: {p_statements_dict[pc]}')

        
        
        return signal_word , pictogram_list, hazard_statements, hazard_codes, p_statements, p_codes_WSU
            
    def _decide_reference_number_GHS(self, reference_options):
        
        if len(reference_options) == 1:
            return reference_options[0]

        references=[]

        RefList = self._full_json['Record']['Reference']
        for ref in RefList:
            if ref['ReferenceNumber'] in reference_options:
                references.append({'number' : ref['ReferenceNumber'], 'chem_name' : ref['Name'], 'source_name' : ref['SourceName']})
        
        # Filter to keep only those with "ECHA" in the source name
        echa_references = [ref for ref in references if "ECHA" in ref['source_name']]
        
        if not echa_references:
            # If no ECHA references, could return an error or the first available ref, or another criteria
            return references[0]['number'] if references else 0

        # If there's exactly one ECHA reference, return it
        if len(echa_references) == 1:
            return echa_references[0]['number']

        # If there are multiple ECHA references, find the one with the chemical name most similar to self.full_name
        best_match = max(echa_references, key=lambda ref: lev.ratio(ref['chem_name'], self.full_name))
        return best_match['number']



    def generate_GHS_html(self, launch_page = True):
        
        pictograms_html = ''.join(f'<img src="pictograms/{pict}.svg" alt="Pictogram">' for pict in self.pictograms)
       
        sorted_hazard_statements = sorted(self.hazard_statements)
        danger_statements = [statement for statement in sorted_hazard_statements if "[Danger" in statement]
        warning_statements = [statement for statement in sorted_hazard_statements if "[Warning" in statement]
        combined_statements = danger_statements + [''] + warning_statements#endure Danger first
        hazard_statements_html = ''.join(
            f'<li>{statement.replace("[Danger", "[<strong>Danger</strong>").replace("[Warning", "[<strong>Warning</strong>")}</li>'
            if statement else '<li style="list-style-type:none;">&nbsp;</li>'
            for statement in combined_statements
        )
                
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Safety Information</title>
    <style>
        body {{ font-family: Arial, sans-serif; }}
        .compound-name {{ font-size: 32px; font-weight: bold; margin-bottom: 20px; }}
        .pictograms {{ display: flex; }}
        img {{ max-width: 100px; max-height: 100px; margin-right: 10px; }}
        .signal-word {{ font-weight: bold; font-size: 24px; }}
        .hazard-statements {{ list-style-type: none; padding: 0; }}
    </style>
</head>
<body>
    <div class="compound-name">
        {self.name}
    </div>
    <div class="pictograms">
        {pictograms_html}
    </div>
    <div class="signal-word">
        Signal Word: {self.signal_word}
    </div>
    <ul class="hazard-statements">
        {hazard_statements_html}
    </ul>
</body>
</html>"""
        sFile = 'safety_information.html'
        with open(sFile, 'w') as file:
            file.write(html_content)
        if launch_page:
            subprocess.run(['start', sFile], shell=True)
        return sFile
    
    def generate_imgstring(self):
        m = rdkChem.MolFromSmiles(self.SMILES)
        img = rdkDraw.MolToImage(m,size=(150,150),fitImage=False)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        return img_str
        

    #--------------------------------------------------------------------------------------------------------------------
    #Now we transition to some Weber-State specific checks that ultimately feed into a Particularly Hazardous designation
    def _check_cid_carcinogenlists(self):  
        OSHA_known = False
        NTP_known = False 
        NTP_anticipated = False 
        IARC = False 
        IARC_classification = 'Not listed'
        lists = []

        current_dir = os.path.dirname(os.path.abspath(__file__))
        data_file_path = os.path.join(current_dir, 'lists', 'OSHA_carcinogens.csv')
        
        df_OSHA = pd.read_csv(data_file_path)
        if self.cid in df_OSHA['cid'].values:
            OSHA_known=True
            lists.append('OSHA known carcinogens')

        data_file_path = os.path.join(current_dir, 'lists', 'OSHA_carcinogens.csv')
        df_known = pd.read_csv(data_file_path)
        if self.cid in df_known['cid'].values:
            NTP_known = True
            lists.append('National Toxicology Program known carcinogens')
        
        data_file_path = os.path.join(current_dir, 'lists', 'NTP_anticipated_carcinogens.csv')
        df_anticipated = pd.read_csv(data_file_path)
        if self.cid in df_anticipated['cid'].values:
            NTP_anticipated = True
            lists.append('National Toxicology Program anticipated carcinogens')

        data_file_path = os.path.join(current_dir, 'lists', 'IARC_classification.csv')
        df_IARC = pd.read_csv(data_file_path)
        IARC_classification = df_IARC.loc[df_IARC['cid'] == self.cid, 'group']
        if not IARC_classification.empty:
            lists.append(f'IARC classification: Group {IARC_classification.iloc[0]}')
            IARC = True

        WSU_CHP_carcinogen = False
        if OSHA_known or NTP_known or IARC == '1':
            sCertainty = "Known carcinogen"
            WSU_CHP_carcinogen = True
        elif NTP_anticipated or IARC == '2A':
            sCertainty = "Probable carcinogen"
            WSU_CHP_carcinogen = True
        elif IARC == '2B':
            sCertainty = "Possible carcinogen"
            WSU_CHP_carcinogen = True
        else:
            sCertainty = "Not a known carcinogen"
        
        return WSU_CHP_carcinogen,sCertainty,lists
    
    def _check_hazardstatements_reproductivetoxicity(self):
        reproductive_toxin = False
        flagged_codes = []
        reproductive_toxin_type = None
        pattern = r"^H[\dA-Za-z]+"
        
        cat1_tox_codes = ["H360","H360F","H360D","H360Fd","H360Df"]
        cat2_tox_codes = ["H361","H361f","H361d","H361fd"]
        bf_tox_codes = ["H362"]

        for hc in self.hazard_codes:
            if hc in cat1_tox_codes:
                reproductive_toxin = True
                reproductive_toxin_type = "Known or presumed human reproductive toxicant (OSHA Category 1)"
                flagged_codes.append(hc)
            if hc in cat2_tox_codes:
                reproductive_toxin = True
                reproductive_toxin_type = "Suspected human reproductive toxicant (OSHA Category 2)"
                flagged_codes.append(hc)
            if hc in bf_tox_codes:
                reproductive_toxin = True
                reproductive_toxin_type = "Possible toxic effects through lactation"
                flagged_codes.append(hc)


        return reproductive_toxin, reproductive_toxin_type, flagged_codes


    def _check_toxicity(self):
        WSU_highly_acute_toxin = False

        LD50s_oral=[]
        LD50s_dermal=[]
        LC50s=[]
        flagged_exposure_paths=[]

        if 'H300' in self.hazard_codes:# or 'H304' in self.hazard_codes:
            WSU_highly_acute_toxin = True
            flagged_exposure_paths.append('oral toxicity')
        
        if 'H310' in self.hazard_codes:
            WSU_highly_acute_toxin = True
            flagged_exposure_paths.append('dermal toxicity')

        if 'H330' in self.hazard_codes:
            WSU_highly_acute_toxin = True
            flagged_exposure_paths.append('inhalation toxicity')

        ToxicityInfo = next((item for item in self._full_json['Record']['Section'] if item['TOCHeading'] == 'Toxicity'), None)
        if ToxicityInfo:
            ToxicologyInfo = next((item for item in ToxicityInfo['Section'] if item['TOCHeading'] == 'Toxicological Information'), None)
            if ToxicologyInfo:
                NHTVs = next((item for item in ToxicologyInfo['Section'] if item['TOCHeading'] == 'Non-Human Toxicity Values'), None)
                if NHTVs:
                    ToxValueInfo = NHTVs['Information']
                    for i in ToxValueInfo:
                        try:
                            info_string = i['Value']['StringWithMarkup'][0]['String']
                            if "LD50" in info_string and "oral" in info_string.lower():
                                value = self._re_parse_LD50_string(info_string)
                                if value is not None:
                                    LD50s_oral.append(value)
                            if "LD50" in info_string and "dermal" in info_string.lower():
                                value = self._re_parse_LD50_string(info_string)
                                if value is not None:
                                    LD50s_dermal.append(value)
                            if "LC50" in info_string and "ppm" in info_string.lower():
                                value = self._re_parse_lc50_string(info_string)
                                if value is not None:
                                    LC50s.append(value)
                        except:
                            pass

        if len(LD50s_oral) > 0:
            self.LD50_oral = min(LD50s_oral)
            if self.LD50_oral <= 50:
                flagged_exposure_paths.append("oral toxicity")
                WSU_highly_acute_toxin = True
        
        if len(LD50s_dermal) > 0:
            self.LD50_dermal = min(LD50s_dermal)
            if self.LD50_dermal <= 200:
                flagged_exposure_paths.append("dermal toxicity")
                WSU_highly_acute_toxin = True
        
        if len(LC50s) > 0:
            self.LC50 = min(LC50s)
            if self.LC50 <= 200: #ppm gasses
                flagged_exposure_paths.append("inhalation toxicity")
                WSU_highly_acute_toxin = True
        

        return WSU_highly_acute_toxin, flagged_exposure_paths
    
    def _re_parse_LD50_string(self,info_string):
        regex = re.compile(r'LD50\s+\w+\s+(oral|dermal)\s+(?P<value>\d+(\.\d+)?)\s*(?P<unit>mg/kg|g/kg|ml/kg|ug/kg)(\s*bw)?', re.IGNORECASE)

        matches = list(regex.finditer(info_string))
        if matches:
            for match in matches:
                value = float(match.group('value'))
                unit = match.group('unit').lower()  # Normalize the unit for consistency

                # Convert the value to mg/kg based on the unit
                if unit == 'g/kg':
                    return value * 1000
                elif unit == 'mg/kg':
                    return value
                elif unit == 'ml/kg':
                    return value * 1000
                elif unit == 'ug/kg':
                    return value / 1000
        else:
            #print(f"Unexpected LD50 unit in {info_string}. Needs to be added to chemical._re_parse_LD50_string")
            return None
            
    def _re_parse_lc50_string(self, info_string):
        # Updated regex to capture LC50 inhalation data with units of ppm or mg/m³ (cu m)
        regex = re.compile(
            r'LC50\s+(Rat|Mouse)\s+(male\s+)?inhalation\s+(\d+(\.\d+)?(\s*-\s*\d+(\.\d+)?)?)\s*(ppm|mg/cu m|mg/m³)/\s*(\d+\s*min|\d+\s*hr)',
            re.IGNORECASE)

        matches = regex.finditer(info_string)
        if len(list(matches))>0:
            for match in matches:
                concentration = match.group(3)
                if '-' in concentration:
                    concentration = concentration.split('-')[0].strip()
                concentration = float(concentration)
                unit = match.group(7).lower()

                if unit == 'ppm':
                    return concentration
        else:
            #print(f"Unexpected LC50 unit in {info_string}. Needs to be added to chemical._re_parse_lc50_string")
            return None

    
    def _check_if_particularly_hazardous(self):
        particularly_hazardous = False
        
        info_dict = {}

        WSU_carcinogen,cancer_certainty,carcinogen_lists = self._check_cid_carcinogenlists()

        if WSU_carcinogen:
            info_dict['cancer info'] = {"certainty":cancer_certainty, "listings": carcinogen_lists}

        WSU_reproductive_toxin, reproductive_toxin_certainty, _ = self._check_hazardstatements_reproductivetoxicity()
        
        if WSU_reproductive_toxin:
            info_dict['reproductive toxin info'] = reproductive_toxin_certainty
        
        WSU_highly_acute_toxin, exposure_paths = self._check_toxicity()

        if WSU_highly_acute_toxin:
            info_dict['acute toxin info'] = exposure_paths

        if WSU_carcinogen or WSU_reproductive_toxin or WSU_highly_acute_toxin:
            particularly_hazardous = True

        if 'GHS_404' in self.hazard_codes:
            WSU_No_GHS = True
            info_dict['info'] = 'No GHS Information Available'
        else:
            WSU_No_GHS = False

        return particularly_hazardous, WSU_carcinogen, WSU_reproductive_toxin, WSU_highly_acute_toxin, WSU_No_GHS, info_dict
    
    def _parse_physical_properties(self):
        boiling_point = None 
        melting_point = None
        flash_point = None
        density = None
        
        PhysicalProp = next((item for item in self._full_json['Record']['Section'] if item['TOCHeading'] == 'Chemical and Physical Properties'), None)
        if PhysicalProp:
            ExperimentalProp = next((item for item in PhysicalProp['Section'] if item['TOCHeading'] == 'Experimental Properties'), None)
            if ExperimentalProp:
                BP = next((item for item in ExperimentalProp['Section'] if item['TOCHeading'] == 'Boiling Point'), None)
                if BP:
                    bp_list = []
                    for bp_i in BP["Information"]:
                        try:
                            bp_string = bp_i['Value']['StringWithMarkup'][0]['String']
                            matches = re.findall(r'(-?\d+\.?\d*)\s*°\s*([CF])', bp_string)
                            for match in matches:
                                temp, unit = float(match[0]), match[1]
                                if unit == 'C':
                                    temp = (temp * 9/5) + 32
                                bp_list.append(temp)
                        except:
                            pass
                    if len(bp_list)>0:
                        boiling_point = sum(bp_list) / len(bp_list)

                FP = next((item for item in ExperimentalProp['Section'] if item['TOCHeading'] == 'Flash Point'), None)
                if FP:
                    fp_list = []
                    for fp_i in FP["Information"]:
                        try:
                            fp_string = fp_i['Value']['StringWithMarkup'][0]['String']
                            matches = re.findall(r'(-?\d+\.?\d*)\s*°\s*([CF])', fp_string)
                            for match in matches:
                                temp, unit = float(match[0]), match[1]
                                if unit == 'C':
                                    temp = (temp * 9/5) + 32
                                fp_list.append(temp)
                        except:
                            pass
                    if len(fp_list) > 0:
                        flash_point = sum(fp_list) / len(fp_list)

                MP = next((item for item in ExperimentalProp['Section'] if item['TOCHeading'] == 'Melting Point'), None)
                if MP:
                    mp_list = []
                    for mp_i in MP["Information"]:
                        try:
                            mp_string = mp_i['Value']['StringWithMarkup'][0]['String']
                            matches = re.findall(r'(-?\d+\.?\d*)\s*°\s*([CF])', mp_string)
                            for match in matches:
                                temp, unit = float(match[0]), match[1]
                                if unit == 'C':
                                    temp = (temp * 9/5) + 32
                                mp_list.append(temp)
                        except:
                            pass
                    if len(mp_list)>0:
                        melting_point = sum(mp_list) / len(mp_list)

                densities = next((item for item in ExperimentalProp['Section'] if item['TOCHeading'] == 'Density'), None)
                if densities:
                    d_list = []
                    for d_i in densities["Information"]:
                        try:
                            d_string = d_i['Value']['StringWithMarkup'][0]['String']
                            matches = re.findall(r'\b(\d+\.\d+)\b', d_string)
                            for match in matches:
                                d = float(match)
                                d_list.append(d)
                        except:
                            pass
                    if len(d_list)>0:
                        density = _fudged_stats(d_list)
                
        return flash_point, boiling_point, melting_point, density
    
    def _parse_flammability_info(self):

        flam_class = None 
        flam_info = None
        flash_point = self.flash_point 
        boiling_point = self.boiling_point

        if flash_point:
            if flash_point < 100:   #class I
                if flash_point < 73: #class IA and IB
                    if boiling_point is None or boiling_point < 100:
                        flam_class = "IA"
                    else:
                        flam_class = "IB"
                else:
                    flam_class = "IC"
                flam_info = "Can ignite at room temperature."
            elif flash_point < 140:
                flam_class = "II"
                flam_info = "Can ignite when heated above 100 °F / 37 °C."
            elif flash_point < 200:
                flam_class = "IIIA"
                flam_info = "Can ignite when heated above 140 °F / 60 °C."
            elif flash_point is not None:
                flam_class = "IIIB"
                flam_info = "Can ignite when heated above 200 °F / 93 °C."

        return flam_class, flam_info

    def _parse_peroxide_info(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        data_file_path = os.path.join(current_dir, 'lists', 'UMN_peroxide_formers.csv')
        
        df_PFC = pd.read_csv(data_file_path)
        # Check if 'cid' exists in the 'CID' column and retrieve the associated peroxide former class
        if self.cid in df_PFC['CID'].values:
            peroxide_former_class = df_PFC.loc[df_PFC['CID'] == self.cid, 'PFC Class'].iloc[0]
        else:
            peroxide_former_class = None 
        peroxide_formation_descriptor = peroxide_former_class 

        if peroxide_former_class == 'A':
            peroxide_formation_descriptor = 'Forms explosive levels of peroxides without concentration'
        elif peroxide_former_class == 'B':
            peroxide_formation_descriptor = 'Forms peroxides on concentration (distillation/evaporation)'
        return peroxide_former_class, peroxide_formation_descriptor
    
    def _parse_hazardous_waste_info(self):
        current_dir = os.path.dirname(os.path.abspath(__file__))
        D_list_path = os.path.join(current_dir, 'lists', 'EPA_D_List.csv')
        P_list_path = os.path.join(current_dir, 'lists', 'EPA_P_List.csv')
        U_list_path = os.path.join(current_dir, 'lists', 'EPA_U_List.csv')
        
        haz_waste = False
        haz_waste_info = None

        df_D = pd.read_csv(D_list_path)
        if self.cid in df_D['cid'].values:
            haz_waste = True
            haz_waste_info = 'D Listed Waste'
        else:
            if self.dp_molecule:
                for e in self.dp_molecule.element_dict:
                    if e in ['As','Ba','Cd','Cr','Pb','Hg','Se','Ag']:
                        haz_waste = True
                        haz_waste_info = 'D Listed Waste'
                        
        df_P = pd.read_csv(P_list_path)
        if self.cid in df_P['cid'].values:
            haz_waste = True
            if haz_waste_info is None:
                haz_waste_info = 'P Listed Waste'
            else:
                haz_waste_info += ', P Listed Waste'

        df_U = pd.read_csv(U_list_path)
        if self.cid in df_U['cid'].values:
            haz_waste = True
            if haz_waste_info is None:
                haz_waste_info = 'U Listed Waste'
            else:
                haz_waste_info += ', U Listed Waste'
        
        return haz_waste, haz_waste_info
    
    def _parse_disposal_info(self):
        disposal_info = [] 

        corrosive_triggers = [
            'Corrosive' in self.pictograms,
            'H290' in self.hazard_codes
        ]

        env_tox_triggers = [
            'Environment' in self.pictograms,
            'H400' in self.hazard_codes,
            'H410' in self.hazard_codes,
            'H411' in self.hazard_codes
        ]


        if any(corrosive_triggers):
            disposal_info.append('Corrosive')
        
        if any(env_tox_triggers):
                disposal_info.append('Environmentally Toxic')

        return disposal_info
    
    def _parse_name_difference(self):
        entered_name = self.name
        pubchem_name = self.full_name

        # Calculate the Levenshtein distance between the two names
        distance = lev.distance(entered_name.lower(), pubchem_name.lower())

        # Consider names different if the distance is more than a small allowance (e.g., greater than 1)
        different = distance > 1

        return different
    
def _fudged_stats(data):
    #takes a list of data and tries to guess what data is real and what is garbage
    # mostly used to take values of BP/MP/density and throw out the random numbers that were incorrectly included (e.g. years)
    
    data = np.array(data)
    
    # Step 1: Compute initial mean and mode
    initial_mean = np.mean(data)
    mode_value, mode_count = stats.mode(data)
    
    # Step 2: Identify and remove outliers using the interquartile range (IQR) method
    Q1 = np.percentile(data, 25)
    Q3 = np.percentile(data, 75)
    IQR = Q3 - Q1
    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR
    cleaned_data = data[(data >= lower_bound) & (data <= upper_bound)]
    
    # Step 3: Recompute mean and mode of cleaned data
    cleaned_mean = np.mean(cleaned_data) if cleaned_data.size > 0 else initial_mean
    cleaned_mode_result = stats.mode(cleaned_data, nan_policy='omit')
    
    if isinstance(cleaned_mode_result.mode, np.ndarray):
        cleaned_mode_value = cleaned_mode_result.mode[0] if cleaned_mode_result.mode.size > 0 else None
        cleaned_mode_count = cleaned_mode_result.count[0] if cleaned_mode_result.count.size > 0 else 0
    else:
        cleaned_mode_value = cleaned_mode_result.mode
        cleaned_mode_count = cleaned_mode_result.count
    
    # Step 4: Return the best value representing the central tendency
    if cleaned_mode_count > 1:  # If mode is repeated, it's a good candidate
        return cleaned_mode_value
    else:
        return cleaned_mean
        

if __name__ == "__main__":
    from molecule import molecule

    chem = chemical('bromine')
    print('')
    print(chem.full_name)
    print(chem.SMILES)
    # if chem.name_difference:
    #     print(chem.name)
    
    # if(chem.dp_molecule):
    #     print(chem.dp_molecule.formula)

    # for hs in chem.hazard_statements:
    #     print(hs)
    
    # if chem.WSU_particularly_hazardous:
    #     print('')
    #     print("PARTICULARLY HAZARDOUS SUBSTANCE")
    #     print("Carcinogen:",chem.WSU_carcinogen)
    #     print("Reproductive toxin:",chem.WSU_reproductive_toxin)
    #     print("Acute toxin:",chem.WSU_highly_acute_toxin, f'{chem.LD50_oral}/50 oral, {chem.LD50_dermal}/200 dermal, {chem.LC50}/200 inhalation')
        
    #     for k in chem.WSU_PHC_info.keys():
    #         print(k,chem.WSU_PHC_info[k])

    # if chem.hazardous_waste:
    #     print('')
    #     print("HAZARDOUS WASTE")
    #     print(chem.hazardous_waste_info)
    
    # if chem.p_statements:
    #     print('')
    #     for p in chem.p_statements:
    #         print(p)
    
    # print('')
    # print('Flam class', chem.flammability_class)
    # print('perox',chem.peroxide_class)

    