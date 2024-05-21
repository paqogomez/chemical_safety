from flask import Flask, render_template, request, redirect, url_for
import os
import shutil
import re
import json
from difflib import SequenceMatcher
from natsort import natsorted
from datetime import date

from chemical_safety.chemical import chemical

# Define the paths
package_config_path = os.path.join(os.path.dirname(__file__), 'config.json')
user_config_dir = os.path.expanduser('~/.chemical_safety')
user_config_path = os.path.join(user_config_dir, 'config.json')

# Ensure the user config directory exists
os.makedirs(user_config_dir, exist_ok=True)

# Copy default config to user directory if it doesn't exist
if not os.path.exists(user_config_path):
    shutil.copyfile(package_config_path, user_config_path)

# Load configuration
with open(user_config_path) as config_file:
    CONFIG = json.load(config_file)

def create_app():
    app = Flask(__name__)

    @app.route('/')
    def home():
        return render_template('index.html')

    @app.route('/chemical_lookup')
    def chemical_lookup():
        return render_template('search_form.html', lookup_type='Chemical', page_title='Chemical Lookup')

    @app.route('/multi_chemical_lookup')
    def multi_chemical_lookup():
        return render_template('search_form.html', lookup_type='multi chemical', page_title='Multiple Chemical Lookup')

    @app.route('/experiment_lookup')
    def experiment_lookup():
        return render_template('search_form.html', lookup_type='experiment', page_title='Experiment Lookup')

    @app.route('/room_lookup')
    def room_lookup():
        # Placeholder for label printing functionality
        return render_template('room_lookup.html')

    @app.route('/course_lookup')
    def course_lookup():
        return render_template('search_form.html', lookup_type='course', page_title='Course Lookup')

    @app.route('/emergency_info')
    def emergency_info():
        # Placeholder for displaying emergency information
        return render_template('emergency_info.html')

    @app.route('/secondary_label')
    def secondary_label():
        return render_template('search_form.html', lookup_type='secondary', page_title='Secondary Container Builder')

    @app.route('/phs_list')
    def phs_list():
        # Placeholder for label printing functionality
        return render_template('phs_list.html')

    @app.route('/lookup', methods=['GET', 'POST'])
    def lookup():
        if request.method == 'POST':
            search_term = request.form['search_term']
            lookup_type = request.form.get('lookup_type', 'chemical')
            # Redirect to the appropriate lookup route with the search term as a query parameter
            if lookup_type == 'Chemical':
                result = chemical(search_term)
                return render_template('chemical_lookup.html', lookup_type=lookup_type, result=result)
            elif lookup_type == 'multi chemical':
                result = [chemical(c) for c in search_term.split(', ')]
                return render_template('multi_chemical_lookup.html', lookup_type=lookup_type, result=result, experiment_name = "Custom Chemical List")
            elif lookup_type == 'experiment':
                chemlist,experiment_name = get_experiment_chem_list(search_term)
                result = [chemical(c) for c in chemlist]
                return render_template('experiment_lookup.html', lookup_type=lookup_type, result=result, experiment_name = experiment_name)
            elif lookup_type == 'course':
                course_data, course_title = build_course_summary(search_term)
                return render_template('course_lookup.html', lookup_type=lookup_type, course_data=course_data, course_title = course_title)
            elif lookup_type == 'secondary':
                result = [chemical(c) for c in search_term.split(', ')]
                return render_template('secondary_label_builder.html', lookup_type=lookup_type, result=result)
            else:
                return "Lookup type not supported", 400
        else:
            # Render the search form with a default lookup type context
            return render_template('search_form.html', lookup_type='Chemical', page_title='Chemical Lookup')
        
    @app.route('/generate_label', methods=['POST'])
    def generate_label():    
        container_name = request.form.get('container_name')
        generator_name = request.form.get('generator_name')
        signal_word = request.form.get('signal_word')
        chemical_cids = request.form.get('chemical_cids').split(',')
        disposal_info_list = list(set(request.form.getlist('disposal[]')))
        haz_waste_list = request.form.getlist('hazwaste[]')
        PHS_list = request.form.getlist('PHS[]')
        PHS_type = request.form.getlist('phs_type[]')

        haz_set = set()
        for hw in haz_waste_list:
            hw_designations = hw.split(', ')
            for hwd in hw_designations:
                haz_set.add(hwd)
        haz_list = list(haz_set)
        if len(haz_list)>0:
            hazwaste_info_string = ', '.join(haz_list)
        else:
            hazwaste_info_string = None

        disposal_info_list = disposal_info_list
        if len(disposal_info_list)>0:
            disposal_info_string = ', '.join(disposal_info_list)
        else:
            disposal_info_string = None

        PHS = any(PHS_list)

        if len(PHS_type) > 0:
            PHS_type = list(set(PHS_type))
        else:
            PHS_type = []

        print(haz_waste_list)
        print(disposal_info_string)
        
        all_selected_pictograms = set()
        all_selected_statements = set()

        for cid in chemical_cids:
            selected_pictograms = request.form.getlist(f'pictograms_{cid}')
            selected_statements = request.form.getlist(f'hazard_statements_{cid}')
            for sp in selected_pictograms:
                all_selected_pictograms.add(sp)
            for ss in selected_statements:
                all_selected_statements.add(ss)

        all_selected_statements = list(all_selected_statements)
        danger_statements = [s.replace("Danger:","<strong>Danger:</strong>").split('(')[0] for s in all_selected_statements if "Danger:" in s]
        warning_statements = [s.replace("Warning:","<strong>Warning:</strong>").split('(')[0] for s in all_selected_statements if "Warning:" in s]
        hazard_statements = danger_statements + warning_statements

        label_dict = {
            'container_name' : container_name,
            'signal_word' : signal_word,
            'generator' : generator_name,
            'pictograms' : list(all_selected_pictograms),
            'hazard_statements' : hazard_statements,
            'date' : date.today().strftime("%B %d, %Y"),
            'disposal' : disposal_info_string,
            'hazwaste' : hazwaste_info_string,
            'PHS' : PHS,
            'PHS_types' : PHS_type
        }

        return render_template('secondary_label.html', data = label_dict)

    return app

def enumerate(sequence, start=0):
    return zip(range(start, len(sequence) + start), sequence)

def build_course_summary(search_term):
    course_list = get_course_list()
    best_course, _ = custom_match(search_term,course_list)[0]
    directory_path = os.path.join('static/courses', best_course.replace(' ', ''))
    exp_names = [f for f in list_experiments(best_course)]

    exp_summary = []

    for file in exp_names:
        with open(os.path.join(directory_path, file+'.txt'), 'r') as f:
            experiment_chemical_data = [chemical(line.strip()) for line in f.readlines()]

        disposal = set()
        phs = False
        ppe = False

        for chem in experiment_chemical_data:
            if chem.WSU_particularly_hazardous:
                phs=True
                ppe = True
            if 'P262' in chem.p_codes:
                ppe = True
            if chem.disposal_info:
                for di in chem.disposal_info:
                    disposal.add(di)
            if chem.hazardous_waste:
                disposal.add(chem.hazardous_waste_info) 
        exp_summary.append({'name' : file, 'PHS' : phs, 'disposal' : list(disposal), 'PPE' : ppe, 'chem_data': experiment_chemical_data})

    return exp_summary, best_course
        
def get_course_list():

    user_static_dir = CONFIG.get('user_courses_dir', 'None')

    # Check if user has set a valid directory
    if user_static_dir == "None" or not os.path.exists(os.path.expanduser(user_static_dir)):
        directory_path = 'static/courses'
    else:
        directory_path = user_static_dir
    
    course_list = []

    pattern = re.compile(r'^([A-Z]{4})(\d{4})$')

    if os.path.exists(directory_path):
        for entry in os.listdir(directory_path):
            if os.path.isdir(os.path.join(directory_path, entry)):
                match = pattern.match(entry)
                if match:
                    course_name = f"{match.group(1)} {match.group(2)}"
                    course_list.append(course_name)
    else:
        print(f"Directory not found: {directory_path}")

    return natsorted(course_list)

def list_experiments(course):
    
    course = course.replace(' ', '')

    user_static_dir = CONFIG.get('user_courses_dir', 'None')

    # Check if user has set a valid directory
    if user_static_dir == "None" or not os.path.exists(os.path.expanduser(user_static_dir)):
        directory_path = os.path.join('static/courses', course)
    else:
        directory_path = os.path.join(os.path.expanduser(user_static_dir), course)
    
    txt_files = []
    if os.path.exists(directory_path):
        for filename in os.listdir(directory_path):
            if filename.endswith('.txt'):
                txt_files.append(os.path.splitext(filename)[0])
    else:
        print(f"Directory not found: {directory_path}")
    return natsorted(txt_files)

def get_experiment_chem_list(search_term):
    course_list = get_course_list()
    best_course = None
    best_distance = float('inf')  # Use infinity as initial comparison value
    
    pattern = re.compile(r'([A-Z]{4})\s*(\d{4})')
    match = pattern.search(search_term)

    if match:
        best_course, _ = custom_match(f"{match.group(0)} {match.group(1)}", course_list)[0]
        search_term = search_term.replace(match.group(0), '')

    user_static_dir = CONFIG.get('user_courses_dir', 'None')

    # Check if user has set a valid directory
    if user_static_dir == "None" or not os.path.exists(os.path.expanduser(user_static_dir)):
        directory_path = 'static/courses'
    else:
        directory_path = os.path.expanduser(user_static_dir)
    
    if best_course:
        directory_path = os.path.join(directory_path, best_course.replace(' ', ''))
        txt_files = [os.path.splitext(f)[0] for f in os.listdir(directory_path) if f.endswith('.txt')]
    else:
        txt_files_dict = {}
        for course in course_list:
            course_directory_path = os.path.join(directory_path, course.replace(' ', ''))
            for f in os.listdir(course_directory_path): 
                if f.endswith('.txt'):
                    f_name = os.path.splitext(f)[0]
                    txt_files_dict[f_name] = course_directory_path
        txt_files = txt_files_dict.keys()

    best_experiment, _ = custom_match(search_term, txt_files)[0]

    if best_course is None:
        directory_path = txt_files_dict[best_experiment]
        with open(os.path.join(directory_path, best_experiment+".txt"), 'r') as f:
            experiment_data = [line.strip() for line in f.readlines()]
        return natsorted(experiment_data), best_experiment
    return [], f"No valid match for {search_term}"    # Return an empty list if no match is found
    

def custom_match(search_str, choices, weight_number=0.7, weight_text=0.3):
    """
    Custom matching function that places more weight on numerical parts of the strings.
    
    :param search_str: The string to search for.
    :param choices: A list of strings to search against.
    :param weight_number: The weight to place on matching numbers. Default is 0.7.
    :param weight_text: The weight to place on the rest of the text. Default is 0.3.
    :return: A list of tuples with the match and its score, sorted by score.
    """
    search_numbers = [int(num) for num in re.findall(r'\d+', search_str)]
    results = []

    for choice in choices:
        choice_numbers = [int(num) for num in re.findall(r'\d+', choice)]
        number_score = 1 if search_numbers == choice_numbers else 0
        text_score = SequenceMatcher(None, search_str, choice).ratio() 

        # Calculate final score with weighted sum of number and text similarities
        final_score = (weight_number * number_score) + (weight_text * text_score)
        results.append((choice, final_score))

    # Sort the results based on the score in descending order
    return sorted(results, key=lambda x: x[1], reverse=True)

def dashboard():
    app = create_app()
    app.run(debug=True)

if __name__ == '__main__':
    dashboard()