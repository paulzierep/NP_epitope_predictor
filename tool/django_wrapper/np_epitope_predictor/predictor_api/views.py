
from django.shortcuts import render, redirect
from django.http import HttpResponse

import pickle

from .forms import SMILES_input

from np_epitope_predictor.settings import PREDICTOR_PATH

import sys
sys.path.insert(0, PREDICTOR_PATH)

from run_prediction_wrapper import prediction2django

import random

import os
#from epitope_prediction import NP_Epitope_Prediction, Molecule_Group_Classifier, Epitope_Predictor
# sys.modules['Molecule_Group_Classifier'] = Molecule_Group_Classifier

from np_epitope_predictor.settings import MEDIA_ROOT

def index(request):

    if not request.method == 'POST':

        form = SMILES_input()
        return render(request, 'index.html', {'form':form})
    else:
        form = SMILES_input(request.POST) #create a form object
        if form.is_valid():
            smiles = form.cleaned_data['smiles']
            results = prediction2django(smiles)

            #simple clean up
            #removes all files in media except for the newest 100
            #avoids server overload of files
            max_files = 100

            files = list(os.listdir(MEDIA_ROOT))
            files = [os.path.join(MEDIA_ROOT, f) for f in files] # add path to each file
            files.sort(key=lambda x: os.path.getmtime(x), reverse=True)

            if len(files) > max_files:
                for file in files[max_files:]:
                    print(file)
                    os.remove(file)

            #unique files for each upload
            file_id = random.randint(0, 100000000)

            output_path = os.path.join(MEDIA_ROOT, '{0}.pickle'.format(file_id))

            with open(output_path, 'wb') as handle:
                pickle.dump(results, handle)

            return redirect('results', file_id = file_id)
            # return render(request, 'results.html', {'results':results})


def results(request, file_id):

    output_path = os.path.join(MEDIA_ROOT, '{0}.pickle'.format(file_id))

    with open(output_path, 'rb') as handle:
        results = pickle.load(handle)

    return render(request, 'results.html', {'results':results})

def help(request):

    return render(request, 'help.html', {})